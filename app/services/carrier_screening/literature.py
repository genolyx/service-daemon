"""
Variant Literature Search Module

PubMed/Europe PMC에서 변이 관련 문헌을 검색하고 SQLite에 영구 캐시합니다.

캐시 설계 원칙:
  - 논문(PMID)은 출판 후 변하지 않으므로 TTL 없이 영구 캐시
  - 검색 결과도 영구 캐시 (재검색은 명시적 요청 시에만)
  - 캐시 키: gene + hgvsc/hgvsp (변이 동일성 기준)

REST API 설계 원칙:
  - 모든 기능을 /api/literature/* 엔드포인트로 노출
  - Portal에 독립적으로 어떤 클라이언트에서도 호출 가능
  - portal/index.html은 이 API를 소비하는 얇은 클라이언트일 뿐

검색 전략 (PTM-platform pubmed.py 패턴 이식, 변이 검색에 맞게 조정):
  Tier 1: gene + hgvsc (예: "BRCA2 c.5266dupC")
  Tier 2: gene + hgvsp (예: "BRCA2 p.Gln1756fs")
  Tier 3: gene + effect 키워드 (예: "BRCA2 frameshift")
  Tier 4: gene alias fallback (MyGene.info)
  Fallback: Europe PMC
"""

import asyncio
import json
import logging
import os
import re
import sqlite3
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional, Tuple

from ...datetime_kst import now_kst_iso
from urllib.parse import quote_plus

logger = logging.getLogger(__name__)

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EUROPEPMC_BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"
MYGENE_BASE = "https://mygene.info/v3"

_EFFECT_KEYWORDS = {
    "frameshift_variant": "frameshift",
    "stop_gained": "nonsense",
    "splice_acceptor_variant": "splice acceptor",
    "splice_donor_variant": "splice donor",
    "missense_variant": "missense",
    "inframe_insertion": "inframe insertion",
    "inframe_deletion": "inframe deletion",
    "start_lost": "start codon",
    "stop_lost": "stop codon",
}


# ══════════════════════════════════════════════════════════════
# SQLite Cache
# ══════════════════════════════════════════════════════════════

_DB_PATH: Optional[str] = None
_DB_LOCK = asyncio.Lock()

# NCBI rate-limit: API key 허용 10 req/sec, 무key 3 req/sec
# 세마포어로 동시 요청 수를 제한하고, 429 시 지수 백오프로 재시도
_NCBI_SEMAPHORE = asyncio.Semaphore(5)
_NCBI_RETRY_DELAYS = (1.0, 2.0, 4.0)  # 최대 3회 재시도


def init_literature_db(db_path: str) -> None:
    """DB 초기화 (앱 시작 시 1회 호출)."""
    global _DB_PATH
    _DB_PATH = db_path
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.executescript("""
        CREATE TABLE IF NOT EXISTS search_cache (
            id           INTEGER PRIMARY KEY AUTOINCREMENT,
            cache_key    TEXT UNIQUE NOT NULL,
            gene         TEXT NOT NULL,
            hgvsc        TEXT,
            hgvsp        TEXT,
            query_terms  TEXT,
            pmids        TEXT NOT NULL,
            total_found  INTEGER DEFAULT 0,
            tiers_used   TEXT,
            cached_at    TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS article_cache (
            pmid         TEXT PRIMARY KEY,
            title        TEXT,
            abstract     TEXT,
            authors      TEXT,
            journal      TEXT,
            pub_date     TEXT,
            doi          TEXT,
            cached_at    TEXT NOT NULL
        );

        CREATE INDEX IF NOT EXISTS idx_search_gene ON search_cache(gene);
        CREATE INDEX IF NOT EXISTS idx_article_cached ON article_cache(cached_at);
    """)
    conn.commit()
    conn.close()
    logger.info(f"Literature DB initialized: {db_path}")


def _db() -> sqlite3.Connection:
    if not _DB_PATH:
        raise RuntimeError("Literature DB not initialized. Call init_literature_db() first.")
    conn = sqlite3.connect(_DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def _make_cache_key(gene: str, hgvsc: str = "", hgvsp: str = "") -> str:
    """캐시 키: gene + 변이 식별자 (hgvsc 우선, 없으면 hgvsp)."""
    variant_id = (hgvsc or hgvsp or "").strip()
    return f"gene:{gene.upper()}:variant:{variant_id}"


def _now_iso() -> str:
    return now_kst_iso()


# ── 캐시 읽기 ──────────────────────────────────────────────

def _get_search_cache(cache_key: str) -> Optional[Dict[str, Any]]:
    conn = _db()
    cur = conn.cursor()
    cur.execute("SELECT * FROM search_cache WHERE cache_key = ?", (cache_key,))
    row = cur.fetchone()
    conn.close()
    if not row:
        return None
    pmids = json.loads(row["pmids"]) if row["pmids"] else []
    return {
        "cache_key": row["cache_key"],
        "gene": row["gene"],
        "hgvsc": row["hgvsc"],
        "hgvsp": row["hgvsp"],
        "pmids": pmids,
        "total_found": row["total_found"],
        "tiers_used": json.loads(row["tiers_used"] or "{}"),
        "cached_at": row["cached_at"],
        "from_cache": True,
    }


def _get_articles_by_pmids(pmids: List[str]) -> List[Dict[str, Any]]:
    if not pmids:
        return []
    conn = _db()
    cur = conn.cursor()
    placeholders = ",".join("?" * len(pmids))
    cur.execute(f"SELECT * FROM article_cache WHERE pmid IN ({placeholders})", pmids)
    rows = cur.fetchall()
    conn.close()
    result = []
    for row in rows:
        result.append({
            "pmid": row["pmid"],
            "title": row["title"],
            "abstract": row["abstract"],
            "authors": json.loads(row["authors"] or "[]"),
            "journal": row["journal"],
            "pub_date": row["pub_date"],
            "doi": row["doi"],
            "cached_at": row["cached_at"],
        })
    return result


def _get_cached_pmids(pmids: List[str]) -> Tuple[List[Dict], List[str]]:
    """PMID 목록 중 캐시된 것과 아닌 것을 분리합니다."""
    cached = _get_articles_by_pmids(pmids)
    cached_pmids = {a["pmid"] for a in cached}
    uncached = [p for p in pmids if p not in cached_pmids]
    return cached, uncached


# ── 캐시 쓰기 ──────────────────────────────────────────────

def _save_search_cache(
    cache_key: str, gene: str, hgvsc: str, hgvsp: str,
    query_terms: str, pmids: List[str], tiers_used: Dict,
) -> None:
    conn = _db()
    cur = conn.cursor()
    cur.execute("""
        INSERT OR REPLACE INTO search_cache
            (cache_key, gene, hgvsc, hgvsp, query_terms, pmids, total_found, tiers_used, cached_at)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        cache_key, gene, hgvsc, hgvsp, query_terms,
        json.dumps(pmids), len(pmids),
        json.dumps(tiers_used), _now_iso(),
    ))
    conn.commit()
    conn.close()


def _save_articles(articles: List[Dict[str, Any]]) -> None:
    if not articles:
        return
    conn = _db()
    cur = conn.cursor()
    now = _now_iso()
    for art in articles:
        cur.execute("""
            INSERT OR IGNORE INTO article_cache
                (pmid, title, abstract, authors, journal, pub_date, doi, cached_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            art.get("pmid", ""),
            art.get("title", ""),
            art.get("abstract", ""),
            json.dumps(art.get("authors", [])),
            art.get("journal", ""),
            art.get("pub_date", ""),
            art.get("doi", ""),
            now,
        ))
    conn.commit()
    conn.close()


# ══════════════════════════════════════════════════════════════
# 메인 Public API
# ══════════════════════════════════════════════════════════════

async def search_variant_literature(
    gene: str,
    hgvsc: str = "",
    hgvsp: str = "",
    effect: str = "",
    max_results: int = 10,
    ncbi_email: str = "service-daemon@example.com",
    ncbi_api_key: str = "",
    ncbi_tool: str = "service-daemon",
    force_refresh: bool = False,
) -> Dict[str, Any]:
    """
    변이에 대한 PubMed 문헌을 검색하고 캐시합니다.

    Returns:
        {
            "gene": str,
            "hgvsc": str,
            "hgvsp": str,
            "pmids": list[str],
            "articles": list[dict],
            "total_found": int,
            "tiers_used": dict,
            "cached_at": str,
            "from_cache": bool,
        }
    """
    if not gene:
        return {"gene": "", "pmids": [], "articles": [], "total_found": 0, "from_cache": False}

    cache_key = _make_cache_key(gene, hgvsc, hgvsp)

    # 캐시 확인
    if not force_refresh and _DB_PATH:
        cached = _get_search_cache(cache_key)
        if cached:
            cached_arts, uncached_pmids = _get_cached_pmids(cached["pmids"])
            if not uncached_pmids:
                # 검색 결과 + 논문 모두 캐시 히트 → PubMed 요청 없음
                cached["articles"] = _score_and_sort(cached_arts, gene, hgvsc, hgvsp)
                logger.info(f"Literature cache hit: {cache_key} ({len(cached_arts)} articles)")
                return cached
            else:
                # 검색 캐시는 있으나 일부 논문 미캐시 → 누락 논문만 EFetch (tier 재검색 없음)
                ncbi_params = {"email": ncbi_email, "api_key": ncbi_api_key, "tool": ncbi_tool}
                new_articles = await _fetch_article_details(uncached_pmids, ncbi_params)
                all_articles = cached_arts + new_articles
                scored = _score_and_sort(all_articles, gene, hgvsc, hgvsp)
                if new_articles and _DB_PATH:
                    async with _DB_LOCK:
                        _save_articles(new_articles)
                logger.info(
                    f"Literature partial cache hit: {cache_key} "
                    f"({len(cached_arts)} cached + {len(new_articles)} new articles)"
                )
                return {
                    "gene": gene, "hgvsc": hgvsc, "hgvsp": hgvsp,
                    "pmids": cached["pmids"], "articles": scored,
                    "total_found": len(scored),
                    "tiers_used": cached.get("tiers_used", {}),
                    "cached_at": cached.get("cached_at", _now_iso()),
                    "from_cache": True,
                }

    # 캐시 없음 → NCBI 전체 tier 검색
    ncbi_params = {"email": ncbi_email, "api_key": ncbi_api_key, "tool": ncbi_tool}

    # 4-tier 검색
    pmids, tiers_used, query_terms = await _multi_tier_search(
        gene, hgvsc, hgvsp, effect, max_results, ncbi_params
    )

    # 논문 상세 조회 (캐시에 없는 것만)
    cached_arts, uncached_pmids = _get_cached_pmids(pmids) if _DB_PATH else ([], pmids)
    new_articles = []
    if uncached_pmids:
        new_articles = await _fetch_article_details(uncached_pmids, ncbi_params)

    all_articles = cached_arts + new_articles

    # 관련성 점수 + 정렬
    scored = _score_and_sort(all_articles, gene, hgvsc, hgvsp)

    # 캐시 저장
    if _DB_PATH:
        async with _DB_LOCK:
            _save_articles(new_articles)
            _save_search_cache(cache_key, gene, hgvsc, hgvsp, query_terms, pmids, tiers_used)

    return {
        "gene": gene,
        "hgvsc": hgvsc,
        "hgvsp": hgvsp,
        "pmids": pmids,
        "articles": scored,
        "total_found": len(scored),
        "tiers_used": tiers_used,
        "cached_at": _now_iso(),
        "from_cache": False,
    }


# ══════════════════════════════════════════════════════════════
# Multi-tier 검색 전략
# ══════════════════════════════════════════════════════════════

async def _multi_tier_search(
    gene: str, hgvsc: str, hgvsp: str, effect: str,
    max_results: int, ncbi_params: Dict,
) -> Tuple[List[str], Dict, str]:
    """
    Variant-specific 우선 검색 전략:
      Tier 1: gene + hgvsc (가장 정밀 → 한도 없이 모두 수집)
      Tier 2: gene + hgvsp
      Tier 3: gene + effect + carrier (변이 유형 기반)
      Tier 4: gene alias + 임상 키워드 (variant-specific 결과가 없을 때만)
      Europe PMC: 보완 검색
      Tier 5: 넓은 범위 (결과 극히 부족할 때만)
    최종 max_results 개 반환 (relevance 점수 순 정렬은 caller에서 수행)
    """
    tiers: Dict[str, int] = {}
    seen: dict = {}

    def add_pmids(tier_name: str, new_pmids: List[str]) -> None:
        added = 0
        for p in new_pmids:
            if p not in seen:
                seen[p] = True
                added += 1
        tiers[tier_name] = added

    variant_specific_found = 0  # Tier1/2 히트 수 추적

    # Tier 1: gene + hgvsc (변이 정확 매칭, 논문이 많을수록 좋음 → 상한 2×max)
    t1_query = ""
    if hgvsc:
        bare = hgvsc.replace("c.", "")
        t1_query = (
            f'"{gene}"[Title/Abstract] AND '
            f'("{hgvsc}"[Title/Abstract] OR "{bare}"[Title/Abstract])'
        )
        t1 = await _esearch(t1_query, max_results * 2, ncbi_params)
        add_pmids("tier1_hgvsc", t1)
        variant_specific_found += len(t1)

    # Tier 2: gene + hgvsp
    if hgvsp:
        t2_query = f'"{gene}"[Title/Abstract] AND "{hgvsp}"[Title/Abstract]'
        t2 = await _esearch(t2_query, max_results * 2, ncbi_params)
        add_pmids("tier2_hgvsp", t2)
        variant_specific_found += len(t2)

    # Tier 3: gene + effect 키워드 + carrier (변이 유형 기반, 여전히 variant-flavored)
    effect_kw = _EFFECT_KEYWORDS.get(effect.split("&")[0].strip(), "")
    if effect_kw and len(seen) < max_results:
        t3_query = (
            f'"{gene}"[Title/Abstract] AND '
            f'"{effect_kw}"[Title/Abstract] AND "carrier"[Title/Abstract]'
        )
        t3 = await _esearch(t3_query, max_results, ncbi_params)
        add_pmids("tier3_effect", t3)

    # Europe PMC: variant-specific 결과 보완
    if len(seen) < max_results:
        epmc = await _search_europe_pmc(gene, hgvsc or hgvsp, max_results)
        add_pmids("europe_pmc", epmc)

    # Tier 4: gene alias — variant-specific 히트가 없을 때만 (너무 넓어지는 것 방지)
    if variant_specific_found == 0 and len(seen) < max_results:
        aliases = await _fetch_gene_aliases(gene)
        if aliases:
            alias_str = " OR ".join(f'"{a}"[Title/Abstract]' for a in aliases[:3])
            t4_query = (
                f'({alias_str}) AND '
                f'("carrier screening" OR "hereditary" OR "pathogenic")'
            )
            t4 = await _esearch(t4_query, max_results, ncbi_params)
            add_pmids("tier4_alias", t4)

    # Tier 5: 최후 넓은 범위 (결과 극히 부족할 때)
    if len(seen) < 3:
        t5_query = (
            f'"{gene}"[Title/Abstract] AND '
            f'("pathogenic variant" OR "germline variant" OR "disease-causing")'
        )
        t5 = await _esearch(t5_query, max_results, ncbi_params)
        add_pmids("tier5_broad", t5)

    final_pmids = list(seen.keys())[:max_results]
    main_query = t1_query or f'"{gene}" variant search'

    logger.info(
        f"Literature search [{gene} {hgvsc or hgvsp}]: "
        f"{len(final_pmids)} PMIDs (variant-specific={variant_specific_found}) | tiers={tiers}"
    )

    return final_pmids, tiers, main_query


# ══════════════════════════════════════════════════════════════
# NCBI E-utilities
# ══════════════════════════════════════════════════════════════

async def _ncbi_get(url: str, params: Dict, timeout: int = 30) -> Optional[object]:
    """NCBI HTTP GET with semaphore + 429 retry backoff."""
    import httpx
    async with _NCBI_SEMAPHORE:
        for attempt, delay in enumerate([0.0] + list(_NCBI_RETRY_DELAYS)):
            if delay:
                await asyncio.sleep(delay)
            try:
                async with httpx.AsyncClient(timeout=timeout) as client:
                    resp = await client.get(url, params=params)
                    if resp.status_code == 429:
                        logger.debug(f"NCBI 429 rate-limit (attempt {attempt+1}), retrying in {_NCBI_RETRY_DELAYS[min(attempt, len(_NCBI_RETRY_DELAYS)-1)]}s…")
                        continue
                    resp.raise_for_status()
                    return resp
            except httpx.HTTPStatusError:
                raise
            except Exception as e:
                if attempt < len(_NCBI_RETRY_DELAYS):
                    continue
                raise
    return None


async def _esearch(query: str, max_results: int, ncbi_params: Dict) -> List[str]:
    """NCBI ESearch → PMID 목록."""
    try:
        params = {
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "retmode": "json",
            "tool": ncbi_params.get("tool", "service-daemon"),
            "email": ncbi_params.get("email", "service-daemon@example.com"),
        }
        if ncbi_params.get("api_key"):
            params["api_key"] = ncbi_params["api_key"]

        resp = await _ncbi_get(f"{NCBI_BASE}/esearch.fcgi", params)
        if resp is None:
            return []
        data = resp.json()
        return data.get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        logger.warning(f"ESearch failed [{query[:60]}...]: {e}")
        return []


async def _fetch_article_details(pmids: List[str], ncbi_params: Dict) -> List[Dict[str, Any]]:
    """NCBI EFetch → 논문 상세 정보."""
    if not pmids:
        return []
    try:
        params = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "xml",
            "tool": ncbi_params.get("tool", "service-daemon"),
            "email": ncbi_params.get("email", "service-daemon@example.com"),
        }
        if ncbi_params.get("api_key"):
            params["api_key"] = ncbi_params["api_key"]

        resp = await _ncbi_get(f"{NCBI_BASE}/efetch.fcgi", params, timeout=60)
        if resp is None:
            return []
        return _parse_pubmed_xml(resp.text)
    except Exception as e:
        logger.warning(f"EFetch failed for {len(pmids)} PMIDs: {e}")
        return []


def _parse_pubmed_xml(xml_text: str) -> List[Dict[str, Any]]:
    """PubMed XML → 논문 딕셔너리 목록."""
    articles = []
    try:
        root = ET.fromstring(xml_text)
        for pa in root.findall(".//PubmedArticle"):
            art = _parse_single_article(pa)
            if art:
                articles.append(art)
    except Exception as e:
        logger.warning(f"PubMed XML parse error: {e}")
    return articles


def _parse_single_article(pa: ET.Element) -> Optional[Dict[str, Any]]:
    try:
        pmid_el = pa.find(".//PMID")
        pmid = pmid_el.text if pmid_el is not None else ""
        if not pmid:
            return None

        article_el = pa.find(".//Article")
        if article_el is None:
            return None

        title_el = article_el.find(".//ArticleTitle")
        title = "".join(title_el.itertext()) if title_el is not None else ""

        abstract_parts = []
        for at in article_el.findall(".//AbstractText"):
            label = at.get("Label", "")
            text = "".join(at.itertext())
            abstract_parts.append(f"{label}: {text}" if label else text)
        abstract = " ".join(abstract_parts)

        journal_el = article_el.find(".//Title")
        journal = journal_el.text if journal_el is not None else ""

        authors = []
        for author in article_el.findall(".//Author"):
            last = author.findtext("LastName", "")
            fore = author.findtext("ForeName", "")
            if last:
                authors.append(f"{last} {fore}".strip())

        pub_date = ""
        pd_el = article_el.find(".//PubDate")
        if pd_el is not None:
            year = pd_el.findtext("Year", "")
            month = pd_el.findtext("Month", "")
            pub_date = f"{year} {month}".strip()

        doi = ""
        for aid in pa.findall(".//ArticleId"):
            if aid.get("IdType") == "doi":
                doi = aid.text or ""
                break

        return {
            "pmid": pmid,
            "title": title,
            "abstract": abstract,
            "authors": authors[:5],
            "journal": journal,
            "pub_date": pub_date,
            "doi": doi,
        }
    except Exception:
        return None


# ══════════════════════════════════════════════════════════════
# Europe PMC fallback
# ══════════════════════════════════════════════════════════════

async def _search_europe_pmc(gene: str, variant_id: str, max_results: int) -> List[str]:
    try:
        import httpx
        query = f'"{gene}" AND "{variant_id}"' if variant_id else f'"{gene}" carrier screening pathogenic'
        async with httpx.AsyncClient(timeout=30) as client:
            resp = await client.get(
                f"{EUROPEPMC_BASE}/search",
                params={"query": query, "format": "json", "resultType": "lite", "pageSize": str(max_results)},
            )
            resp.raise_for_status()
            data = resp.json()
            return [
                str(r["pmid"])
                for r in data.get("resultList", {}).get("result", [])
                if r.get("pmid")
            ]
    except Exception as e:
        logger.debug(f"Europe PMC fallback failed: {e}")
        return []


# ══════════════════════════════════════════════════════════════
# Gene Aliases (MyGene.info)
# ══════════════════════════════════════════════════════════════

async def _fetch_gene_aliases(gene: str) -> List[str]:
    try:
        import httpx
        async with httpx.AsyncClient(timeout=10) as client:
            resp = await client.get(
                f"{MYGENE_BASE}/query",
                params={"q": gene, "fields": "alias,symbol", "size": "1"},
            )
            if resp.status_code != 200:
                return []
            hits = resp.json().get("hits", [])
            if not hits:
                return []
            hit = hits[0]
            aliases = hit.get("alias", [])
            if isinstance(aliases, str):
                aliases = [aliases]
            symbol = hit.get("symbol", "")
            if symbol and symbol != gene:
                aliases.append(symbol)
            return [a for a in aliases if a and a != gene][:5]
    except Exception:
        return []


# ══════════════════════════════════════════════════════════════
# 관련성 점수
# ══════════════════════════════════════════════════════════════

def _score_and_sort(
    articles: List[Dict[str, Any]],
    gene: str, hgvsc: str, hgvsp: str,
) -> List[Dict[str, Any]]:
    for art in articles:
        art["relevance_score"] = _calculate_relevance(art, gene, hgvsc, hgvsp)
    return sorted(articles, key=lambda a: a.get("relevance_score", 0), reverse=True)


def _calculate_relevance(art: Dict, gene: str, hgvsc: str, hgvsp: str) -> int:
    score = 0
    title = (art.get("title") or "").lower()
    abstract = (art.get("abstract") or "").lower()
    gene_l = gene.lower()

    if gene_l in title:
        score += 30
    elif gene_l in abstract:
        score += 20

    for v in [hgvsc, hgvsp]:
        if not v:
            continue
        v_l = v.lower()
        if v_l in title:
            score += 25
            break
        elif v_l in abstract:
            score += 15
            break

    for kw in ["carrier", "pathogenic", "hereditary", "germline", "variant"]:
        if kw in title:
            score += 5
        elif kw in abstract:
            score += 2

    return min(max(score, 0), 100)


# ══════════════════════════════════════════════════════════════
# 캐시 관리 (REST API에서 호출)
# ══════════════════════════════════════════════════════════════

def list_cached_articles(
    cursor: int = 0,
    count: int = 50,
    search: str = "",
    sort_by: str = "cached_at",
) -> Dict[str, Any]:
    """캐시된 논문 목록 반환 (페이지네이션 + 검색)."""
    conn = _db()
    cur = conn.cursor()
    cur.execute("SELECT * FROM article_cache ORDER BY cached_at DESC")
    rows = cur.fetchall()
    conn.close()

    articles = [
        {
            "pmid": r["pmid"],
            "title": r["title"],
            "abstract": r["abstract"],
            "authors": json.loads(r["authors"] or "[]"),
            "journal": r["journal"],
            "pub_date": r["pub_date"],
            "doi": r["doi"],
            "cached_at": r["cached_at"],
        }
        for r in rows
    ]

    if sort_by == "pub_date":
        articles.sort(key=lambda a: a.get("pub_date") or "", reverse=True)

    if search:
        q = search.lower()
        articles = [
            a for a in articles
            if q in (a.get("title") or "").lower()
            or q in (a.get("abstract") or "").lower()
            or q in (a.get("pmid") or "").lower()
            or q in " ".join(a.get("authors") or []).lower()
        ]

    total = len(articles)
    page = articles[cursor : cursor + count]
    return {
        "articles": page,
        "total": total,
        "cursor": cursor,
        "next_cursor": cursor + count if cursor + count < total else total,
        "has_more": cursor + count < total,
    }


def get_cached_article(pmid: str) -> Optional[Dict[str, Any]]:
    conn = _db()
    cur = conn.cursor()
    cur.execute("SELECT * FROM article_cache WHERE pmid = ?", (pmid,))
    row = cur.fetchone()
    conn.close()
    if not row:
        return None
    return {
        "pmid": row["pmid"],
        "title": row["title"],
        "abstract": row["abstract"],
        "authors": json.loads(row["authors"] or "[]"),
        "journal": row["journal"],
        "pub_date": row["pub_date"],
        "doi": row["doi"],
        "cached_at": row["cached_at"],
    }


def delete_cached_article(pmid: str) -> bool:
    conn = _db()
    cur = conn.cursor()
    cur.execute("DELETE FROM article_cache WHERE pmid = ?", (pmid,))
    cur.execute("UPDATE search_cache SET pmids = json_remove(pmids, json_each.key) WHERE json_each.value = ? ", (pmid,))
    conn.commit()
    deleted = cur.rowcount > 0
    conn.close()
    return deleted


def get_literature_stats() -> Dict[str, Any]:
    """캐시 통계."""
    conn = _db()
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) as cnt FROM article_cache")
    article_count = cur.fetchone()["cnt"]
    cur.execute("SELECT COUNT(*) as cnt FROM search_cache")
    search_count = cur.fetchone()["cnt"]
    cur.execute("SELECT DISTINCT gene FROM search_cache ORDER BY gene")
    genes = [r["gene"] for r in cur.fetchall()]
    conn.close()
    return {
        "total_articles": article_count,
        "total_searches": search_count,
        "unique_genes": len(genes),
        "genes": genes[:50],
    }


def list_cached_searches(
    gene: str = "",
    cursor: int = 0,
    count: int = 50,
) -> Dict[str, Any]:
    """캐시된 검색 목록 (variant 별)."""
    conn = _db()
    cur = conn.cursor()
    if gene:
        cur.execute(
            "SELECT * FROM search_cache WHERE gene = ? ORDER BY cached_at DESC LIMIT ? OFFSET ?",
            (gene.upper(), count, cursor),
        )
    else:
        cur.execute(
            "SELECT * FROM search_cache ORDER BY cached_at DESC LIMIT ? OFFSET ?",
            (count, cursor),
        )
    rows = cur.fetchall()
    cur.execute("SELECT COUNT(*) as cnt FROM search_cache" + (" WHERE gene=?" if gene else ""),
                (gene.upper(),) if gene else ())
    total = cur.fetchone()["cnt"]
    conn.close()

    searches = [
        {
            "cache_key": r["cache_key"],
            "gene": r["gene"],
            "hgvsc": r["hgvsc"],
            "hgvsp": r["hgvsp"],
            "total_found": r["total_found"],
            "cached_at": r["cached_at"],
        }
        for r in rows
    ]
    return {"searches": searches, "total": total, "cursor": cursor, "has_more": cursor + count < total}


def clear_literature_cache() -> Dict[str, int]:
    """전체 캐시 삭제 (관리자 기능)."""
    conn = _db()
    cur = conn.cursor()
    cur.execute("DELETE FROM article_cache")
    articles_deleted = cur.rowcount
    cur.execute("DELETE FROM search_cache")
    searches_deleted = cur.rowcount
    conn.commit()
    conn.close()
    logger.warning(f"Literature cache cleared: {articles_deleted} articles, {searches_deleted} searches")
    return {"articles_deleted": articles_deleted, "searches_deleted": searches_deleted}
