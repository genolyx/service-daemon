"""
Carrier Screening Service Plugin Package

Nextflow 파이프라인 실행 → VCF Annotation → 리뷰 데이터 생성 → 리포트 생성까지
Carrier Screening 서비스의 전체 워크플로우를 담당합니다.

모듈 구성:
    plugin.py       - ServicePlugin 구현체 (메인 오케스트레이터)
    annotator.py    - VCF annotation 엔진 (ClinVar, gnomAD, snpEff, ClinGen, ACMG)
    vcf_parser.py   - VCF 파싱 및 변이 추출
    review.py       - 리뷰 페이지 JSON 생성
    report.py       - 최종 PDF 리포트 생성
    acmg.py         - ACMG 변이 분류 (Lite + AI 모드)
"""

from .plugin import CarrierScreeningPlugin


def create_plugin():
    """플러그인 팩토리 함수 (service registry에서 호출)"""
    return CarrierScreeningPlugin()
