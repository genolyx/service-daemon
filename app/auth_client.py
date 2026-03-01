"""
Authentication Client

플랫폼 API 인증을 관리합니다.
기존 nipt-daemon의 auth_client.py를 일반화하여 재사용합니다.
"""

import time
import httpx
import asyncio
import logging
from typing import Dict, Optional
from dataclasses import dataclass

from .config import settings

logger = logging.getLogger(__name__)


@dataclass
class TokenInfo:
    """토큰 정보"""
    token: str
    expiry: float

    @property
    def is_valid(self) -> bool:
        """토큰 유효 여부 (30초 여유)"""
        return time.time() < self.expiry - 30

    @property
    def is_expired(self) -> bool:
        return not self.is_valid


@dataclass
class AuthTokens:
    """인증 토큰 세트"""
    access_token: TokenInfo
    refresh_token: Optional[str] = None

    @property
    def is_access_token_valid(self) -> bool:
        return self.access_token.is_valid

    @property
    def has_refresh_token(self) -> bool:
        return self.refresh_token is not None


class AuthClient:
    """플랫폼 인증 클라이언트 (accessToken/refreshToken 지원)"""

    def __init__(self):
        self._auth_tokens: Optional[AuthTokens] = None
        self._refresh_lock = asyncio.Lock()

    async def get_token(self, force_refresh: bool = False) -> str:
        """유효한 accessToken 획득"""
        if not self._auth_tokens or force_refresh:
            return await self._perform_login()

        if self._auth_tokens.is_access_token_valid:
            return self._auth_tokens.access_token.token

        if self._auth_tokens.has_refresh_token:
            async with self._refresh_lock:
                if self._auth_tokens and self._auth_tokens.is_access_token_valid:
                    return self._auth_tokens.access_token.token
                success = await self._refresh_access_token()
                if success and self._auth_tokens:
                    return self._auth_tokens.access_token.token

        return await self._perform_login()

    async def _perform_login(self) -> str:
        """username/password 로그인 (재시도 포함)"""
        logger.info("Performing login with username/password")
        max_retries = 3
        last_error = None

        for attempt in range(max_retries):
            try:
                async with httpx.AsyncClient() as client:
                    response = await client.post(
                        (settings.auth_url or "").strip(),
                        headers={"Content-Type": "application/json"},
                        json={
                            "username": (settings.api_username or "").strip(),
                            "password": (settings.api_password or "").strip()
                        },
                        timeout=10.0
                    )
                    response.raise_for_status()

                token_data = response.json()
                access_token = self._extract_access_token(token_data)
                expires_in = self._extract_expires_in(token_data)
                refresh_token = self._extract_refresh_token(token_data)

                self._auth_tokens = AuthTokens(
                    access_token=TokenInfo(
                        token=access_token,
                        expiry=time.time() + expires_in
                    ),
                    refresh_token=refresh_token
                )

                if attempt > 0:
                    logger.info(f"NETWORK_RETRY_SUCCESS: Login succeeded after {attempt + 1} attempts")

                logger.info(
                    f"Login successful (access expires in {expires_in}s, "
                    f"refresh: {'Yes' if refresh_token else 'No'})"
                )
                return self._auth_tokens.access_token.token

            except (httpx.ConnectError, httpx.TimeoutException, OSError) as e:
                last_error = e
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt
                    logger.warning(
                        f"NETWORK_RETRY: {type(e).__name__} during login "
                        f"(attempt {attempt + 1}/{max_retries}): {e}. "
                        f"Retrying in {wait_time}s..."
                    )
                    await asyncio.sleep(wait_time)
                else:
                    logger.error(f"NETWORK_RETRY_FAILED: Login failed after {max_retries} attempts: {e}")
                    self._auth_tokens = None
                    raise

            except Exception as e:
                logger.error(f"Login failed: {e}")
                self._auth_tokens = None
                raise

        if last_error:
            self._auth_tokens = None
            raise last_error
        raise RuntimeError("Unexpected error in login")

    async def _refresh_access_token(self) -> bool:
        """refreshToken으로 accessToken 갱신"""
        if not self._auth_tokens or not self._auth_tokens.has_refresh_token:
            logger.warning("No refresh token available for token refresh")
            return False

        logger.info("Refreshing access token using refresh token")
        max_retries = 3

        for attempt in range(max_retries):
            try:
                refresh_url = f"{settings.platform_api_base}/auth/refresh"
                async with httpx.AsyncClient() as client:
                    response = await client.post(
                        refresh_url,
                        headers={"Content-Type": "application/json"},
                        json={"refreshToken": self._auth_tokens.refresh_token},
                        timeout=10.0
                    )
                    response.raise_for_status()

                token_data = response.json()
                access_token = self._extract_access_token(token_data)
                expires_in = self._extract_expires_in(token_data)
                new_refresh_token = self._extract_refresh_token(token_data)

                if new_refresh_token:
                    self._auth_tokens.refresh_token = new_refresh_token

                self._auth_tokens.access_token = TokenInfo(
                    token=access_token,
                    expiry=time.time() + expires_in
                )

                if attempt > 0:
                    logger.info(f"NETWORK_RETRY_SUCCESS: Token refresh succeeded after {attempt + 1} attempts")

                logger.info(f"Token refresh successful (expires in {expires_in}s)")
                return True

            except (httpx.ConnectError, httpx.TimeoutException, OSError) as e:
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt
                    logger.warning(
                        f"NETWORK_RETRY: {type(e).__name__} during token refresh "
                        f"(attempt {attempt + 1}/{max_retries}): {e}. "
                        f"Retrying in {wait_time}s..."
                    )
                    await asyncio.sleep(wait_time)
                else:
                    logger.error(f"NETWORK_RETRY_FAILED: Token refresh failed after {max_retries} attempts: {e}")
                    self._auth_tokens = None
                    return False

            except Exception as e:
                logger.error(f"Token refresh failed: {e}")
                self._auth_tokens = None
                return False

        return False

    async def auth_request(self, method: str, url: str, **kwargs) -> httpx.Response:
        """인증이 포함된 HTTP 요청 (401 재시도 + 네트워크 재시도)"""
        max_network_retries = 3

        for attempt in range(max_network_retries):
            try:
                response = await self._make_request(method, url, **kwargs)

                if response.status_code == 401:
                    logger.info("Received 401, attempting token refresh and retry")
                    async with self._refresh_lock:
                        try:
                            await self.get_token(force_refresh=True)
                            response = await self._make_request(method, url, **kwargs)
                        except Exception as e:
                            logger.error(f"Token refresh and retry failed: {e}")
                            response.raise_for_status()
                            return response

                if attempt > 0:
                    logger.info(
                        f"NETWORK_RETRY_SUCCESS: Request to {url} recovered after {attempt + 1} attempts"
                    )

                response.raise_for_status()
                return response

            except (httpx.ConnectError, httpx.TimeoutException, OSError) as e:
                if attempt < max_network_retries - 1:
                    wait_time = 2 ** attempt
                    logger.warning(
                        f"NETWORK_RETRY: {type(e).__name__} for {method} {url} "
                        f"(attempt {attempt + 1}/{max_network_retries}): {e}. "
                        f"Retrying in {wait_time}s..."
                    )
                    await asyncio.sleep(wait_time)
                else:
                    logger.error(
                        f"NETWORK_RETRY_FAILED: {type(e).__name__} for {method} {url} "
                        f"after {max_network_retries} attempts: {e}"
                    )
                    raise

        raise RuntimeError(f"Unexpected error in auth_request for {url}")

    async def _make_request(self, method: str, url: str, **kwargs) -> httpx.Response:
        """실제 HTTP 요청 수행"""
        token = await self.get_token()
        headers = kwargs.pop("headers", {})
        headers["Authorization"] = f"Bearer {token}"
        timeout = kwargs.pop("timeout", 10.0)

        async with httpx.AsyncClient(timeout=timeout) as client:
            response = await client.request(method, url, headers=headers, **kwargs)

        if response.status_code >= 400:
            logger.error(f"[{response.status_code}] Error Response from server: {response.text}")

        return response

    async def get_auth_headers(self) -> Dict[str, str]:
        """인증 헤더 생성"""
        token = await self.get_token()
        return {
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json",
            "User-Agent": "Service-Daemon/1.0"
        }

    async def logout(self):
        """세션 종료"""
        if not self._auth_tokens:
            return
        try:
            logout_url = f"{settings.platform_api_base}/auth/logout"
            async with httpx.AsyncClient(timeout=10.0) as client:
                await client.post(
                    logout_url,
                    headers={"Authorization": f"Bearer {self._auth_tokens.access_token.token}"}
                )
            logger.info("Successfully logged out from server")
        except Exception as e:
            logger.warning(f"Logout error: {e}")
        finally:
            self._auth_tokens = None

    def clear_token(self):
        """토큰 강제 클리어"""
        self._auth_tokens = None

    @property
    def has_valid_token(self) -> bool:
        return self._auth_tokens is not None and self._auth_tokens.is_access_token_valid

    @property
    def token_expires_in(self) -> Optional[int]:
        if not self._auth_tokens:
            return None
        return max(0, int(self._auth_tokens.access_token.expiry - time.time()))


# 전역 인스턴스
_auth_client = AuthClient()


async def get_auth_client() -> AuthClient:
    """전역 AuthClient 인스턴스 반환"""
    return _auth_client


async def auth_request(method: str, url: str, **kwargs) -> httpx.Response:
    """편의 함수: 인증 포함 HTTP 요청"""
    return await _auth_client.auth_request(method, url, **kwargs)
