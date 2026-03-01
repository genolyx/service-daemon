"""
Service Plugin Registry

서비스 플러그인을 동적으로 로딩하고 관리합니다.
"""

import logging
from typing import Dict, Optional

from .base import ServicePlugin

logger = logging.getLogger(__name__)

# 전역 플러그인 레지스트리
_registry: Dict[str, ServicePlugin] = {}


def register_plugin(service_code: str, plugin: ServicePlugin):
    """플러그인 등록"""
    _registry[service_code] = plugin
    logger.info(f"Registered service plugin: {service_code} ({plugin.__class__.__name__})")


def get_plugin(service_code: str) -> Optional[ServicePlugin]:
    """service_code로 플러그인 조회"""
    return _registry.get(service_code)


def get_all_plugins() -> Dict[str, ServicePlugin]:
    """등록된 모든 플러그인 반환"""
    return dict(_registry)


def list_service_codes() -> list:
    """등록된 서비스 코드 목록"""
    return list(_registry.keys())


def load_plugins(enabled_services: list):
    """
    활성화된 서비스의 플러그인을 동적으로 로딩합니다.
    
    각 서비스 플러그인은 app/services/ 디렉토리 내에 모듈로 존재하며,
    모듈 내에 create_plugin() 팩토리 함수를 구현해야 합니다.
    """
    import importlib

    for service_code in enabled_services:
        module_name = f"app.services.{service_code}"
        try:
            module = importlib.import_module(module_name)
            if hasattr(module, "create_plugin"):
                plugin = module.create_plugin()
                register_plugin(service_code, plugin)
            else:
                logger.warning(
                    f"Service module '{module_name}' does not have create_plugin(). Skipping."
                )
        except ImportError as e:
            logger.warning(f"Could not load service plugin '{service_code}': {e}")
        except Exception as e:
            logger.error(f"Error loading service plugin '{service_code}': {e}", exc_info=True)

    logger.info(f"Loaded {len(_registry)} service plugin(s): {list(_registry.keys())}")
