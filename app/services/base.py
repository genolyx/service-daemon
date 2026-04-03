"""
Service Plugin Base Class

모든 서비스 플러그인이 구현해야 하는 추상 인터페이스를 정의합니다.
"""

import logging
from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional

from app.models import Job, OutputFile

logger = logging.getLogger(__name__)


class ServicePlugin(ABC):
    """
    서비스 플러그인 추상 기반 클래스.
    
    새로운 서비스를 추가하려면 이 클래스를 상속받아
    모든 추상 메소드를 구현하고, 모듈 레벨에서 create_plugin() 팩토리 함수를 제공합니다.
    
    Example:
        class MyServicePlugin(ServicePlugin):
            ...
        
        def create_plugin() -> ServicePlugin:
            return MyServicePlugin()
    """

    @property
    @abstractmethod
    def service_code(self) -> str:
        """서비스 코드 (예: 'carrier_screening', 'nipt')"""
        ...

    @property
    @abstractmethod
    def display_name(self) -> str:
        """서비스 표시 이름 (예: 'Carrier Screening')"""
        ...

    @abstractmethod
    async def prepare_inputs(self, job: Job) -> bool:
        """
        분석에 필요한 입력 파일을 준비합니다.
        
        - FASTQ 파일 다운로드 (URL인 경우)
        - 로컬 경로 확인 및 심볼릭 링크 생성
        - 분석 디렉토리 구조 생성
        
        Args:
            job: 작업 정보
            
        Returns:
            True: 준비 성공, False: 준비 실패
        """
        ...

    @abstractmethod
    async def get_pipeline_command(self, job: Job) -> str:
        """
        실행할 파이프라인의 전체 쉘 명령어를 생성합니다.
        
        Args:
            job: 작업 정보
            
        Returns:
            실행할 쉘 명령어 문자열
        """
        ...

    @abstractmethod
    async def check_completion(self, job: Job) -> bool:
        """
        파이프라인의 성공적인 완료 여부를 확인합니다.
        
        - 결과 파일 존재 여부 확인
        - 로그 파일에서 에러 확인
        - exit code 확인
        
        Args:
            job: 작업 정보
            
        Returns:
            True: 성공적으로 완료, False: 실패 또는 미완료
        """
        ...

    @abstractmethod
    async def process_results(self, job: Job) -> bool:
        """
        파이프라인 완료 후 결과를 후처리합니다.
        
        - Annotation 수행 (ClinVar, gnomAD, dbSNP 등)
        - 리뷰 페이지용 JSON 생성
        - 리포트 PDF 생성
        - 이미지 파일 정리
        
        Args:
            job: 작업 정보
            
        Returns:
            True: 후처리 성공, False: 후처리 실패
        """
        ...

    @abstractmethod
    async def get_output_files(self, job: Job) -> List[OutputFile]:
        """
        플랫폼에 업로드할 최종 결과 파일 목록을 반환합니다.
        
        Args:
            job: 작업 정보
            
        Returns:
            OutputFile 객체 목록
        """
        ...

    def get_pipeline_cwd(self, job: Job) -> Optional[str]:
        """
        ``subprocess`` 실행 시 작업 디렉터리.
        None 이면 runner 가 ``job.analysis_dir`` / ``analysis_base_dir`` 를 사용합니다.
        """
        return None

    def get_progress_stages(self) -> Dict[str, int]:
        """
        파이프라인 진행률 계산을 위한 단계별 퍼센트 매핑.
        
        서비스별로 오버라이드하여 커스텀 진행률을 제공할 수 있습니다.
        
        Returns:
            {process_name: progress_percent} 딕셔너리
        """
        return {}

    async def on_job_start(self, job: Job):
        """작업 시작 시 호출되는 훅 (선택적 오버라이드)"""
        logger.info(f"[{self.service_code}] Job started: {job.order_id}")

    async def on_job_complete(self, job: Job):
        """작업 완료 시 호출되는 훅 (선택적 오버라이드)"""
        logger.info(f"[{self.service_code}] Job completed: {job.order_id}")

    async def on_job_failed(self, job: Job, error: str):
        """작업 실패 시 호출되는 훅 (선택적 오버라이드)"""
        logger.error(f"[{self.service_code}] Job failed: {job.order_id} - {error}")

    async def generate_report(self, job: Job, confirmed_variants: list, reviewer_info: dict) -> bool:
        """
        리뷰어 확정 후 최종 리포트를 생성합니다 (선택적 오버라이드).
        
        Portal에서 리뷰어가 변이를 확정한 후 호출됩니다.
        report.json과 PDF를 생성하여 output 디렉토리에 저장합니다.
        
        Args:
            job: 작업 정보
            confirmed_variants: 리뷰어가 확정한 변이 목록
            reviewer_info: 리뷰어 정보
            
        Returns:
            True: 성공, False: 실패
        """
        logger.warning(f"[{self.service_code}] generate_report not implemented")
        return False

    async def cleanup(self, job: Job):
        """
        작업 완료/실패 후 임시 파일 정리 (선택적 오버라이드).
        
        기본 구현은 아무것도 하지 않습니다.
        """
        pass

    def validate_params(self, params: Dict[str, Any], strict: bool = True) -> tuple:
        """
        서비스별 파라미터 유효성 검사 (선택적 오버라이드).

        strict=False: Save(초안 저장) 단계 — 필수 항목 중 나중에 채울 수 있는 것은 건너뜀.
        strict=True : Submit/Run 단계 — 모든 필수 항목 검사.

        Returns:
            (is_valid: bool, error_message: str)
        """
        return True, ""

    def sync_is_complete(self, job: Job) -> bool:
        """
        daemon 재시작 복구 시 파이프라인이 실제로 완료됐는지 동기적으로 확인 (선택적 오버라이드).

        on-disk 상태만 보고 True/False 반환. 기본 구현은 False.
        True를 반환하면 _restore_from_store 에서 FAILED 대신 REPORT_READY 로 복구합니다.
        """
        return False
