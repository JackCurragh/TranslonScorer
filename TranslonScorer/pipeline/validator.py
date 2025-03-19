
from .config import Config
from ..utils.logging import log_error, log_info
import click 

def validate_config(config: Config) -> None:
    """
    Validate the configuration, raising appropriate exceptions for invalid settings.
    This is a wrapper around Config.validate() that logs errors.
    """
    try:
        config.validate()
        log_info("Configuration validated successfully.")
        return True
    except (ValueError, FileNotFoundError) as e:
        log_error(f"Configuration error: {str(e)}")
        raise
    except Exception as e:
        log_error(f"Unexpected error validating configuration: {str(e)}")
        raise