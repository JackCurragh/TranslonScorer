"""
Logging configuration for TranslonScorer.

This module provides consistent logging functionality across the package,
including error handling and informational messages.
"""

import logging
import sys
from typing import Optional, Type


def setup_logging(level: int = logging.INFO) -> None:
    """
    Set up logging configuration for the package.

    Args:
        level: The logging level to use (default: INFO)
    """
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
    )


def log_info(message: str) -> None:
    """
    Log an informational message.

    Args:
        message: The message to log
    """
    logging.info(message)


def log_warning(message: str) -> None:
    """
    Log a warning message.

    Args:
        message: The warning message to log
    """
    logging.warning(message)


def log_error(
    message: str,
    raise_exception: bool = True,
    exception_type: Optional[Type[Exception]] = RuntimeError,
) -> None:
    """
    Log an error message and optionally raise an exception.

    Args:
        message: The error message to log
        raise_exception: Whether to raise an exception (default: True)
        exception_type: Type of exception to raise (default: RuntimeError)

    Raises:
        The specified exception type with the error message if raise_exception is True
    """
    logging.error(message)
    if raise_exception:
        raise exception_type(message)


# Set up logging when the module is imported
setup_logging()
