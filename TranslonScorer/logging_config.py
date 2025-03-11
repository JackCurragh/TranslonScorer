"""Logging configuration for TranslonScorer."""

import logging
import os
from typing import Optional

# Create the base logger
logger = logging.getLogger("TranslonScorer")
logger.setLevel(logging.INFO)

# Create formatters
CONSOLE_FORMAT = '%(message)s'
FILE_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

console_formatter = logging.Formatter(CONSOLE_FORMAT)
file_formatter = logging.Formatter(FILE_FORMAT)

# Console handler
console_handler = logging.StreamHandler()
console_handler.setFormatter(console_formatter)
console_handler.setLevel(logging.INFO)
logger.addHandler(console_handler)

# Global flag for logging state
LOGGING_ENABLED = True

def configure_logging(enable: bool = True, log_file: Optional[str] = None, log_level: str = "INFO") -> None:
    """
    Configure the logging system.

    Args:
        enable (bool): Whether to enable logging
        log_file (str, optional): Path to log file. If provided, logs will be written to this file
        log_level (str): Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    """
    global LOGGING_ENABLED
    LOGGING_ENABLED = enable
    
    # Set log level
    level = getattr(logging, log_level.upper())
    logger.setLevel(level)
    
    if enable:
        # Clear existing handlers
        logger.handlers = []
        
        # Always add console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(console_formatter)
        console_handler.setLevel(level)
        logger.addHandler(console_handler)
        
        # Add file handler if log file specified
        if log_file:
            # Create log directory if it doesn't exist
            log_dir = os.path.dirname(log_file)
            if log_dir and not os.path.exists(log_dir):
                os.makedirs(log_dir)
                
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(file_formatter)
            file_handler.setLevel(level)
            logger.addHandler(file_handler)
    else:
        logger.setLevel(logging.WARNING)

def log_info(message: str) -> None:
    """Log an info message if logging is enabled."""
    if LOGGING_ENABLED:
        logger.info(message)

def log_warning(message: str) -> None:
    """Log a warning message if logging is enabled."""
    if LOGGING_ENABLED:
        logger.warning(message)

def log_error(message: str, raise_exception: bool = True, exception_type: type = RuntimeError) -> None:
    """
    Log an error message and optionally raise an exception.

    Args:
        message: The error message to log
        raise_exception: Whether to raise an exception after logging (default: True)
        exception_type: Type of exception to raise (default: RuntimeError)

    Raises:
        The specified exception type with the error message if raise_exception is True
    """
    if LOGGING_ENABLED:
        logger.error(message)
    if raise_exception:
        raise exception_type(message)

def log_debug(message: str) -> None:
    """Log a debug message if logging is enabled."""
    if LOGGING_ENABLED:
        logger.debug(message) 