"""
Error Sanitization

Prevents internal exception details from leaking to API consumers.
Logs full errors server-side, returns generic messages to clients.
"""

import logging

logger = logging.getLogger(__name__)


def safe_error_detail(error: Exception, user_message: str = "An internal error occurred") -> str:
    """
    Log full error server-side and return a safe message for the client.

    Args:
        error: The caught exception
        user_message: Generic message to return to the client

    Returns:
        The user_message (never the actual exception details)
    """
    logger.error("Internal error: %s: %s", type(error).__name__, error)
    return user_message
