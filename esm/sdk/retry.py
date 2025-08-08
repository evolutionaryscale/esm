import inspect
from contextvars import ContextVar
from functools import wraps

import httpx
from tenacity import (
    retry,
    retry_if_exception_type,
    retry_if_result,
    stop_after_attempt,
    wait_incrementing,
)

from esm.sdk.api import ESMProteinError

skip_retries_var = ContextVar("skip_retries", default=False)


def retry_if_specific_error(exception):
    """
    We only retry on specific errors.
    """
    return isinstance(exception, ESMProteinError) and exception.error_code in {
        429,
        500,
        502,
        504,
        500,
    }


def log_retry_attempt(retry_state):
    print(
        f"Retrying... Attempt {retry_state.attempt_number} after {retry_state.next_action.sleep}s due to: {retry_state.outcome.result()}"
    )


def retry_decorator(func):
    """
    A static method that returns a retry decorator. This decorator uses the
    instance's retry settings.
    """

    @wraps(func)
    async def async_wrapper(instance, *args, **kwargs):
        if skip_retries_var.get():
            return await func(instance, *args, **kwargs)
        retry_decorator = retry(
            retry=retry_if_result(retry_if_specific_error)
            | retry_if_exception_type(httpx.ConnectTimeout),  # ADDED
            wait=wait_incrementing(
                increment=1, start=instance.min_retry_wait, max=instance.max_retry_wait
            ),
            stop=stop_after_attempt(instance.max_retry_attempts),
            before_sleep=log_retry_attempt,
        )
        # Apply the retry decorator to the function
        return await retry_decorator(func)(instance, *args, **kwargs)

    @wraps(func)
    def wrapper(instance, *args, **kwargs):
        if skip_retries_var.get():
            return func(instance, *args, **kwargs)
        retry_decorator = retry(
            retry=retry_if_result(retry_if_specific_error)
            | retry_if_exception_type(httpx.ConnectTimeout),  # ADDED
            wait=wait_incrementing(
                increment=1, start=instance.min_retry_wait, max=instance.max_retry_wait
            ),
            stop=stop_after_attempt(instance.max_retry_attempts),
            before_sleep=log_retry_attempt,
        )
        # Apply the retry decorator to the function
        return retry_decorator(func)(instance, *args, **kwargs)

    return async_wrapper if inspect.iscoroutinefunction(func) else wrapper
