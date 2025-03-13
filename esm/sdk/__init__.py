import os

from esm.sdk.forge import ESM3ForgeInferenceClient
from esm.utils.forge_context_manager import ForgeBatchExecutor

# Note: please do not import ESM3SageMakerClient here since that requires AWS SDK.


def client(
    model="esm3-sm-open-v1",
    url="https://forge.evolutionaryscale.ai",
    token=os.environ.get("ESM_API_KEY", ""),
    request_timeout=None,
):
    """
    Args:
        model: Name of the model to use.
        url: URL of a forge server.
        token: User's API token.
        request_timeout: Amount of time to wait for a request to finish.
            Default is wait indefinitely.
    """
    return ESM3ForgeInferenceClient(model, url, token, request_timeout)


def batch_executor(max_attempts: int = 10):
    """
    Args:
        max_attempts: Maximum number of attempts to make before giving up.

    Usage:
        with batch_executor() as executor:
            executor.execute_batch(fn, **kwargs)
    """
    return ForgeBatchExecutor(max_attempts)
