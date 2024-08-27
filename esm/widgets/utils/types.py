from typing import Any, Callable

from esm.sdk.api import ESM3InferenceClient


class ClientInitContainer:
    client_init_callback: Callable[[], ESM3InferenceClient] | None = None

    def __call__(self, *args: Any, **kwds: Any) -> ESM3InferenceClient:
        if self.client_init_callback is None:
            raise ValueError("Client not initialized.")
        return self.client_init_callback()
