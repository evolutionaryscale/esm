from urllib.parse import urljoin

import requests

from esm.sdk.forge import ESM3ForgeInferenceClient


class ESM3SageMakerClient(ESM3ForgeInferenceClient):
    # TODO(jungong) : get rid of url parameter after we switch to SageMaker Runtime.
    def __init__(self, model: str, endpoint_name: str, url: str):
        # Dummy URL and token to make ESM3ForgeInferenceClient happy.
        super().__init__(model=model, url=url, token="sagemaker")

        self._endpoint_name = endpoint_name

    def __post(self, endpoint, request, potential_sequence_of_concern):
        request["potential_sequence_of_concern"] = potential_sequence_of_concern

        invocations_request = {
            # Duplicate these fields at the top level to make Forge requests consistent.
            "model": request["model"],
            "request_id": request["request_id"],
            "user_id": request["user_id"],
            # Invocation data bits.
            "api_ver": "v1",  # Must be v1 right now.
            "endpoint": endpoint,
            # Wrapped request.
            endpoint: request,
        }

        # Note: this does NOT actually work on SageMaker.
        # This temporary implementation is for testing only for now.
        # TODO(jungong) : switch to sagemaker_runtime.invoke_endpoint()
        # after I figure out how to test this on SageMaker for real.
        response = requests.post(
            urljoin(self.url, "/invocations"),  # SageMaker's one & only entrypoint.
            json=invocations_request,
            headers=self.headers,
            timeout=self.request_timeout,
        )

        if not response.ok:
            raise RuntimeError(f"Failure in {endpoint}: {response.text}")

        data = response.json()

        # Response must match request.
        assert data["endpoint"] == endpoint

        # Get the actual responses under the endpoint key.
        data = data[endpoint]

        return data
