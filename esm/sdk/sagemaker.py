import json

import boto3

from esm.sdk.forge import (
    ESM3ForgeInferenceClient,
    SequenceStructureForgeInferenceClient,
)


class SequenceStructureSageMakerClient(SequenceStructureForgeInferenceClient):
    def __init__(self, endpoint_name: str):
        """SequenceStructure (folding and inverse folding) client that talks to a SageMaker endpoint.

        Args:
            endpoint_name: Name of the SageMaker endpoint.
        """
        # Dummy URL and token to make SequenceStructureForgeInferenceClient happy.
        super().__init__(url="", token="dummy")

        self._endpoint_name = endpoint_name

        self._client = boto3.client(service_name="sagemaker-runtime")

    def _post(self, endpoint, request, potential_sequence_of_concern):
        request["potential_sequence_of_concern"] = potential_sequence_of_concern
        request["model"] = request.get("model", None)
        invocations_request = {
            # Duplicate these fields at the top level to make Forge requests consistent.
            "model": request["model"],
            "request_id": "",  # Forge specific field.
            "user_id": "",  # Forge specific field.
            # Invocation data bits.
            "api_ver": "v1",  # Must be v1 right now.
            "endpoint": endpoint,
            # Wrapped request.
            endpoint: request,
        }

        try:
            response = self._client.invoke_endpoint(
                EndpointName=self._endpoint_name,
                ContentType="application/json",
                Body=json.dumps(invocations_request),
            )
        except Exception as e:
            raise RuntimeError(f"Failure in {endpoint}: {e}") from e

        data = json.loads(response["Body"].read().decode())

        # Response must match request.
        assert (
            data["endpoint"] == endpoint
        ), f"Response endpoint is {data['endpoint']} but request is {endpoint}"

        # Get the actual responses under the endpoint key.
        data = data[endpoint]

        return data


class ESM3SageMakerClient(ESM3ForgeInferenceClient):
    def __init__(self, endpoint_name: str, model: str):
        """ESM3 client that talks to a SageMaker endpoint.

        Args:
            endpoint_name: Name of the SageMaker endpoint.
            model: Name of the ESM3 model.
        """
        # Dummy URL and token to make ESM3ForgeInferenceClient happy.
        super().__init__(model=model, url="", token="dummy")

        self._endpoint_name = endpoint_name
        self._model = model

        self._client = boto3.client(service_name="sagemaker-runtime")

    def _post(self, endpoint, request, potential_sequence_of_concern):
        request["potential_sequence_of_concern"] = potential_sequence_of_concern

        invocations_request = {
            # Duplicate these fields at the top level to make Forge requests consistent.
            "model": request["model"],
            "request_id": "",  # Forge specific field.
            "user_id": "",  # Forge specific field.
            # Invocation data bits.
            "api_ver": "v1",  # Must be v1 right now.
            "endpoint": endpoint,
            # Wrapped request.
            endpoint: request,
        }

        try:
            response = self._client.invoke_endpoint(
                EndpointName=self._endpoint_name,
                ContentType="application/json",
                Body=json.dumps(invocations_request),
            )
        except Exception as e:
            raise RuntimeError(f"Failure in {endpoint}: {e}")

        data = json.loads(response["Body"].read().decode())

        # Response must match request.
        assert data["endpoint"] == endpoint

        # Get the actual responses under the endpoint key.
        data = data[endpoint]

        return data
