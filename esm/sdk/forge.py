import base64
from concurrent.futures import ThreadPoolExecutor
from functools import wraps
from typing import Sequence
from urllib.parse import urljoin

import requests
import torch
from tenacity import retry, retry_if_result, stop_after_attempt, wait_exponential

from esm.sdk.api import (
    ESM3InferenceClient,
    ESMProtein,
    ESMProteinError,
    ESMProteinTensor,
    ForwardAndSampleOutput,
    ForwardTrackData,
    GenerationConfig,
    InverseFoldingConfig,
    LogitsConfig,
    LogitsOutput,
    ProteinType,
    SamplingConfig,
    SamplingTrackConfig,
)
from esm.utils.misc import (
    deserialize_tensors,
    maybe_list,
    maybe_tensor,
)
from esm.utils.sampling import validate_sampling_config
from esm.utils.types import FunctionAnnotation


def _list_to_function_annotations(l) -> list[FunctionAnnotation] | None:
    if l is None or len(l) <= 0:
        return None
    return [FunctionAnnotation(*t) for t in l]


def retry_if_specific_error(exception):
    """
    We only retry on specific errors.
    Currently we retry for 502 (bad gateway) and 429 (rate limit)
    """
    return isinstance(exception, ESMProteinError) and exception.error_code in {
        429,
        502,
        504,
    }


def log_retry_attempt(retry_state):
    print(
        f"Retrying... Attempt {retry_state.attempt_number} after {retry_state.next_action.sleep}s due to: {retry_state.outcome.result()}"
    )


def _validate_protein_tensor_input(input):
    if not isinstance(input, ESMProteinTensor):
        raise ValueError(
            "Input must be an ESMProteinTensor instance. "
            "Use encode() API to encode an ESMProtein into ESMProteinTensor."
        )


class SequenceStructureForgeInferenceClient:
    def __init__(
        self,
        url: str = "https://forge.evolutionaryscale.ai",
        token: str = "",
        request_timeout: int | None = None,
    ):
        if token == "":
            raise RuntimeError(
                "Please provide a token to connect to Forge via token=YOUR_API_TOKEN_HERE"
            )
        self.url = url
        self.token = token
        self.headers = {"Authorization": f"Bearer {self.token}"}
        self.request_timeout = request_timeout

    def fold(
        self,
        sequence: str,
        potential_sequence_of_concern: bool,
        model_name: str | None = None,
    ) -> ESMProtein | ESMProteinError:
        request = {"sequence": sequence}
        if model_name is not None:
            request["model"] = model_name
        try:
            data = self._post("fold", request, potential_sequence_of_concern)
        except ESMProteinError as e:
            return e

        return ESMProtein(
            sequence=sequence,
            coordinates=maybe_tensor(data["coordinates"], convert_none_to_nan=True),
        )

    def inverse_fold(
        self,
        coordinates: torch.Tensor,
        config: InverseFoldingConfig,
        potential_sequence_of_concern: bool,
        model_name: str | None = None,
    ) -> ESMProtein | ESMProteinError:
        inverse_folding_config = {
            "invalid_ids": config.invalid_ids,
            "temperature": config.temperature,
        }
        request = {
            "coordinates": maybe_list(coordinates, convert_nan_to_none=True),
            "inverse_folding_config": inverse_folding_config,
        }
        if model_name is not None:
            request["model"] = model_name
        try:
            data = self._post("inverse_fold", request, potential_sequence_of_concern)
        except ESMProteinError as e:
            return e

        return ESMProtein(sequence=data["sequence"])

    def _post(self, endpoint, request, potential_sequence_of_concern):
        request["potential_sequence_of_concern"] = potential_sequence_of_concern

        response = requests.post(
            urljoin(self.url, f"/api/v1/{endpoint}"),
            json=request,
            headers=self.headers,
            timeout=self.request_timeout,
        )

        if not response.ok:
            raise ESMProteinError(
                error_code=response.status_code,
                error_msg=f"Failure in {endpoint}: {response.text}",
            )

        data = response.json()
        # Nextjs puts outputs dict under "data" key.
        # Lift it up for easier downstream processing.
        if "outputs" not in data and "data" in data:
            data = data["data"]

        # Print warning message if there is any.
        if "warning_messages" in data and data["warning_messages"] is not None:
            for msg in data["warning_messages"]:
                print("\033[31m", msg, "\033[0m")

        return data


class ESM3ForgeInferenceClient(ESM3InferenceClient):
    def __init__(
        self,
        model: str,
        url: str = "https://forge.evolutionaryscale.ai",
        token: str = "",
        request_timeout: int | None = None,
        min_retry_wait: int = 1,
        max_retry_wait: int = 10,
        max_retry_attempts: int = 5,
    ):
        if token == "":
            raise RuntimeError(
                "Please provide a token to connect to Forge via token=YOUR_API_TOKEN_HERE"
            )
        self.model = model  # Name of the model to run.
        self.url = url
        self.token = token
        self.headers = {"Authorization": f"Bearer {self.token}"}
        self.request_timeout = request_timeout
        self.min_retry_wait = min_retry_wait
        self.max_retry_wait = max_retry_wait
        self.max_retry_attempts = max_retry_attempts

    @staticmethod
    def retry_decorator(func):
        """
        A static method that returns a retry decorator. This decorator uses the
        instance's retry settings.
        """

        @wraps(func)
        def wrapper(instance, *args, **kwargs):
            retry_decorator = retry(
                retry=retry_if_result(retry_if_specific_error),
                wait=wait_exponential(
                    multiplier=1,
                    min=instance.min_retry_wait,
                    max=instance.max_retry_wait,
                ),
                stop=stop_after_attempt(instance.max_retry_attempts),
                before_sleep=log_retry_attempt,
            )
            # Apply the retry decorator to the function
            return retry_decorator(func)(instance, *args, **kwargs)

        return wrapper

    @retry_decorator
    def generate(self, input: ProteinType, config: GenerationConfig) -> ProteinType:
        if isinstance(input, ESMProtein):
            output = self.__generate_protein(input, config)
        elif isinstance(input, ESMProteinTensor):
            output = self.__generate_protein_tensor(input, config)
        else:
            return ESMProteinError(
                error_code=500, error_msg=f"Unknown input type {type(input)}"
            )

        if (
            isinstance(output, ESMProtein)
            and isinstance(input, ESMProtein)
            and config.track not in ["function", "residue_annotations"]
        ):
            # Function and residue annotation encoding/decoding is lossy
            # There is no guarantee that decoding encoded tokens will yield the same input
            output.function_annotations = input.function_annotations

        return output

    def batch_generate(
        self, inputs: Sequence[ProteinType], configs: Sequence[GenerationConfig]
    ) -> Sequence[ProteinType]:
        """Forge supports auto-batching. So batch_generate() for the Forge client
        is as simple as running a collection of generate() in parallel using asyncio.
        """
        with ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(self.generate, protein, config)
                for protein, config in zip(inputs, configs)
            ]
            results = []
            for future in futures:
                try:
                    results.append(future.result())
                except Exception as e:
                    results.append(ESMProteinError(500, str(e)))
        return results

    def __generate_protein(
        self, input: ESMProtein, config: GenerationConfig
    ) -> ESMProtein | ESMProteinError:
        req = {}
        req["sequence"] = input.sequence
        req["secondary_structure"] = input.secondary_structure
        req["sasa"] = input.sasa
        if input.function_annotations is not None:
            req["function"] = [x.to_tuple() for x in input.function_annotations]
        req["coordinates"] = maybe_list(input.coordinates, convert_nan_to_none=True)

        request = {
            "model": self.model,
            "inputs": req,
            "track": config.track,
            "invalid_ids": config.invalid_ids,
            "schedule": config.schedule,
            "num_steps": config.num_steps,
            "temperature": config.temperature,
            "top_p": config.top_p,
            "condition_on_coordinates_only": config.condition_on_coordinates_only,
        }
        try:
            data = self._post("generate", request, input.potential_sequence_of_concern)
        except ESMProteinError as e:
            return e

        return ESMProtein(
            sequence=data["outputs"]["sequence"],
            secondary_structure=data["outputs"]["secondary_structure"],
            sasa=data["outputs"]["sasa"],
            function_annotations=_list_to_function_annotations(
                data["outputs"]["function"]
            ),
            coordinates=maybe_tensor(
                data["outputs"]["coordinates"], convert_none_to_nan=True
            ),
            plddt=maybe_tensor(data["outputs"]["plddt"]),
            ptm=maybe_tensor(data["outputs"]["ptm"]),
        )

    def __generate_protein_tensor(
        self, input: ESMProteinTensor, config: GenerationConfig
    ) -> ESMProteinTensor | ESMProteinError:
        req = {}
        req["sequence"] = maybe_list(input.sequence)
        req["structure"] = maybe_list(input.structure)
        req["secondary_structure"] = maybe_list(input.secondary_structure)
        req["sasa"] = maybe_list(input.sasa)
        req["function"] = maybe_list(input.function)
        req["coordinates"] = maybe_list(input.coordinates, convert_nan_to_none=True)
        req["residue_annotation"] = maybe_list(input.residue_annotations)

        request = {
            "model": self.model,
            "inputs": req,
            "track": config.track,
            "invalid_ids": config.invalid_ids,
            "schedule": config.schedule,
            "num_steps": config.num_steps,
            "temperature": config.temperature,
            "top_p": config.top_p,
            "condition_on_coordinates_only": config.condition_on_coordinates_only,
        }

        try:
            data = self._post(
                "generate_tensor", request, input.potential_sequence_of_concern
            )
        except ESMProteinError as e:
            return e

        def _field_to_tensor(field, convert_none_to_nan: bool = False):
            if field not in data["outputs"]:
                return None
            return maybe_tensor(
                data["outputs"][field], convert_none_to_nan=convert_none_to_nan
            )

        output = ESMProteinTensor(
            sequence=_field_to_tensor("sequence"),
            structure=_field_to_tensor("structure"),
            secondary_structure=_field_to_tensor("secondary_structure"),
            sasa=_field_to_tensor("sasa"),
            function=_field_to_tensor("function"),
            residue_annotations=_field_to_tensor("residue_annotation"),
            coordinates=_field_to_tensor("coordinates", convert_none_to_nan=True),
        )

        return output

    @retry_decorator
    def forward_and_sample(
        self, input: ESMProteinTensor, sampling_configuration: SamplingConfig
    ) -> ForwardAndSampleOutput | ESMProteinError:
        _validate_protein_tensor_input(input)
        validate_sampling_config(sampling_configuration, on_invalid="raise")

        req = {}
        sampling_config = {}
        embedding_config = {
            "sequence": sampling_configuration.return_mean_embedding,
            "per_residue": sampling_configuration.return_per_residue_embeddings,
        }

        req["sequence"] = maybe_list(input.sequence)
        req["structure"] = maybe_list(input.structure)
        req["secondary_structure"] = maybe_list(input.secondary_structure)
        req["sasa"] = maybe_list(input.sasa)
        req["function"] = maybe_list(input.function)
        req["coordinates"] = maybe_list(input.coordinates, convert_nan_to_none=True)
        req["residue_annotation"] = maybe_list(input.residue_annotations)

        def do_track(t: str):
            track: SamplingTrackConfig | None
            if (track := getattr(sampling_configuration, t, None)) is None:
                sampling_config[t] = None
            else:
                sampling_config[t] = {
                    "temperature": track.temperature,
                    "top_p": track.top_p,
                    "only_sample_masked_tokens": track.only_sample_masked_tokens,
                    "invalid_ids": track.invalid_ids,
                    "topk_logprobs": track.topk_logprobs,
                }

        do_track("sequence")
        do_track("structure")
        do_track("secondary_structure")
        do_track("sasa")
        do_track("function")

        request = {
            "model": self.model,
            "inputs": req,
            "sampling_config": sampling_config,
            "embedding_config": embedding_config,
        }
        try:
            data = self._post(
                "forward_and_sample", request, input.potential_sequence_of_concern
            )
        except ESMProteinError as e:
            return e

        def get(k, field):
            if data[k] is None:
                return None
            v = data[k][field]
            return torch.tensor(v) if v is not None else None

        tokens = ESMProteinTensor(
            sequence=get("sequence", "tokens"),
            structure=get("structure", "tokens"),
            secondary_structure=get("secondary_structure", "tokens"),
            sasa=get("sasa", "tokens"),
            function=get("function", "tokens"),
        )

        def get_track(field):
            return ForwardTrackData(
                sequence=get("sequence", field),
                structure=get("structure", field),
                secondary_structure=get("secondary_structure", field),
                sasa=get("sasa", field),
                function=get("function", field),
            )

        def operate_on_track(track: ForwardTrackData, fn):
            apply = lambda x: fn(x) if x is not None else None
            return ForwardTrackData(
                sequence=apply(track.sequence),
                structure=apply(track.structure),
                secondary_structure=apply(track.secondary_structure),
                sasa=apply(track.sasa),
                function=apply(track.function),
            )

        logprob = get_track("logprobs")
        output = ForwardAndSampleOutput(
            protein_tensor=tokens,
            logprob=logprob,
            prob=operate_on_track(logprob, torch.exp),
            entropy=get_track("entropy"),
            topk_logprob=get_track("topk_logprobs"),
            topk_tokens=get_track("topk_tokens"),
            per_residue_embedding=data["embeddings"]["per_residue"],
            mean_embedding=data["embeddings"]["sequence"],
        )
        return output

    @retry_decorator
    def encode(self, input: ESMProtein) -> ESMProteinTensor | ESMProteinError:
        tracks = {}
        tracks["sequence"] = input.sequence
        tracks["secondary_structure"] = input.secondary_structure
        tracks["sasa"] = input.sasa
        if input.function_annotations is not None:
            tracks["function"] = [x.to_tuple() for x in input.function_annotations]
        tracks["coordinates"] = maybe_list(input.coordinates, convert_nan_to_none=True)

        request = {"inputs": tracks, "model": self.model}

        try:
            data = self._post("encode", request, input.potential_sequence_of_concern)
        except ESMProteinError as e:
            return e

        return ESMProteinTensor(
            sequence=maybe_tensor(data["outputs"]["sequence"]),
            structure=maybe_tensor(data["outputs"]["structure"]),
            coordinates=maybe_tensor(
                data["outputs"]["coordinates"], convert_none_to_nan=True
            ),
            secondary_structure=maybe_tensor(data["outputs"]["secondary_structure"]),
            sasa=maybe_tensor(data["outputs"]["sasa"]),
            function=maybe_tensor(data["outputs"]["function"]),
            residue_annotations=maybe_tensor(data["outputs"]["residue_annotation"]),
        )

    @retry_decorator
    def decode(self, input: ESMProteinTensor) -> ESMProtein | ESMProteinError:
        _validate_protein_tensor_input(input)

        tokens = {}
        tokens["sequence"] = maybe_list(input.sequence)
        tokens["structure"] = maybe_list(input.structure)
        tokens["secondary_structure"] = maybe_list(input.secondary_structure)
        tokens["sasa"] = maybe_list(input.sasa)
        tokens["function"] = maybe_list(input.function)
        tokens["residue_annotation"] = maybe_list(input.residue_annotations)
        tokens["coordinates"] = maybe_list(input.coordinates, convert_nan_to_none=True)

        request = {"model": self.model, "inputs": tokens}

        try:
            data = self._post("decode", request, input.potential_sequence_of_concern)
        except ESMProteinError as e:
            return e

        return ESMProtein(
            sequence=data["outputs"]["sequence"],
            secondary_structure=data["outputs"]["secondary_structure"],
            sasa=data["outputs"]["sasa"],
            function_annotations=_list_to_function_annotations(
                data["outputs"]["function"]
            ),
            coordinates=maybe_tensor(
                data["outputs"]["coordinates"], convert_none_to_nan=True
            ),
            plddt=maybe_tensor(data["outputs"]["plddt"]),
            ptm=maybe_tensor(data["outputs"]["ptm"]),
        )

    @retry_decorator
    def logits(
        self,
        input: ESMProteinTensor,
        config: LogitsConfig = LogitsConfig(),
        return_bytes: bool = True,
    ) -> LogitsOutput | ESMProteinError:
        _validate_protein_tensor_input(input)

        # Note: using raw model forwards is discouraged because of the byte size
        # of the logits.
        # Please use forward_and_sample instead.
        req = {}
        req["sequence"] = maybe_list(input.sequence)
        req["structure"] = maybe_list(input.structure)
        req["secondary_structure"] = maybe_list(input.secondary_structure)
        req["sasa"] = maybe_list(input.sasa)
        req["function"] = maybe_list(input.function)
        req["coordinates"] = maybe_list(input.coordinates, convert_nan_to_none=True)
        req["residue_annotation"] = maybe_list(input.residue_annotations)

        logits_config = {
            "sequence": config.sequence,
            "structure": config.structure,
            "secondary_structure": config.secondary_structure,
            "sasa": config.sasa,
            "function": config.function,
            "residue_annotations": config.residue_annotations,
            "return_embeddings": config.return_embeddings,
            "return_hidden_states": config.return_hidden_states,
            "ith_hidden_layer": config.ith_hidden_layer,
        }

        request = {"model": self.model, "inputs": req, "logits_config": logits_config}
        try:
            data = self._post(
                "logits",
                request,
                input.potential_sequence_of_concern,
                return_bytes=return_bytes,
            )
        except ESMProteinError as e:
            return e

        def _maybe_logits(track: str):
            if "logits" in data and track in data["logits"]:
                return maybe_tensor(data["logits"][track])
            return None

        def _maybe_b64_decode(obj):
            return (
                deserialize_tensors(base64.b64decode(obj))
                if return_bytes and obj is not None
                else obj
            )

        logits = _maybe_b64_decode(data["logits"])
        data["logits"] = dict(logits) if logits is not None else logits
        data["embeddings"] = _maybe_b64_decode(data["embeddings"])
        data["hidden_states"] = _maybe_b64_decode(data["hidden_states"])

        output = LogitsOutput(
            logits=ForwardTrackData(
                sequence=_maybe_logits("sequence"),
                structure=_maybe_logits("structure"),
                secondary_structure=_maybe_logits("secondary_structure"),
                sasa=_maybe_logits("sasa"),
                function=_maybe_logits("function"),
            ),
            embeddings=maybe_tensor(data["embeddings"]),
            residue_annotation_logits=_maybe_logits("residue_annotation"),
            hidden_states=maybe_tensor(data["hidden_states"]),
        )

        return output

    def _post(
        self,
        endpoint,
        request,
        potential_sequence_of_concern,
        return_bytes: bool = False,
    ):
        request["potential_sequence_of_concern"] = potential_sequence_of_concern
        headers = dict(self.headers)
        if return_bytes:
            headers["return-bytes"] = "true"
        response = requests.post(
            urljoin(self.url, f"/api/v1/{endpoint}"),
            json=request,
            headers=headers,
            timeout=self.request_timeout,
        )

        if not response.ok:
            raise ESMProteinError(
                error_code=response.status_code,
                error_msg=f"Failure in {endpoint}: {response.text}",
            )

        data = response.json()
        # Nextjs puts outputs dict under "data" key.
        # Lift it up for easier downstream processing.
        if "outputs" not in data and "data" in data:
            data = data["data"]

        return data

    @property
    def raw_model(self):
        raise NotImplementedError(
            f"Can not get underlying remote model {self.model} from a Forge client."
        )
