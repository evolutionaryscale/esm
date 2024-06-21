from esm.models.esm3 import ESM3
from esm.sdk.api import (
    ESM3InferenceClient,
    ESMProtein,
    GenerationConfig,
    SamplingConfig,
    SamplingTrackConfig,
)
from esm.utils.structure.protein_chain import ProteinChain
from esm.utils.types import FunctionAnnotation


def get_sample_protein() -> ESMProtein:
    protein = ProteinChain.from_rcsb("1utn")
    protein = ESMProtein.from_protein_chain(protein)
    protein.function_annotations = [
        # Peptidase S1A, chymotrypsin family: https://www.ebi.ac.uk/interpro/structure/PDB/1utn/
        FunctionAnnotation(label="peptidase", start=100, end=114),
        FunctionAnnotation(label="chymotrypsin", start=190, end=202),
    ]
    return protein


def main(client: ESM3InferenceClient):
    # Single step decoding
    protein = get_sample_protein()
    protein.function_annotations = None
    protein = client.encode(protein)
    single_step_protein = client.forward_and_sample(
        protein,
        SamplingConfig(structure=SamplingTrackConfig(topk_logprobs=2)),
    )
    single_step_protein.protein_tensor.sequence = protein.sequence
    single_step_protein = client.decode(single_step_protein.protein_tensor)

    # Folding
    protein = get_sample_protein()
    sequence_length = len(protein.sequence)  # type: ignore
    num_steps = int(sequence_length / 16)
    protein.coordinates = None
    protein.function_annotations = None
    protein.sasa = None
    folded_protein = client.generate(
        protein,
        GenerationConfig(track="structure", schedule="cosine", num_steps=num_steps),
    )
    folded_protein.to_pdb("./sample_folded.pdb")

    # Inverse Folding
    protein = get_sample_protein()
    protein.sequence = None
    protein.sasa = None
    protein.function_annotations = None
    inv_folded_protein = client.generate(
        protein,
        GenerationConfig(track="sequence", schedule="cosine", num_steps=num_steps),
    )
    inv_folded_protein.to_pdb("./sample_inv_folded.pdb")

    # Chain of Thought (Function -> Secondary Structure -> Structure -> Sequence)
    cot_protein = get_sample_protein()
    cot_protein.sequence = "_" * len(cot_protein.sequence)  # type: ignore
    cot_protein.coordinates = None
    cot_protein.sasa = None
    cot_protein_tensor = client.encode(cot_protein)
    for cot_track in ["secondary_structure", "structure", "sequence"]:
        cot_protein_tensor = client.generate(
            cot_protein_tensor,
            GenerationConfig(track=cot_track, schedule="cosine", num_steps=10),
        )
    cot_protein = client.decode(cot_protein_tensor)
    cot_protein.to_pdb("./sample_cot.pdb")


if __name__ == "__main__":
    main(ESM3.from_pretrained("esm3_sm_open_v1"))
