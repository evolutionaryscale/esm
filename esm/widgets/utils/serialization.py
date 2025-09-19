import base64
import json
from io import StringIO
from typing import Literal

from ipywidgets import widgets

from esm.sdk.api import ESMProtein


def protein_to_pdb_buffer(protein: ESMProtein) -> bytes:
    pdb_buffer = StringIO()
    protein.to_pdb(pdb_buffer)
    pdb_buffer.seek(0)
    return pdb_buffer.read().encode()


def create_download_button_from_buffer(
    buffer: bytes,
    filename: str,
    description: str = "Download",
    type: Literal["json", "bytes"] = "bytes",
) -> widgets.HTML:
    b64 = base64.b64encode(buffer).decode()
    if type == "json":
        payload = f"data:text/json;base64,{b64}"
    elif type == "bytes":
        payload = f"data:application/octet-stream;base64,{b64}"
    html_buttons = f"""
    <html>
    <head>
    <meta name="viewport" content="width=device-width, initial-scale=1">
    </head>
    <body>
    <a download="{filename}" href="{payload}" download>
    <button class="p-Widget jupyter-widgets jupyter-button widget-button">{description}</button>
    </a>
    </body>
    </html>
    """
    download_link = widgets.HTML(html_buttons)
    return download_link


def create_download_results_button(
    protein_list: list[ESMProtein], filename: str
) -> widgets.HTML:
    serialized_proteins = [serialize_protein(p) for p in protein_list]
    serialized_data = json.dumps(serialized_proteins, indent=4)
    return create_download_button_from_buffer(
        buffer=serialized_data.encode(),
        filename=filename,
        type="json",
        description="Download As JSON",
    )


def serialize_protein(protein: ESMProtein) -> str:
    protein_dict = {
        "sequence": protein.sequence,
        "coordinates": protein.coordinates.tolist()
        if protein.coordinates is not None
        else None,
        "secondary_structure": protein.secondary_structure,
        "sasa": protein.sasa,
        "function_annotations": [
            (annotation.label, annotation.start, annotation.end)
            for annotation in protein.function_annotations
        ]
        if protein.function_annotations is not None
        else None,
        "plddt": protein.plddt.tolist() if protein.plddt is not None else None,
        "ptm": protein.ptm.tolist() if protein.ptm is not None else None,
    }
    return json.dumps(protein_dict, indent=4)
