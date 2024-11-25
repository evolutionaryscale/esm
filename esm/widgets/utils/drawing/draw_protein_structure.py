import py3Dmol
from IPython.display import clear_output
from ipywidgets import widgets

from esm.utils.structure.protein_chain import ProteinChain


def draw_protein_structure(
    output: widgets.Output,
    protein_chain: ProteinChain,
    highlighted_ranges: list[tuple[int, int, str]] = [],
):
    pdb_str = protein_chain.to_pdb_string()
    with output:
        clear_output(wait=True)
        view = py3Dmol.view(width=500, height=500)
        view.addModel(pdb_str, "pdb")
        view.setStyle({"cartoon": {"color": "gray"}})

        for start, end, color in highlighted_ranges:
            view.setStyle(
                {"resi": str(start) + "-" + str(end)}, {"cartoon": {"color": color}}
            )

        view.zoomTo()
        view.show()
