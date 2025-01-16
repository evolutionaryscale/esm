We welcome community contributions to help make this package better!

## Contributing

_to be written_

## Testing

Improving test coverage and automating the process is high priority but not done yet. For now, to test:

0. Ensure you are testing in a clean environment - we use micromamba for this

```bash
micromamba create -n esm
micromamba activate esm
micromamba install -c conda-forge python=3.10

# in root level of repo
pip install -e .
pip install examples/requirements.txt
python -c 'from huggingface_hub import login; login()'


```

```bash

# Ensure package can correctly interact with forge.
# This will require an API key from forge.evolutionaryscale.ai
ESM_API_KEY=$ESM_API_KEY PYTHONPATH='.' python examples/api/esm3.py



```

1. Ensure the following scripts run without errors. Most have a pip install command installing the published `esm` package - comment this out so your release candidate version is tested and not the already published version.

```bash
ESM_API_KEY=$ESM_API_KEY PYTHONPATH='.' python examples/api/esm3.py
pip install treon
treon examples/tutorials/1_esmprotein.ipynb
treon examples/tutorials/4_gfp_design.ipynb

# requires a GPU
treon examples/tutorials/2_generate.ipynb
python examples/api/esm3.py
python examples/local/raw_forwards.py


```

`examples/tutorials/1_esmprotein.ipynb` works. Remember to skip running the first cell - it will reinstall stock esm instead of your deployed version.

3. Ensure `examples/tutorials/3_generate.ipynb` works. Note this notebook will require a node with a GPU that can fit ESM3 small open.

4. Ensure
