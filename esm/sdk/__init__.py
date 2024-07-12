import os

from esm.sdk.forge import ESM3ForgeInferenceClient


def client(
    model="esm3-sm-open-v1",
    url="https://forge.evolutionaryscale.ai",
    token=os.environ.get("ESM_API_KEY", ""),
):
    return ESM3ForgeInferenceClient(model, url, token)
