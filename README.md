# ESM3

[ESM3](https://www.evolutionaryscale.ai/papers/esm3-simulating-500-million-years-of-evolution-with-a-language-model) is a frontier generative model for biology, able to jointly reason across three fundamental biological properties of proteins: sequence, structure, and function. These three data modalities are represented as tracks of discrete tokens at the input and output of ESM3. You can present the model with a combination of partial inputs across the tracks, and ESM3 will provide output predictions for all the tracks.

ESM3 is a _generative_ masked language model. You can prompt it with partial sequence, structure, and function keywords, and iteratively sample masked positions until all positions are unmasked. This iterative sampling is what the `.generate()` function does.

<!--![ESM3 Diagram](_assets/esm3_diagram.png)-->
<img src="_assets/esm3_diagram.png" alt="ESM3 Diagram" width="400" />

The ESM3 architecture is highly scalable due to its transformer backbone and all-to-all reasoning over discrete token sequences. At its largest scale, ESM3 was trained with 1.07e24 FLOPs on 2.78 billion proteins and 771 billion unique tokens, and has 98 billion parameters.
Learn more by reading the [blog post](https://www.evolutionaryscale.ai/blog/esm3-release) and [the pre-print (Hayes et al., 2024)](https://www.evolutionaryscale.ai/papers/esm3-simulating-500-million-years-of-evolution-with-a-language-model).

Here we present `esm3-open-small`. With 1.4B parameters it is the smallest and fastest model in the family.
ESM3-open is available under a [non-commercial license](https://www.evolutionaryscale.ai/policies/community-license-agreement), reproduced under `LICENSE.md`.
Visit our [Discussions page](https://github.com/evolutionaryscale/esm/discussions) to get in touch, provide feedback, ask questions or share your experience with ESM3!

## Quickstart for ESM3-open

```
pip install esm
```

In order to download the weights, we require users to accept our non-commercial license.
The weights are stored on HuggingFace Hub under [HuggingFace/EvolutionaryScale/esm3](https://huggingface.co/EvolutionaryScale/esm3).
Please create an account and accept the license.

```py
from huggingface_hub import login
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig

# Will instruct you how to get an API key from huggingface hub, make one with "Read" permission.
login()

# This will download the model weights and instantiate the model on your machine.
model: ESM3InferenceClient = ESM3.from_pretrained("esm3-open").to("cuda") # or "cpu"

# Generate a completion for a partial Carbonic Anhydrase (2vvb)
prompt = "___________________________________________________DQATSLRILNNGHAFNVEFDDSQDKAVLKGGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVVDVLDSIKTKGKSADFTNFDPRGLLPESLDYWTYPGSLTTPP___________________________________________________________"
protein = ESMProtein(sequence=prompt)
# Generate the sequence, then the structure. This will iteratively unmask the sequence track.
protein = model.generate(protein, GenerationConfig(track="sequence", num_steps=8, temperature=0.7))
# We can show the predicted structure for the generated sequence.
protein = model.generate(protein, GenerationConfig(track="structure", num_steps=8))
protein.to_pdb("./generation.pdb")
# Then we can do a round trip design by inverse folding the sequence and recomputing the structure
protein.sequence = None
protein = model.generate(protein, GenerationConfig(track="sequence", num_steps=8))
protein.coordinates = None
protein = model.generate(protein, GenerationConfig(track="structure", num_steps=8))
protein.to_pdb("./round_tripped.pdb")
```

Congratulations, you just generated your first proteins with ESM3!
Let's explore some more advanced prompting with the help of our [notebooks and scripts](examples/).

`generate.ipynb` will walk through two prompting examples (scaffolding and secondary structure editing) using the open model:
[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/evolutionaryscale/esm/blob/main/examples/generate.ipynb)

`gfp_design.ipynb` will walk through the more complex generation procedure we used to design esmGFP:
[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/evolutionaryscale/esm/blob/main/examples/gfp_design.ipynb)

We also provide example scripts that show common workflows under `examples/`:

- [local_generate.py](./examples/local_generate.py) shows how simple and elegant common tasks are: it shows folding, inverse folding and chain of thought generation, all by calling just `model.generate()` for iterative decoding.
- [seqfun_struct.py](./examples/seqfun_struct.py) shows direct use of the model as a standard pytorch model with a simple model `forward` call.

## Forge: Access to larger ESM3 models

You can apply for beta access to the full family of larger and higher capability ESM3 models at [EvolutionaryScale Forge](https://forge.evolutionaryscale.ai).

We encourage users to interact with the Forge API through the python `esm` library instead of the command line.
The python interface enables you to interactively load proteins, build prompts, and inspect generated proteins
with the `ESMProtein` and config classes used to interact with the local model.

In any example script try to replace a local `ESM3` model with a Forge API client:

```py
# Instead of loading the model locally on your machine:
model: ESM3InferenceClient = ESM3.from_pretrained("esm3_sm_open_v1").to("cuda") # or "cpu"
# just replace the line with this:
model: ESM3InferenceClient = esm.sdk.client("esm3-medium-2024-08", token="<your forge token>")
# and now you're interfacing with the model running on our remote servers.
...
```

and the exact same code will work.
This enables a seamless transition from smaller and faster models, to our large 98B protein language models for protein design work.

## Responsible Development

EvolutionaryScale is a public benefit company. Our mission is to develop artificial intelligence to understand biology for the benefit of human health and society, through partnership with the scientific community, and open, safe, and responsible research. Inspired by the history of our field as well as [new principles and recommendations](https://responsiblebiodesign.ai/), we have created a Responsible Development Framework to guide our work towards our mission with transparency and clarity.

The core tenets of our framework are

- We will communicate the benefits and risks of our research
- We will proactively and rigorously evaluate the risk of our models before public deployment
- We will adopt risk mitigation strategies and precautionary guardrails
- We will work with stakeholders in government, policy, and civil society to keep them informed

With this in mind, we have performed a variety of mitigations for `esm3-sm-open-v1`, detailed in our [paper](https://www.evolutionaryscale.ai/papers/esm3-simulating-500-million-years-of-evolution-with-a-language-model)

## License

**The Big Picture:**

1. The EvolutionaryScale AI Model is **only** available under this Community License Agreement for **non-commercial use** by **individuals** or **non-commercial organizations** (including universities, non-profit organizations and research institutes, educational and government bodies).

2. You **may not** use the EvolutionaryScale AI Model or any derivative works of the EvolutionaryScale AI Model or its outputs:

   1. in connection with **any commercial activities**, for example, any activities **by, on behalf of or for a commercial entity** or to develop **any product or service** such as hosting the AI Model behind an API; or

   2. without attribution to EvolutionaryScale and this Community License Agreement; or

   3. to **train** a AI-powered third party model **similar to EvolutionaryScale’s AI Model**, even for non-commercial usage. You may, however, create **Derivative Works** of ESM3, for example by finetuning or adding model layers.

3. You **can publish, share and adapt** the EvolutionaryScale AI Model and its outputs for **non-commercial purposes** in accordance with the Community License Agreement, including a **non-commercial restriction** on the adapted model.

Please read our non-commercial [Community License Agreement](https://www.evolutionaryscale.ai/policies/community-license-agreement) reproduced under [./LICENSE.md](LICENSE.md) before using ESM3.
