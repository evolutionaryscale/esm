- [Installation ](#installation)
- [ESM C](#esm-c)
  - [ESM C 300M and 600M via GitHub](#esm-c-github)
  - [ESM C via Forge API for Free Non-Commercial Use](#esm-c-forge)
  - [ESM C via SageMaker for Commercial Use](#esm-c-sagemaker)
  - [ESM C Example Usage](#esmc-example)
- [ESM 3](#esm3)
  - [Quickstart for ESM3-open](#esm3-quickstart)
  - [Forge: Access to larger ESM3 models](#esm3-forge)
  - [ESM 3 Example Usage](#esm3-example)
- [Responsible Development ](#responsible-development)
- [Licenses](#licenses)


## Installation <a name="installation"></a>

To get started with ESM, install the python library using pip:

```bash
pip install esm
```

## ESM C <a name="esm-c"></a>
[ESM Cambrian](https://www.evolutionaryscale.ai/blog/esm-cambrian) is a parallel model family to our flagship ESM3 generative models. While ESM3 focuses on controllable generation of proteins for therapeutic and many other applications, ESM C focuses on creating representations of the underlying biology of proteins.

ESM C comes with major performance benefits over ESM2. The 300M parameter ESM C delivers similar performance to ESM2 650M with dramatically reduced memory requirements and faster inference. The 600M parameter ESM C rivals the 3B parameter ESM2 and approaches the capabilities of the 15B model, delivering frontier performance with far greater efficiency. The 6B parameter ESM C sets a new benchmark, outperforming the best ESM2 models by a wide margin.

ESM C models are available immediately for academic and commercial use under a new license structure designed to promote openness and enable scientists and builders. You can find the high level take-away of the license structure in the [Licenses](#licenses) section of this page, and complete license details in [LICENSE.md](LICENSE.md).

You can use the following guides to start using ESM C models today, either running the model locally, [the Forge API](https://forge.evolutionaryscale.ai/) and [AWS SageMaker](https://aws.amazon.com/marketplace/seller-profile?id=seller-iw2nbscescndm).

### ESM C Local Models via GitHub <a name="esm-c-github"></a>
The code and weights for the ESM C 300M model are available under the Cambrian Open [license agreement](#licenses). The weights for the ESM C 600M model are available under the Cambrian Non-Commercial [license agreement](#licenses).

When running the code below, a pytorch model will be instantiated locally on your machine, with the weights downloaded from the [HuggingFace hub](https://huggingface.co/EvolutionaryScale).
```py
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig

protein = ESMProtein(sequence="AAAAA")
client = ESMC.from_pretrained("esmc_300m").to("cuda") # or "cpu"
protein_tensor = client.encode(protein)
logits_output = client.logits(
   protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
)
print(logits_output.logits, logits_output.embeddings)
```

### ESM C via Forge API for Free Non-Commercial Use  <a name="esm-c-forge"></a>

The ESM C model family, including ESMC 6B, are accessible via EvolutionaryScale Forge for free [non-commercial use](#licenses).
Apply for access and copy the API token from the console by first visiting https://forge.evolutionaryscale.ai.

With the code below, a local python client talks to the model inference server hosted by EvolutionaryScale.

```py
from esm.sdk.forge import ESM3ForgeInferenceClient
from esm.sdk.api import ESMProtein, LogitsConfig

# Apply for forge access and get an access token
forge_client = ESM3ForgeInferenceClient(model="esmc-6b-2024-12", url="https://forge.evolutionaryscale.ai", token="<your forge token>")
protein_tensor = forge_client.encode(protein)
logits_output = forge_client.logits(
   protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
)
print(logits_output.logits, logits_output.embeddings)
```

Remember to replace `<your forge token>` with your actual Forge access token.

### ESM C via SageMaker for Commercial Use  <a name="esm-c-sagemaker"></a>

ESM C models are also available on Amazon SageMaker under the Cambrian Inference Clickthrough License Agreement.
Under this license agreement models are available for broad commercial use to commercial entities.

You will need an admin AWS access to an AWS account to follow these instructions. To deploy, first we need to deploy the AWS package:

1. Find the ESM C model version you want to subscribe to. All of our offerings are visible [here](https://aws.amazon.com/marketplace/seller-profile?id=seller-iw2nbscescndm).
2. Click the name of the model version you are interested in, review pricing information and the end user license agreement (EULA), then click "Continue to Subscribe".
3. Once you have subscribed, you should be able to see our model under your [marketplace subscriptions](https://us-east-1.console.aws.amazon.com/marketplace/home#/subscriptions).
4. Click the product name and then from the "Actions" dropdown select "Configure".
5. You will next see the "Configure and Launch" UI. There are multiple deployment paths - we recommend using "AWS CloudFormation".
6. The default value for "Service Access" may or may not work. We recommend clicking "Create and use a new service role".
7. Click "Launch CloudFormation Template".  This takes 15 to 25 minutes depending on model size.
8. On the "Quick create stack" page, ensure the stack name and endpoint names are not already used. You can check existing stack names [here](https://console.aws.amazon.com/cloudformation/home/stacks) and existing endpoint names [here](https://us-east-1.console.aws.amazon.com/sagemaker/home?region=us-east-1#/endpoints).

The Sagemaker deployment of the model now lives on a dedicated GPU instance inside your AWS environment, and will be billed directly to your AWS account.
Make sure to remember to shut down the instance after you stop using it. Find the CloudFormation stack you created [here](https://us-east-1.console.aws.amazon.com/cloudformation/home), select it, and then click "Delete" to clean up all resources.

After creating the endpoint, you can create a sagemaker client and use it the same way as a forge client. They share the same API.
The local python client talks to the Sagemaker endpoint you just deployed, which runs on an instance with a GPU to run model inference.


Ensure that the code below runs in an environment that has AWS credentials available for the account which provisioned SageMaker resources.  Learn more about general AWS credential options [here](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-authentication.html#cli-chap-authentication-precedence).

```py
from esm.sdk.sagemaker import ESM3SageMakerClient
from esm.sdk.api import ESMProtein, LogitsConfig

sagemaker_client = ESM3SageMakerClient(
   # E.g. "Endpoint-ESMC-6B-1"
   endpoint_name=SAGE_ENDPOINT_NAME,
   # E.g. "esmc-6b-2024-12". Same model names as in Forge.
   model=MODEL_NAME,
)

protein = ESMProtein(sequence="AAAAA")
protein_tensor = sagemaker_client.encode(protein)
logits_output = sagemaker_client.logits(
   protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
)
print(logits_output.logits, logits_output.embeddings)
```

### ESM C Example Usage
 <a name="esmc-example"></a>
Look at [esmc_examples.py](./examples/esmc_examples.py) for the standard usage (extracting embeddings and model amino acid prediction).

More coming soon.

## ESM 3  <a name="esm3"></a>

[ESM3](https://www.evolutionaryscale.ai/papers/esm3-simulating-500-million-years-of-evolution-with-a-language-model) is a frontier generative model for biology, able to jointly reason across three fundamental biological properties of proteins: sequence, structure, and function. These three data modalities are represented as tracks of discrete tokens at the input and output of ESM3. You can present the model with a combination of partial inputs across the tracks, and ESM3 will provide output predictions for all the tracks.

ESM3 is a _generative_ masked language model. You can prompt it with partial sequence, structure, and function keywords, and iteratively sample masked positions until all positions are unmasked. This iterative sampling is what the `.generate()` function does.

<!--![ESM3 Diagram](_assets/esm3_diagram.png)-->
<img src="_assets/esm3_diagram.png" alt="ESM3 Diagram" width="400" />

The ESM3 architecture is highly scalable due to its transformer backbone and all-to-all reasoning over discrete token sequences. At its largest scale, ESM3 was trained with 1.07e24 FLOPs on 2.78 billion proteins and 771 billion unique tokens, and has 98 billion parameters.
Learn more by reading the [blog post](https://www.evolutionaryscale.ai/blog/esm3-release) and [the pre-print (Hayes et al., 2024)](https://www.evolutionaryscale.ai/papers/esm3-simulating-500-million-years-of-evolution-with-a-language-model).

Here we present `esm3-open-small`. With 1.4B parameters it is the smallest and fastest model in the family.
ESM3-open is available under the [Cambrian non-commercial license agreement](https://www.evolutionaryscale.ai/policies/cambrian-non-commercial-license-agreement), as outlined in `LICENSE.md` (note: updated with ESM C release).
Visit our [Discussions page](https://github.com/evolutionaryscale/esm/discussions) to get in touch, provide feedback, ask questions or share your experience with ESM3!

### Quickstart for ESM3-open <a name="esm3-quickstart"></a>

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

### EvolutionaryScale Forge: Access to larger ESM3 models
 <a name="esm3-forge"></a>

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

### ESM3 Example Usage
 <a name="esm3-example"></a>
Let's explore some more advanced prompting with the help of our [notebooks and scripts](examples/).

`generate.ipynb` will walk through two prompting examples (scaffolding and secondary structure editing) using the open model:
[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/evolutionaryscale/esm/blob/main/examples/generate.ipynb)

`gfp_design.ipynb` will walk through the more complex generation procedure we used to design esmGFP:
[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/evolutionaryscale/esm/blob/main/examples/gfp_design.ipynb)

We also provide example scripts that show common workflows under `examples/`:

- [local_generate.py](./examples/local_generate.py) shows how simple and elegant common tasks are: it shows folding, inverse folding and chain of thought generation, all by calling just `model.generate()` for iterative decoding.
- [seqfun_struct.py](./examples/seqfun_struct.py) shows direct use of the model as a standard pytorch model with a simple model `forward` call.

## Responsible Development <a name="responsible-development"></a>

EvolutionaryScale is a public benefit company. Our mission is to develop artificial intelligence to understand biology for the benefit of human health and society, through partnership with the scientific community, and open, safe, and responsible research. Inspired by the history of our field as well as [new principles and recommendations](https://responsiblebiodesign.ai/), we have created a Responsible Development Framework to guide our work towards our mission with transparency and clarity.

The core tenets of our framework are

- We will communicate the benefits and risks of our research
- We will proactively and rigorously evaluate the risk of our models before public deployment
- We will adopt risk mitigation strategies and precautionary guardrails
- We will work with stakeholders in government, policy, and civil society to keep them informed

With this in mind, we have performed a variety of mitigations for `esm3-sm-open-v1`, detailed in our [paper](https://www.evolutionaryscale.ai/papers/esm3-simulating-500-million-years-of-evolution-with-a-language-model)

## Licenses  <a name="licenses"></a>
The code and model weights of ESM3 and ESM C are available under a mixture of non-commercial and permissive commercial licenses.
This summary provides a high-level overview. For complete license details, see [LICENSE.md](./LICENSE.md).

### How can I access the models and which licenses apply?

The models can be accessed in three different ways, each with its own licensing terms.

1. **Code and weights** via GitHub and HuggingFace are available under either a [non-commercial](https://www.evolutionaryscale.ai/policies/cambrian-non-commercial-license-agreement) (ESM C 600M, ESM3-small-open) or an [open license](https://www.evolutionaryscale.ai/policies/cambrian-open-license-agreement) (codebase, ESM C 300M).
    1. **Building with ESM encouraged**: You can use embeddings, model predictions, fine-tune the models and use components of both the models and code. We strongly encourage anyone to build on ESM C and ESM3! Just remember to maintain the same license terms and release under the ESM name.
2. **Free non-commercial inference API** via Forge. All models are available this way, with free credits granted to students and researchers. We want to enable academics under [non-commercial Terms of Use](https://www.evolutionaryscale.ai/policies/terms-of-use), which mirrors the non-commercial license.
3. **Paid commercial Inference API** for commercial use via SageMaker (Forge coming soon). All ESM C models are available this way to commercial entities for commercial use under a [clickthrough license agreement](https://www.evolutionaryscale.ai/policies/cambrian-inference-clickthrough-license-agreement) with few restrictions. 
    1. In broad strokes: standard commercial use like developing molecules and developing downstream ML models and methods with the model is allowed, while training competing models on the API outputs is not.
    2. Note: For ESM3 commercial use, reach out to [bd@evolutionaryscale.ai](mailto:bd@evolutionaryscale.ai) 

### What changed with the release of ESM C?

We introduced a [clickthrough license agreement](https://www.evolutionaryscale.ai/policies/cambrian-inference-clickthrough-license-agreement) to enable frictionless commercial use of ESM C. 

We introduced the new [Cambrian Open License](https://www.evolutionaryscale.ai/policies/cambrian-open-license-agreement) for ESM C 300M, and at the same time moved all code in the [`esm` repo](https://github.com/evolutionaryscale/esm) under that permissive license.

The [Cambrian non-commercial license](https://www.evolutionaryscale.ai/policies/cambrian-non-commercial-license-agreement) is largely based on the original [ESM3 Community License Agreement](https://www.evolutionaryscale.ai/policies/community-license-agreement), but removed the clause that restricted drug development, added the naming requirement, and extended the meaning of “Derivative Work” to allow training on model outputs. Just remember to release models and methods built on ESM under the same license.
These changes are meant to remove potential gray areas and points of friction for researchers building with ESM.

Finally, The ESM3-open-small model has been moved under the Cambrian non-commercial license.
