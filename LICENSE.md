# License Information

Here are the different licenses that govern access to the ESM codebase and the models inclusive of weights:

| License                                                    | Component                                                           |
|------------------------------------------------------------|-------------------------------------------------------------------|
| [Cambrian Open License Agreement](https://www.evolutionaryscale.ai/policies/cambrian-open-license-agreement)        | Code on GitHub (excluding model weights)                         |
| [Cambrian Open License Agreement](https://www.evolutionaryscale.ai/policies/cambrian-open-license-agreement)    | ESM C 300M (incl weights)                         |
| [Cambrian Non-Commercial License Agreement](https://www.evolutionaryscale.ai/policies/cambrian-non-commercial-license-agreement)    | ESM-3 Open Model (incl weights)        |
| [Cambrian Non-Commercial License Agreement](https://www.evolutionaryscale.ai/policies/cambrian-non-commercial-license-agreement)   | ESM C 600M (incl weights)                            |
| Governed by API Agreements (See Below)             | API access to all models, including API-only models (ESM3 family, ESM C 6B) |
| [Forge API Terms of Use](https://www.evolutionaryscale.ai/policies/terms-of-use)                                          | Free non-commercial API access via Forge to all models including API-only models (ESM3 family, ESM C 6B) |
| [Cambrian Inference Clickthrough License Agreement](https://www.evolutionaryscale.ai/policies/cambrian-inference-clickthrough-license-agreement)    | Commercial Inference via SageMaker for all ESM C models |


# How can I access the models and which licenses apply?

The models can be accessed in three different ways, each with its own licensing terms.

1. **Code and weights** via GitHub and HuggingFace are available under either a [non-commercial](https://www.evolutionaryscale.ai/policies/cambrian-non-commercial-license-agreement) (ESM C 600M, ESM3-small-open) or an [open license](https://www.evolutionaryscale.ai/policies/cambrian-open-license-agreement) (codebase, ESM C 300M).
    * **Building with ESM encouraged**: You can use embeddings, model predictions, fine-tune the models and use components of both the models and code. We strongly encourage anyone to build on ESM C and ESM3! Just remember to maintain the same license terms and release under the ESM name.
2. **Free non-commercial inference API** via Forge. All models are available this way, with free credits granted to students and researchers. We want to enable academics under [non-commercial Terms of Use](https://www.evolutionaryscale.ai/policies/terms-of-use), which mirrors the non-commercial license.
3. **Paid commercial Inference API** for commercial use via SageMaker (Forge coming soon). All ESM C models are available this way to commercial entities for commercial use under a [clickthrough license agreement](https://www.evolutionaryscale.ai/policies/cambrian-inference-clickthrough-license-agreement) with few restrictions.
    * In broad strokes: standard commercial use like developing molecules and developing downstream ML models and methods with the model is allowed, while training competing models on the API outputs is not.
    * Note: For ESM3 commercial use, reach out to [bd@evolutionaryscale.ai](mailto:bd@evolutionaryscale.ai)

### What changed with the release of ESM C?

We introduced a [clickthrough license agreement](https://www.evolutionaryscale.ai/policies/cambrian-inference-clickthrough-license-agreement) to enable frictionless commercial use of ESM C.

We introduced the new [Cambrian Open License](https://www.evolutionaryscale.ai/policies/cambrian-open-license-agreement) for ESM C 300M, and at the same time moved all code in the [`esm` repo](https://github.com/evolutionaryscale/esm) under that permissive license.

The [Cambrian non-commercial license](https://www.evolutionaryscale.ai/policies/cambrian-non-commercial-license-agreement) is largely based on the original [ESM3 Community License Agreement](https://www.evolutionaryscale.ai/policies/community-license-agreement), but removed the clause that restricted drug development, added the naming requirement, and extended the meaning of “Derivative Work” to allow training on model outputs. Just remember to release models and methods built on ESM under the same license.
These changes are meant to remove potential gray areas and points of friction for researchers building with ESM.

Finally, The ESM3-open-small model has been moved under the [Cambrian non-commercial license](https://www.evolutionaryscale.ai/policies/cambrian-non-commercial-license-agreement).
