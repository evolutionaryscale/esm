from tokenizers import Tokenizer
from tokenizers.models import BPE
from tokenizers.processors import TemplateProcessing
from transformers import PreTrainedTokenizerFast

from esm.tokenization.tokenizer_base import EsmTokenizerBase
from esm.utils.constants import esm3 as C


class EsmSequenceTokenizer(PreTrainedTokenizerFast, EsmTokenizerBase):
    """
    Constructs an ESM tokenizer.
    """

    model_input_names = ["sequence_tokens", "attention_mask"]

    def __init__(
        self,
        unk_token="<unk>",
        cls_token="<cls>",
        pad_token="<pad>",
        mask_token="<mask>",
        eos_token="<eos>",
        chain_break_token="|",
        **kwargs,
    ):
        all_tokens = C.SEQUENCE_VOCAB
        token_to_id = {tok: ind for ind, tok in enumerate(all_tokens)}

        # a character-level tokenizer is the same as BPE with no token merges
        bpe = BPE(token_to_id, merges=[], unk_token=unk_token)
        tokenizer = Tokenizer(bpe)
        special_tokens = [
            cls_token,
            pad_token,
            mask_token,
            eos_token,
            chain_break_token,
        ]
        self.cb_token = chain_break_token
        additional_special_tokens = [chain_break_token]

        tokenizer.add_special_tokens(special_tokens)

        # This is where we configure the automatic addition of special tokens when we call
        # tokenizer(text, add_special_tokens=True). Note that you can also configure how two
        # sequences are merged if you want.
        tokenizer.post_processor = TemplateProcessing(  # type: ignore
            single="<cls> $A <eos>",
            special_tokens=[
                ("<cls>", tokenizer.token_to_id("<cls>")),
                ("<eos>", tokenizer.token_to_id("<eos>")),
            ],
        )
        super().__init__(
            tokenizer_object=tokenizer,
            unk_token=unk_token,
            cls_token=cls_token,
            pad_token=pad_token,
            mask_token=mask_token,
            eos_token=eos_token,
            additional_special_tokens=additional_special_tokens,
            **kwargs,
        )

    # These are a footgun, we never use the `bos` token anywhere so we're just overriding it here.
    @property
    def bos_token(self):
        return self.cls_token

    @property
    def bos_token_id(self):
        return self.cls_token_id

    @property
    def chain_break_token(self):
        return self.cb_token

    @property
    def chain_break_token_id(self):
        return self.convert_tokens_to_ids(self.chain_break_token)

    @property
    def all_token_ids(self):
        return list(range(self.vocab_size))

    @property
    def special_token_ids(self):
        return self.all_special_ids
