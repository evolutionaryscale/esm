import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor

# NOT CURRENTLY USED


class SwiGLU(nn.Module):
    def __init__(self) -> None:
        super().__init__()

    def forward(self, x: Tensor) -> Tensor:
        x1, x2 = x.chunk(2, dim=-1)
        hidden = F.silu(x1) * x2
        return hidden


class FFN(nn.Module):
    def __init__(self, in_proj, activation, out_proj) -> None:
        super().__init__()
        self.in_proj = in_proj
        self.activation = activation
        self.out_proj = out_proj

    def forward(self, x: Tensor) -> Tensor:
        x = self.in_proj(x)
        x = self.activation(x)
        x = self.out_proj(x)
        return x
