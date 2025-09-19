import numpy as np
import torch
import torch.distributed as dist
import torch.nn as nn
import torch.nn.functional as F


class EMACodebook(nn.Module):
    def __init__(
        self,
        n_codes,
        embedding_dim,
        no_random_restart=True,
        restart_thres=1.0,
        ema_decay=0.99,
    ):
        super().__init__()
        self.register_buffer("embeddings", torch.randn(n_codes, embedding_dim))
        self.register_buffer("N", torch.zeros(n_codes))
        self.register_buffer("z_avg", self.embeddings.data.clone())  # pyright: ignore[reportCallIssue]

        self.n_codes = n_codes
        self.embedding_dim = embedding_dim
        self._need_init = True
        self.no_random_restart = no_random_restart
        self.restart_thres = restart_thres
        self.freeze_codebook = False
        self.ema_decay = ema_decay

    def reset_parameters(self):
        # For meta init
        pass

    def _tile(self, x):
        d, ew = x.shape
        if d < self.n_codes:
            n_repeats = (self.n_codes + d - 1) // d
            std = 0.01 / np.sqrt(ew)
            x = x.repeat(n_repeats, 1)
            x = x + torch.randn_like(x) * std
        return x

    def _init_embeddings(self, z):
        # z: [b, t, c]
        self._need_init = False
        flat_inputs = z.view(-1, self.embedding_dim)
        y = self._tile(flat_inputs)

        y.shape[0]
        _k_rand = y[torch.randperm(y.shape[0])][: self.n_codes]
        if dist.is_initialized():
            dist.broadcast(_k_rand, 0)
        self.embeddings.data.copy_(_k_rand)  # pyright: ignore[reportCallIssue]
        self.z_avg.data.copy_(_k_rand)  # pyright: ignore[reportCallIssue]
        self.N.data.copy_(torch.ones(self.n_codes))  # pyright: ignore[reportCallIssue]

    def forward(self, z):
        # z: [b, t, c]
        if self._need_init and self.training and not self.freeze_codebook:
            self._init_embeddings(z)
        # z is of shape [batch_size, sequence length, channels]
        flat_inputs = z.view(-1, self.embedding_dim)
        distances = (
            (flat_inputs**2).sum(dim=1, keepdim=True)
            - 2 * flat_inputs @ self.embeddings.t()  # pyright: ignore[reportCallIssue]
            + (self.embeddings.t() ** 2).sum(dim=0, keepdim=True)  # pyright: ignore[reportCallIssue]
        )  # [bt, c]

        encoding_indices = torch.argmin(distances, dim=1)
        encoding_indices = encoding_indices.view(*z.shape[:2])  # [b, t, ncode]

        embeddings = F.embedding(
            encoding_indices,
            self.embeddings,  # pyright: ignore[reportArgumentType]
        )  # [b, t, c] # pyright: ignore[reportArgumentType]

        commitment_loss = 0.25 * F.mse_loss(z, embeddings.detach())

        # EMA codebook update
        if self.training and not self.freeze_codebook:
            assert False, "Not implemented"
        embeddings_st = (embeddings - z).detach() + z

        return embeddings_st, encoding_indices, commitment_loss

    def dictionary_lookup(self, encodings):
        embeddings = F.embedding(encodings, self.embeddings)  # pyright: ignore[reportArgumentType]
        return embeddings

    def soft_codebook_lookup(self, weights: torch.Tensor) -> torch.Tensor:
        return weights @ self.embeddings  # pyright: ignore[reportOperatorIssue]
