from __future__ import annotations

import typing as T
from dataclasses import dataclass

import torch
from typing_extensions import Self

from esm.utils.misc import fp32_autocast_context


@T.runtime_checkable
class Rotation(T.Protocol):
    @classmethod
    def identity(cls, shape: tuple[int, ...], **tensor_kwargs) -> Self: ...

    @classmethod
    def random(cls, shape: tuple[int, ...], **tensor_kwargs) -> Self: ...

    def __getitem__(self, idx: T.Any) -> Self: ...

    @property
    def tensor(self) -> torch.Tensor:
        # We claim that this should be zero-cost abstraction that returns the raw tensor backing this
        # object. The raw tensor should always have exactly 1 more dim than self.shape, which should be
        # implemented using reshaping
        ...

    @property
    def shape(self) -> torch.Size:
        # The "shape" of the rotation, as if it was a torch.tensor object
        # This means that 1x4 quaternions are treated as size (1,) for example
        ...

    def as_matrix(self) -> RotationMatrix: ...

    def compose(self, other: Self) -> Self:
        # To be safe, we force users to explicitly convert between rotation types.
        ...

    def convert_compose(self, other: Self) -> Self:
        # This function will automatically convert between types of rotations
        ...

    def apply(self, p: torch.Tensor) -> torch.Tensor:
        # rotates points by this rotation object
        ...

    def invert(self) -> Self: ...

    @property
    def dtype(self) -> torch.dtype:
        return self.tensor.dtype

    @property
    def device(self) -> torch.device:
        return self.tensor.device

    @property
    def requires_grad(self) -> bool:
        return self.tensor.requires_grad

    @classmethod
    def _from_tensor(cls, t: torch.Tensor) -> Self:
        # This function exists to simplify the below functions, esp type signatures
        # Its implementation is different from Affine3D.from_tensor and does not
        # autodetect rotation types.
        return cls(t)  # type: ignore

    def to(self, **kwargs) -> Self:
        return self._from_tensor(self.tensor.to(**kwargs))

    def detach(self, *args, **kwargs) -> Self:
        return self._from_tensor(self.tensor.detach(**kwargs))

    def tensor_apply(self, func) -> Self:
        # Applys a function to the underlying tensor
        return self._from_tensor(
            torch.stack([func(x) for x in self.tensor.unbind(dim=-1)], dim=-1)
        )


class RotationMatrix(Rotation):
    def __init__(self, rots: torch.Tensor):
        if rots.shape[-1] == 9:
            rots = rots.unflatten(-1, (3, 3))
        assert rots.shape[-1] == 3
        assert rots.shape[-2] == 3
        # Force full precision
        self._rots = rots.to(torch.float32)

    @classmethod
    def identity(cls, shape, **tensor_kwargs):
        rots = torch.eye(3, **tensor_kwargs)
        rots = rots.view(*[1 for _ in range(len(shape))], 3, 3)
        rots = rots.expand(*shape, -1, -1)
        return cls(rots)

    @classmethod
    def random(cls, shape, **tensor_kwargs):
        v1 = torch.randn((*shape, 3), **tensor_kwargs)
        v2 = torch.randn((*shape, 3), **tensor_kwargs)
        return cls(_graham_schmidt(v1, v2))

    def __getitem__(self, idx: T.Any) -> RotationMatrix:
        indices = (idx,) if isinstance(idx, int) or idx is None else tuple(idx)
        return RotationMatrix(self._rots[indices + (slice(None), slice(None))])

    @property
    def shape(self) -> torch.Size:
        return self._rots.shape[:-2]

    def as_matrix(self) -> RotationMatrix:
        return self

    def compose(self, other: RotationMatrix) -> RotationMatrix:
        with fp32_autocast_context(self._rots.device.type):
            return RotationMatrix(self._rots @ other._rots)

    def convert_compose(self, other: Rotation):
        return self.compose(other.as_matrix())

    def apply(self, p: torch.Tensor) -> torch.Tensor:
        with fp32_autocast_context(self.device.type):
            if self._rots.shape[-3] == 1:
                # This is a slight speedup over einsum for batched rotations
                return p @ self._rots.transpose(-1, -2).squeeze(-3)
            else:
                # einsum way faster than bmm!
                return torch.einsum("...ij,...j", self._rots, p)

    def invert(self) -> RotationMatrix:
        return RotationMatrix(self._rots.transpose(-1, -2))

    @property
    def tensor(self) -> torch.Tensor:
        return self._rots.flatten(-2)

    def to_3x3(self) -> torch.Tensor:
        return self._rots

    @staticmethod
    def from_graham_schmidt(
        x_axis: torch.Tensor, xy_plane: torch.Tensor, eps: float = 1e-12
    ) -> RotationMatrix:
        # A low eps here is necessary for good stability!
        return RotationMatrix(_graham_schmidt(x_axis, xy_plane, eps))


@dataclass(frozen=True)
class Affine3D:
    trans: torch.Tensor
    rot: Rotation

    def __post_init__(self):
        assert self.trans.shape[:-1] == self.rot.shape

    @staticmethod
    def identity(
        shape_or_affine: T.Union[tuple[int, ...], "Affine3D"],
        rotation_type: T.Type[Rotation] = RotationMatrix,
        **tensor_kwargs,
    ):
        # Creates a new identity Affine3D object with a specified shape
        # or the same shape as another Affine3D object.
        if isinstance(shape_or_affine, Affine3D):
            kwargs = {"dtype": shape_or_affine.dtype, "device": shape_or_affine.device}
            kwargs.update(tensor_kwargs)
            shape = shape_or_affine.shape
            rotation_type = type(shape_or_affine.rot)
        else:
            kwargs = tensor_kwargs
            shape = shape_or_affine
        return Affine3D(
            torch.zeros((*shape, 3), **kwargs), rotation_type.identity(shape, **kwargs)
        )

    @staticmethod
    def random(
        shape: tuple[int, ...],
        std: float = 1,
        rotation_type: T.Type[Rotation] = RotationMatrix,
        **tensor_kwargs,
    ) -> "Affine3D":
        return Affine3D(
            trans=torch.randn((*shape, 3), **tensor_kwargs).mul(std),
            rot=rotation_type.random(shape, **tensor_kwargs),
        )

    def __getitem__(self, idx: T.Any) -> "Affine3D":
        indices = (idx,) if isinstance(idx, int) or idx is None else tuple(idx)
        return Affine3D(trans=self.trans[indices + (slice(None),)], rot=self.rot[idx])

    @property
    def shape(self) -> torch.Size:
        return self.trans.shape[:-1]

    @property
    def dtype(self) -> torch.dtype:
        return self.trans.dtype

    @property
    def device(self) -> torch.device:
        return self.trans.device

    @property
    def requires_grad(self) -> bool:
        return self.trans.requires_grad

    def to(self, **kwargs) -> "Affine3D":
        return Affine3D(self.trans.to(**kwargs), self.rot.to(**kwargs))

    def detach(self, *args, **kwargs) -> "Affine3D":
        return Affine3D(self.trans.detach(**kwargs), self.rot.detach(**kwargs))

    def tensor_apply(self, func) -> "Affine3D":
        # Applys a function to the underlying tensor
        return self.from_tensor(
            torch.stack([func(x) for x in self.tensor.unbind(dim=-1)], dim=-1)
        )

    def as_matrix(self):
        return Affine3D(trans=self.trans, rot=self.rot.as_matrix())

    def compose(self, other: "Affine3D", autoconvert: bool = False):
        rot = self.rot
        new_rot = (rot.convert_compose if autoconvert else rot.compose)(other.rot)
        new_trans = rot.apply(other.trans) + self.trans
        return Affine3D(trans=new_trans, rot=new_rot)

    def compose_rotation(self, other: Rotation, autoconvert: bool = False):
        return Affine3D(
            trans=self.trans,
            rot=(self.rot.convert_compose if autoconvert else self.rot.compose)(other),
        )

    def scale(self, v: torch.Tensor | float):
        return Affine3D(self.trans * v, self.rot)

    def mask(self, mask: torch.Tensor, with_zero=False):
        # Returns a transform where True positions in mask is identity
        if with_zero:
            tensor = self.tensor
            return Affine3D.from_tensor(
                torch.zeros_like(tensor).where(mask[..., None], tensor)
            )
        else:
            identity = self.identity(
                self.shape,
                rotation_type=type(self.rot),
                device=self.device,
                dtype=self.dtype,
            ).tensor
            return Affine3D.from_tensor(identity.where(mask[..., None], self.tensor))

    def apply(self, p: torch.Tensor) -> torch.Tensor:
        return self.rot.apply(p) + self.trans

    def invert(self):
        inv_rot = self.rot.invert()
        return Affine3D(trans=-inv_rot.apply(self.trans), rot=inv_rot)

    @property
    def tensor(self) -> torch.Tensor:
        return torch.cat([self.rot.tensor, self.trans], dim=-1)

    @staticmethod
    def from_tensor(t: torch.Tensor) -> "Affine3D":
        match t.shape[-1]:
            case 4:
                # Assume tensor 4x4 for backward compat with alphafold
                trans = t[..., :3, 3]
                rot = RotationMatrix(t[..., :3, :3])
            case 12:
                trans = t[..., -3:]
                rot = RotationMatrix(t[..., :-3].unflatten(-1, (3, 3)))
            case _:
                raise RuntimeError(
                    f"Cannot detect rotation fromat from {t.shape[-1] -3}-d flat vector"
                )
        return Affine3D(trans, rot)

    @staticmethod
    def from_tensor_pair(t: torch.Tensor, r: torch.Tensor) -> "Affine3D":
        return Affine3D(t, RotationMatrix(r))

    @staticmethod
    def from_graham_schmidt(
        neg_x_axis: torch.Tensor,
        origin: torch.Tensor,
        xy_plane: torch.Tensor,
        eps: float = 1e-10,
    ):
        # The arguments of this function is for parity with AlphaFold
        x_axis = origin - neg_x_axis
        xy_plane = xy_plane - origin
        return Affine3D(
            trans=origin, rot=RotationMatrix.from_graham_schmidt(x_axis, xy_plane, eps)
        )

    @staticmethod
    def cat(affines: list["Affine3D"], dim: int = 0):
        if dim < 0:
            dim = len(affines[0].shape) + dim
        return Affine3D.from_tensor(torch.cat([x.tensor for x in affines], dim=dim))


def _graham_schmidt(x_axis: torch.Tensor, xy_plane: torch.Tensor, eps: float = 1e-12):
    # A low eps here is necessary for good stability!
    with fp32_autocast_context(x_axis.device.type):
        e1 = xy_plane

        denom = torch.sqrt((x_axis**2).sum(dim=-1, keepdim=True) + eps)
        x_axis = x_axis / denom
        dot = (x_axis * e1).sum(dim=-1, keepdim=True)
        e1 = e1 - x_axis * dot
        denom = torch.sqrt((e1**2).sum(dim=-1, keepdim=True) + eps)
        e1 = e1 / denom
        e2 = torch.cross(x_axis, e1, dim=-1)

        rots = torch.stack([x_axis, e1, e2], dim=-1)

        return rots


def build_affine3d_from_coordinates(
    coords: torch.Tensor,  # (N, CA, C).
) -> tuple[Affine3D, torch.Tensor]:
    _MAX_SUPPORTED_DISTANCE = 1e6
    coord_mask = torch.all(
        torch.all(torch.isfinite(coords) & (coords < _MAX_SUPPORTED_DISTANCE), dim=-1),
        dim=-1,
    )

    def atom3_to_backbone_affine(bb_positions: torch.Tensor) -> Affine3D:
        N, CA, C = bb_positions.unbind(dim=-2)
        return Affine3D.from_graham_schmidt(C, CA, N)

    coords = coords.clone().float()
    coords[~coord_mask] = 0

    # NOTE(thayes): If you have already normalized the coordinates, then
    # the black hole affine translations will be zeros and the rotations will be
    # the identity.
    average_per_n_ca_c = coords.masked_fill(~coord_mask[..., None, None], 0).sum(1) / (
        coord_mask.sum(-1)[..., None, None] + 1e-8
    )
    affine_from_average = atom3_to_backbone_affine(
        average_per_n_ca_c.float()
    ).as_matrix()

    B, S, _, _ = coords.shape
    assert isinstance(B, int)
    assert isinstance(S, int)
    affine_rot_mats = affine_from_average.rot.tensor[..., None, :].expand(B, S, 9)
    affine_trans = affine_from_average.trans[..., None, :].expand(B, S, 3)

    # We use the identity rotation whereever we have no coordinates. This is
    # important because otherwise the rotation matrices will be all zeros, which
    # will cause collapse in the distance/direction attention mechanism.
    identity_rot = RotationMatrix.identity(
        (B, S), dtype=torch.float32, device=coords.device, requires_grad=False
    )
    affine_rot_mats = affine_rot_mats.where(
        coord_mask.any(-1)[..., None, None], identity_rot.tensor
    )
    black_hole_affine = Affine3D(affine_trans, RotationMatrix(affine_rot_mats))

    affine = atom3_to_backbone_affine(coords.float())
    affine = Affine3D.from_tensor(
        affine.tensor.where(coord_mask[..., None], black_hole_affine.tensor)
    )

    return affine, coord_mask
