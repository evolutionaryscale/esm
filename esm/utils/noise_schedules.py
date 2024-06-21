import math

import torch


def cosine_schedule(t: torch.Tensor):
    # t is a tensor of size (batch_size,) with values between 0 and 1. This is the
    # schedule used in the MaskGIT paper
    return torch.cos(t * math.pi * 0.5)


def cubic_schedule(t):
    return 1 - t**3


def linear_schedule(t):
    return 1 - t


def square_root_schedule(t):
    return 1 - torch.sqrt(t)


def square_schedule(t):
    return 1 - t**2


NOISE_SCHEDULE_REGISTRY = {
    "cosine": cosine_schedule,
    "linear": linear_schedule,
    "square_root_schedule": square_root_schedule,
    "cubic": cubic_schedule,
    "square": square_schedule,
}
