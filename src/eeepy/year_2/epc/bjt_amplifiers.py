#!/usr/bin/env python3
from eeepy.utils.Constants import VT as V_T
from dataclasses import dataclass


def calculate_g_m(I_C: float):
    return I_C / V_T


def calculate_r_o(V_A: float, I_C: float):
    return V_A / I_C


def calculate_r_pi(I_B: float):
    return V_T / I_B


@dataclass
class SSParameters:
    """Class containing all the Small-Signal Parameters."""
    V_A: float
    I_C: float
    I_B: float

    def g_m(self) -> float:
        return calculate_g_m(self.I_C)

    def r_o(self) -> float:
        return calculate_r_o(self.V_A, self.I_C)

    def r_pi(self) -> float:
        return calculate_r_pi(self.I_B)
