#!/usr/bin/env python3
"""Functions and equations related to BJT amplifiers."""
from dataclasses import dataclass

from eeepy.utils.constants import VT as V_T


def calculate_g_m(I_C: float):
    """Calculates the g_m of the Small-Signal Parameters."""
    return I_C / V_T


def calculate_r_o(V_A: float, I_C: float):
    """Calculates the r_o of the Small-Signal Parameters."""
    return V_A / I_C


def calculate_r_pi(I_B: float):
    """Calculates the r_pi of the Small-Signal Parameters."""
    return V_T / I_B


@dataclass
class SSParameters:
    """Class containing all the Small-Signal Parameters."""
    V_A: float
    I_C: float
    I_B: float

    def g_m(self) -> float:
        """The g_m of the Small-Signal Parameters."""
        return calculate_g_m(self.I_C)

    def r_o(self) -> float:
        """The r_o of the Small-Signal Parameters."""
        return calculate_r_o(self.V_A, self.I_C)

    def r_pi(self) -> float:
        """The r_pi of the Small-Signal Parameters."""
        return calculate_r_pi(self.I_B)
