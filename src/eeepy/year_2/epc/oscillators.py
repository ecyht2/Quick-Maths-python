#!/usr/bin/env python3
"""Functions and equations related to oscilators."""
from eeepy.utils.Constants import pi
from eeepy.utils.helper import sqrt


def wein_bridge_resonance_freq(R: float, C: float) -> float:
    """Calculates the resonance frequency of a Wein Bridge Oscillator."""


def phase_shift_freq(R: float, C: float, N: int) -> float:
    """Calculates the resonance frequency of a Phase Shift Oscillator."""
    return 1 / (2 * pi * sqrt(2 * N) * R * C)


def colpitts_freq(L: float, C_T: float) -> float:
    """Calculates the oscillation frequency of a Colpitts Oscillator."""
    return 1 / (2 * pi * sqrt(L * C_T))


def colpitts_ct(C_1: float, C_2: float) -> float:
    """Calculates the Ct of a Colpitts Oscillator."""
    return (C_1 * C_2) / (C_1 + C_2)


def colpitts_b(C_1: float, C_2: float) -> float:
    """Calculates the B needed for a Colpitts Oscillator."""
    return C_2 / C_1


def colpitts_gain(C_1: float, C_2: float) -> float:
    """Calculates the voltage gain needed for a Colpitts Oscillator."""
    return C_1 / C_2
