#!/usr/bin/env python3
"""Functions and Equations related to AC synchronous machines."""
import math
from dataclasses import dataclass


class SynchronizedFrequency:
    """The frequency produced by synchrnonous generator.

    f = P N_s / 120
    --------------------------------------------
    f is the frequency produced.
    P is the number of poles.
    N_s is the synchronous speed of the machine.

    f = P ɷ_s / (4π)
    ----------------------------
    f is the frequency produced.
    P is the number of poles.
    ɷ_s is the synchronous speed of the machine.
    """
    @staticmethod
    def frequency(P: int, N_s: float) -> float:
        """Calculates the frequency.

        :param P: The number of poles.
        :param N_s: The synchronous speed of the machine.
        :returns: The frequency produced.
        """
        return P * N_s / 120

    @staticmethod
    def frequency_omega(P: int, w_s: float) -> float:
        """Calculates the frequency.

        :param P: The number of poles.
        :param w_s: The synchronous speed of the machine.
        :returns: The frequency produced.
        """
        return P * w_s / 4 / math.pi

    @staticmethod
    def number_of_poles(f: float, N_s: float) -> int:
        """Calculates the number of poles.

        :param f: The frequency produced.
        :param N_s: The synchronous speed of the machine.
        :returns: The number of poles.
        """
        return 120 * f / N_s

    @staticmethod
    def number_of_poles_omega(f: float, w_s: float) -> int:
        """Calculates the number of poles.

        :param f: The frequency produced.
        :param w_s: The synchronous speed of the machine.
        :returns: The number of poles.
        """
        return 4 * math.pi * f / w_s

    @staticmethod
    def synchronous_speed(f: float, P: float) -> int:
        """Calculates the synchronous speed.

        :param f: The frequency produced.
        :param P: The number of poles.
        :returns: The number of poles.
        """
        return 120 * f / P


@dataclass
class EMFInduced:
    """Calculates the EMF induced."""


class VoltageRegulation:
    """Equation to calculate voltage regulation.

    regulation = E - V / V
    ----------------------
    regulation is the voltage regulation.
    E is the EMF induced.
    V is the terminal voltage.
    """
    @staticmethod
    def regulation(E: float, V: float) -> float:
        """Calculates the voltage regulation.

        :param E: The EMF induced.
        :param V: The terminal voltage.
        :returns: The voltage regulation.
        """
        return (E - V) / V

    @staticmethod
    def regulation_percent(E: float, V: float) -> float:
        """Calculates the voltage regulation percentage.

        :param E: The EMF induced.
        :param V: The terminal voltage.
        :returns: The voltage regulation in percentage.
        """
        return 100 * (E - V) / V
