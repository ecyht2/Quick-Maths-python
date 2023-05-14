#!/usr/bin/env python3
"""Functions and equations relating to magnetic circuits."""
import math

from eeepy.utils.constants import mu0


class AmperesLaw:
    """The Ampere's Law.

    Formula: :math:`MMF = N I`

    MMF is the magneto motive force.

    N is the number of turns.

    I is the current through the turns.
    """
    @staticmethod
    def mmf(N: int, I: float) -> float:
        """Calculates the MMF (Magneto Motive Force).

        :param N: The number of turns.
        :param I: The current through the turn.
        :returns: The magneto motive force.
        """
        return N * I

    @staticmethod
    def number_of_coils(mmf: float, I: float) -> int:
        """Calculates the number of coils.

        :param mmf: The magneto motive force.
        :param I: The current.
        :returns: The number of coils (rounded up).
        """
        return math.ceil(mmf / I)

    @staticmethod
    def current(mmf: float, N: int) -> float:
        """Calculates the current.

        :param mmf: The magneto motive force.
        :param N: The number of coils.
        :returns: The current through the coils.
        """
        return mmf / N


class HopkinsonLaw:
    """The Hopkinson's law equation.

    Formula: :math:`MMF = \\Phi R`

    MMF is the magneto motive force.

    R is the reluctance.

    :math:`\\Phi` is the total flux.
    """
    @staticmethod
    def reluctance(mmf: float, total_flux: float) -> float:
        """Calculates the reluctance from mmf (MMF) and
        total_flux (:math:`\\Phi`)."""
        return mmf / total_flux

    @staticmethod
    def total_flux(mmf: float, reluctance: float) -> float:
        """Calculates the total flux (:math:`\\Phi`) from mmf (MMF) and
        reluctance (R)."""
        return mmf / reluctance

    @staticmethod
    def mmf(total_flux: float, reluctance: float) -> float:
        """Calculates the MMF from total_flux (:math:`\\Phi`) and
        reluctance (R)."""
        return total_flux * reluctance


class ReluctanceEquation:
    """The reluctance of a magnetic circuit.

    Formula: :math:`R = l / (\\mu A)`

    R is the reluctance.

    l is the length of the circuit.

    :math:`\\mu` is the permeability of the material
        (:math:`\\mu = {\\mu}_0 {\\mu}_r`).

    A is the cross-sectional area of the circuit.
    """
    @staticmethod
    def reluctance(l: float, mu: float, A: float) -> float:
        """Calculates the reluctance from permeability.

        :param l: The length of the circuit.
        :param mu: The permeability of the material.
        :param A: The cross-sectional area of the circuit.
        :returns: The reluctance.
        """
        return l / (mu * A)

    @staticmethod
    def reluctance_relative_permeability(l: float,
                                         mu_r: float, A: float) -> float:
        """Calculates the reluctance from permeability.

        :param l: The length of the circuit.
        :param mu_r: The relative permeability of the material.
        :param A: The cross-sectional area of the circuit.
        :returns: The reluctance.
        """
        return l / (mu_r * mu0 * A)

    @staticmethod
    def length(reluctance: float, mu: float, A: float) -> float:
        """Calculates the length of the circuit.

        :param reluctance: The reluctance.
        :param mu: The permeability of the material.
        :param A: The cross-sectional area of the circuit.
        :returns: The length of the circuit.
        """
        return reluctance * mu * A
