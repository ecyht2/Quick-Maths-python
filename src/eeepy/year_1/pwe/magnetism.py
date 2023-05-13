#!/usr/bin/env python3
"""Module containing equation related to magnetism and magnetic circuits."""
from math import pi

from eeepy.utils.constants import mu0
from eeepy.year_1.pwe import imp_ind


class FluxIntensity:
    """Equation for magnetic flux intensity.

    H = N I / l
    -----------
    H is the magnetic flux intensity.
    N is the number of coils.
    I is the current through the coils.
    l is the average path of the magnetic flux.
    """
    @staticmethod
    def H_field_circular(N: int, I: float, r: float) -> float:
        """Calculates the magnetic flux intensity of with circular core.

        :param N: The number of coils.
        :param I: The current through the coils.
        :param r: The radius of the circular core.
        :returns: The magnetic flux intensity.
        """
        return N * I / (2 * pi * r)

    @staticmethod
    def H_field(N: int, i: float, l: float) -> float:
        """Calculates the magnetic flux intensity.

        :param N: The number of coils.
        :param I: The current through the coils.
        :param l: The average path of the magnetic flux.
        :returns: The magnetic flux intensity.
        """
        return N * i / l

    @staticmethod
    def current(H: float, N: int, l: float) -> float:
        """Calculates the current needed to produce a magnetic flux intensity.

        :param H: The magnetic flux intensity.
        :param N: The number of coils.
        :param l: The average path of the magnetic flux.
        :returns: The current through the coils.
        """
        return H * l / N


class FluxDensity:
    """Equation for manetic flux density.

    B = μ_0 μ_r H
    -------------
    B is the magnetifc flux density.
    μ_0 is the permeability of free space.
    μ_r is the relative permeability of the medium.
    H is the magnetic flux intensity.
    """
    @staticmethod
    def B_field_from_current(N: int, I: float, l: float,
                             mur: float = 1) -> float:
        """Calculates the magnetic field density of a coil given the current
        and length of coil.

        :param N: The number of coils.
        :param I: The current through the coils.
        :param l: The average path of the magnetic flux.
        :param mur: The relative permeability of the medium.
        :returns: The magnetifc flux density.
        """
        return mu0 * mur * N * I / l

    @staticmethod
    def B_field(H: float, mur: float = 1) -> float:
        """Calculates the magnetic field density field of a coil given the
        magnetic flux intensity field.

        :param mur: The relative permeability of the medium.
        :param H: The magnetic flux intensity.
        :returns: The magnetifc flux density.
        """
        return mu0 * mur * H

    @staticmethod
    def flux_intensity(B: float, mur: float) -> float:
        """Calculates the magnetic flux intensity.

        :param mur: The relative permeability of the medium.
        :param B: The magnetic flux density.
        :returns: The magnetic flux intensity.
        """
        return B / mur / mu0


def magnetic_force_wire(i: float, B: float, l: float) -> float:
    """
    Calculates the magnetic force on a current carrying wire of length l
    """
    return abs(i * B * l)


def flux(B: float, A: float) -> float:
    """
    Calculates the flux (φ) in a single coil
    """
    return B * A


def flux_linkage(N: int, phi: float) -> float:
    """
    Calculates the flux linkage (λ) of a coil
    """
    return N * phi


def inductance(i: float, linkage: float = 0,
               N: int = 0, phi: float = 0) -> float:
    """
    Calculate the inductance of a coil
    """
    if linkage > 0:
        L = linkage / i
    elif N > 0 and phi > 0:
        L = N * phi / i
    return L


def flux_linkage_two_coil(L: float, i1: float,
                          i2: float, M: float,
                          strengthen: bool = True) -> float:
    """
    Calculates the flux linkage of a two coil system
    """
    if strengthen:
        linkage = L * i1 + M * i2
    else:
        linkage = L * i1 - M * i2
    return linkage


def voltage_two_coil(R: float, I1: float, L: float, dI1: float,
                     I2: float, M: float, strengthen: bool = True) -> float:
    """
    Calculates the flux linkage of a two coil system
    """
    # pylint: disable = too-many-arguments
    V = R * I1
    if strengthen:
        V += L * dI1 + M * I2
    else:
        V += L * I1 - M * I2
    return V


def voltage_two_coil_phasor(R: float, L: float, I: tuple, M: float,
                            omega: float, f: float,
                            strengthen: bool = True) -> float:
    """
    Calculates the flux linkage of a two coil system
    """
    # pylint: disable = too-many-arguments
    if omega > 0:
        pass
    elif f > 0:
        omega = 2 * pi * f
    else:
        raise ValueError("No omega or f or T given")

    Z = R + imp_ind(L, omega=omega)
    V = I[0] * Z
    if strengthen:
        V += M * 1j * omega * I[1]
    else:
        V -= M * 1j * omega * I[1]
    return V
