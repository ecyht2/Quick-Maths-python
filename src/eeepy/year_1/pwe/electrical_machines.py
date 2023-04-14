#!/usr/bin/env python3
"""Module containing functions and equations relating to electrical
 machines."""
from math import pi

from eeepy.utils.units import rad_per_sec


# Electrical Machine
def induced_EMF(P: float, phi: float, Z: float, a: float,
                omega: float = 0, N: float = 0) -> float:
    """
    Calculate the induced EMF for a DC machine
    """
    # pylint: disable = too-many-arguments
    if omega > 0:
        Ea = P * Z * phi * omega / (2 * pi * a)
    elif N > 0:
        Ea = P * Z * phi * N / (60 * a)
    else:
        raise ValueError("No omega or N given")
    return Ea


def machine_constant(P: float, Z: float, a: float) -> float:
    """
    Calculates the machine constant of DC machines
    """
    return Z * P / (2 * pi * a)


def induced_EMF_Km(Km: float, phi: float, omega: float) -> float:
    """
    Calculate the induced EMF for a DC machine using machine constant
    """
    return Km * phi * omega


def armature_power(Ta: float, rpm: float = 0) -> float:
    """
    Calculates the mechanical power developed by a DC motor armature
    """
    omega = rad_per_sec(rpm)
    return Ta * omega


def dc_motor_dev_torque_Km(Km: float, Ia: float, phi: float) -> float:
    """
    calculates torque produce (developed) by a dc motor using machine constant
    """
    return Km * Ia * phi


def dc_motor_dev_torque_pdev(Pdev: float,
                             rpm: float = 0) -> float:
    """
    calculates torque produce (developed) by a dc motor using machine constant
    """
    omega = rad_per_sec(rpm)
    return Pdev / omega


def dc_motor_torque(Km: float, Ia: float, phi: float) -> float:
    """
    calculates torque produce by a dc motor
    """
    return Km * Ia * phi

# def developed_power_shunt(Vdc: float, Im: float,
#                           Rf: float, Ra: float,
#                           Vbrush: float = 0) -> float:
#     """
#     Calculates the devolped power of a shunt DC motor
#     """
#     If = Vdc / RF
#     Ia = Im - If
#     Eb = Vdc - Ia * Ra - Vbrush
#     Pdev = electrical_power(Eb, Ia)
#     return Pdev
