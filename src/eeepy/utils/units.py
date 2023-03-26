#!/usr/bin/env python3
"""Useful functions to convert one unit to another."""
from math import log10, pi


# Power and Energy
# RMS
def RMS_sinusoidal(value: float) -> float:
    """Converts a sinusoidal sigal peak value into RMS."""
    return value / 2**0.5


def RMS_sinusoidal_reverse(value: float) -> float:
    """Converts an RMS value into sinusoidal sigal peak value."""
    return value * 2**0.5


# Angular Frequency
def angular_frequency(f: float = 0, T: float = 0) -> float:
    """Calculates the angular frequency (2πf) of a sinusoidal signal."""
    omega = 0
    if f > 0:
        omega = 2 * pi * f
    elif T > 0:
        omega = 2 * pi / T
    else:
        raise ValueError("No f or T given")
    return omega


def angular_frequency_reverse(omega: float, frequency: float = True) -> float:
    """Calculates the frequency or period of a sinusoidal signal given the
    angular frequency.
    """
    if frequency:
        ret = omega / (2 * pi)
    else:
        ret = 2 * pi / omega
    return ret


# Horse Power
def hp_to_W(hp):
    """Calculates the power in watts from horse power."""
    W = 746 * hp
    return W


def W_to_hp(W):
    """Calculates the power in horse power from watts."""
    hp = W / 746
    return hp


# RPM
def rpm(omega: float = 0) -> float:
    """Calculates the angular velocity in rpm from radians per second(ɷ)."""
    return omega * 60 / (2 * pi)


def rad_per_sec(rpm: float = 0) -> float:
    """Calculates the angular velocity in radians per second(ɷ) from rpm."""
    return rpm * 2 * pi / 60


# Engineering Maths
# db Functions
# Power
def db_power(initial, output):
    """Calculates the power in decibels."""
    return 10 * log10(output / initial)


def db_power_reverse(db, initial):
    """Calculates the output power from decibels."""
    db = db / 10
    return 10**db * initial


def db_power_find_initial(db, output):
    """Calculates the input power from decibels."""
    db = db / 10
    return output / 10**db


# Volts/Current
def db_volts(initial, output):
    """Calculates the voltage in decibels."""
    return 20 * log10(output / initial)


def db_volts_reverse(db, initial):
    """Calculates the output voltage from decibels."""
    db = db / 20
    return 10**db * initial


def db_volts_find_initial(db, output):
    """Calculates the input voltage from decibels."""
    db = db / 20
    return output / 10**db


# Info and System
def C_to_kelvin(T: float) -> float:
    """Converts tempreture from Celcius(°C) to Kelvin (K)."""
    return T + 273


def kelvin_to_C(T: float) -> float:
    """Converts tempreture from Kelvin (K) to Celcius(°C)."""
    return T - 273
