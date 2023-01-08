#!/usr/bin/env python3
"""Functions related to Power Electronics."""


class TimeAreaRuleL(float):
    """Time Area Rule equation for an inductor. (dI = VTA/L)"""
    def __new__(cls, area: float, L: float):
        """Find the ΔI of a given VTA and L.

        Parameters
        ----------
        area: float
            The time Area VTA. (voltage time area)
        L: float
            The inductance.
        """
        return super().__new__(cls, area / L)

    def __init__(self, area: float, value: float):
        """Find the ΔI of a given VTA and L.

        Parameters
        ----------
        area: float
            The time Area VTA. (voltage time area)
        L: float
            The inductance.
        """
        super().__init__()

    @staticmethod
    def find_inductance(VTA: float, dI: float):
        """Find the inductance needed for a given VTA and ΔI. (L = VTA / ΔI)

        Prameters
        ---------
        VTA: float
            The VTA. (voltage time area)
        dI: float
            The change in current ΔI.
        """
        return VTA / dI

    @staticmethod
    def find_VTA(dI: float, L: float):
        """Find the VTA for a given ΔI and L. (VTA = ΔI * L)

        Parameters
        ----------
        dI: float
            The change in current ΔI.
        L: float
            The inductance.
        """
        return dI * L


class TimeAreaRuleC(float):
    """Time Area Rule equation for an capacitor. (dV = ITA/C)"""
    def __new__(cls, area: float, C: float):
        """Find the ΔV of a given ITA and C.

        Parameters
        ----------
        area: float
            The ITA. (current time area)
        C: float
            The capacitance.
        """
        return super().__new__(cls, area / C)

    def __init__(self, area: float, value: float):
        """Find the ΔV of a given ITA and C.

        Parameters
        ----------
        area: float
            The ITA. (current time area)
        C: float
            The capacitance.
        """
        super().__init__()

    @staticmethod
    def find_capacitance(ITA: float, dV: float):
        """Find the capacitance needed for a given ITA and ΔV. (C = ITA / ΔV)

        Prameters
        ---------
        ITA: float
            The ITA. (current time area)
        dV: float
            The change in voltage ΔV.
        """
        return ITA / dV

    @staticmethod
    def find_ITA(dV: float, C: float):
        """Find the ITA for a given ΔV and C. (ITA = ΔV * C)

        Parameters
        ----------
        dV: float
            The change in voltage ΔV.
        C: float
            The capacitance.
        """
        return dV * C


def average_voltage(duty_cycle: float, voltage: float) -> float:
    """Finds the average voltage."""
    return duty_cycle * voltage


def average_current(duty_cycle: float, voltage: float,
                    resistance: float) -> float:
    """Finds the average current."""
    return duty_cycle * voltage / resistance


def peak_and_min(value: float, change: float) -> (float, float):
    """Calculates The peak and minimum value given the change and average.

    :param value: The average value.
    :param change: The total change. (difference between peak and min)
    :return: A tuple of the minimum and peak value in the form (min, peak).
    """
    peak = value + change / 2
    minimum = value - change / 2
    return (minimum, peak)
