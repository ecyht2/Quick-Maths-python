#!/usr/bin/env python3
"""Functions related to Power Electronics."""


class TimeAreaRule(float):
    """Time Area Rule Equations (dI = VTA/L, dV = ITA/C)."""
    def __new__(cls, area: float, value: float):
        """Find the change in I/V

        Parameters
        ----------
        area: float
            Time Area could be VTA (voltage time area) or
            ITA (current time area)
        value: float
            The capcitance or inductance
        """
        return super().__new__(cls, area / value)

    def __init__(self, area: float, value: float):
        """Find the change in I/V

        Parameters
        ----------
        area: float
            Time Area could be VTA (voltage time area) or
            ITA (current time area)
        value: float
            The capcitance or inductance
        """
        super().__init__()

    @classmethod
    def from_inductor(cls, VTA, L):
        return cls(VTA, L)

    @classmethod
    def from_capacitor(cls, ITA, C):
        return cls(ITA, C)


def average_voltage(duty_cycle: float, voltage: float) -> float:
    """Finds the average voltage."""
    return duty_cycle * voltage


def average_current(duty_cycle: float, voltage: float,
                    resistance: float) -> float:
    """Finds the average current."""
    return duty_cycle * voltage / resistance
