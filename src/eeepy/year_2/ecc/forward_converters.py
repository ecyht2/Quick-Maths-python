#!/usr/bin/env python3
"""Functions used to calculate forward converter related values."""


# Non-Isolated
class DutyCycle(float):
    """Duty cycle equation for a forward converter.
    d = V_o / V_s
    """
    def __new__(cls, V_o: float, V_s: float):
        return super().__new__(cls, V_o / V_s)

    def __init__(self, V_o: float, V_s: float):
        """Duty cycle equation for a forward converter.

        :param V_o: Average voltage output of the converter
        :param V_s: Supply voltage of the converter
        """
        super().__init__()

    @staticmethod
    def V_o(d: float, V_s: float) -> float:
        """Calculates the V_o (average voltage) of a forward converter.

        :param d: Duty Cycle of the forward converter
        :param V_s: Supply voltage of the forward converter
        :return:
        """
        return d * V_s

    @staticmethod
    def V_s(d: float, V_o: float) -> float:
        """Calculates the V_s (supply voltage) of a forward converter.

        :param d: Duty Cycle of the forward converter
        :param V_o: Average voltage of the forward converter
        :return:
        """
        return V_o / d


class VTAInductor(float):
    def __new__(cls, d: float, T: float, V_o):
        return super().__new__(cls, (1 - d) * T * V_o)

    def __init__(cls, d: float, T: float, V_o):
        super().__init__()


class VoltageRipple(float):
    """Duty cycle equation for a forward converter.
    ripple = dV / V_o
    """
    def __new__(cls, dV: float):
        return super().__new__(cls, dV)

    def __init__(self, dV: float):
        """Duty cycle equation for a forward converter.

        :param V_o: Average voltage output of the converter
        :param V_s: Supply voltage of the converter
        """
        super().__init__()

    @staticmethod
    def capacitance(dI: float, f: float, dV: float):
        return dI / (8 * f * dV)


class LCRIT(float):
    """Finds the L_crit of a forward converter."""
    def __new__(cls, L_crit: float):
        return super().__new__(cls, L_crit)

    def __init__(self, L_crit: float):
        super().__init__()

    @classmethod
    def from_VTA(cls, VTA: float, dI: float):
        """"""
        L_crit = VTA / dI
        return cls(L_crit)

    @classmethod
    def from_R(cls, d: float, R: float, f_sw: float):
        """"""
        L_crit = (1 - d) * R / (2 * f_sw)
        return cls(L_crit)

    @classmethod
    def from_duty_cycle(cls, d: float, V_o: float, T: float, dI: float):
        """"""
        VTA = VTAInductor(d, T, V_o)
        L_crit = VTA / dI
        return cls(L_crit)


class CurrentRipple(float):
    def __new__(cls, d: float, T: float, V_o: float, L: float):
        VTA = VTAInductor(d, T, V_o)
        dI = VTA / L
        super().__new__(cls, dI)

    def __init__(self, d: float, V: float, T):
        super().__init__()

    @staticmethod
    def inductance(dI: float, d: float, T: float, V_o: float):
        VTA = VTAInductor(d, T, V_o)
        L = VTA / dI
        return L


# Isolated
def v_o_isolated(d: float, V_s: float, N2: float, N1: float):
    """Caclculates the average output voltage of an isolated forward converter.
    """


class RMSCurrentInductor(float):
    # sqrt((I_1**2 + I_1 * I_2 + I_2**2) * d / 3)
    def __new__(cls, I_1, I_2, d):
        ...

    def __init__(self, I_1, I_2, d):
        ...
