#!/usr/bin/env python3
"""Functions used to calculate forward converter related values."""
import math


class DutyCycle(float):
    ...


# Non-Isolated
class DutyCyleNonIsolated(DutyCycle):
    """Duty cycle equation for a forward converter.

    Formula: d = \frac {V_{o}}{V_{s}}
    """
    def __new__(cls, V_o: float, V_s: float):
        return super().__new__(cls, V_o / V_s)

    def __init__(self, V_o: float, V_s: float):
        """Duty cycle equation for a forward converter.

        :param V_o: Average voltage output of the converter.
        :param V_s: Supply voltage of the converter.
        """
        super().__init__()

    @staticmethod
    def V_o(d: DutyCycle, V_s: float) -> float:
        """Calculates the V_o (average voltage) of a forward converter.

        :param d: Duty Cycle of the forward converter.
        :param V_s: Supply voltage of the forward converter.
        :return: The average voltage output of the converter.
        """
        return d * V_s

    @staticmethod
    def V_s(d: DutyCycle, V_o: float) -> float:
        """Calculates the V_s (supply voltage) of a forward converter.

        :param d: Duty Cycle of the forward converter.
        :param V_o: Average voltage of the forward converter.
        :return: The supply voltage of the converter.
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


class ForwardConverter:
    def __init__(self, Vs: float, Vo: float, f: float, d: DutyCycle,
                 I_o_max: float, I_o_min: float):
        self.Vs: float = Vs
        self.Vo: float = Vo
        self.f: float = f
        self.d: DutyCycle = d
        self.I_o_max: float = I_o_max
        self.I_o_min: float = I_o_min

    @property
    def T(self):
        """The period of the forward converter."""
        self._T = 1 / self.f


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


class IsolatedForwardConverter:
    def __init__(self):
        ...


def i_RMS_Q(i1: float, i2: float, d: DutyCycle):
    """Calculates the peak to peak RMS value of the transistor.

    Formula: \\sqrt {\frac{i_{1}^{2} + i_{1} * i_{2} + i_{2}^{2}}{3}

    Parameters
    ----------
    i1: float
        The first current.
    i2: float
        The second current.
    d: DutyCycle
        The duty cycle.
    """
    ans = (i1**2 + i1 * i2 + i2**2) * d
    ans /= 3
    return math.sqrt(ans)
