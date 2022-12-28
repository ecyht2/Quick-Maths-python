#!/usr/bin/env python3
"""Functions used to calculate forward converter related values."""


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
class DutyCycleIsolated(DutyCycle):
    """Duty cycle equation for an isolated forward converter.

    Formula: d = \frac {N_{1}}{N_{2}} \frac {V_{o}}{V_{s}}
    """
    def __new__(cls, V_o: float, V_s: float,
                N1: float, N2: float):
        return super().__new__(cls, (V_o * N1) / (V_s * N2))

    def __init__(self, V_o: float, V_s: float,
                 N1: float, N2: float):
        """Duty cycle equation for an isolated forward converter.

        :param V_o: Average voltage output of the converter.
        :param V_s: Supply voltage of the converter.
        :param N1: The number of primary windings in the transformer.
        :param N2: The number of secondary windings in the transformer.
        """
        super().__init__()

    @classmethod
    def from_turns_ratio(cls, V_o: float, V_s: float,
                         turns_ratio: float):
        """Calculates the duty cycle using the turns ratio instead of number of
        windings.

        :param V_o: Average voltage output of the converter.
        :param V_s: Supply voltage of the converter.
        :param turns_ratio: The turns ratio of the transformer.
        """
        return cls(V_o, V_s, turns_ratio, 1)

    @staticmethod
    def V_o(d: DutyCycle, V_s: float, N1: float, N2: float) -> float:
        """Calculates the V_o (average voltage) of the forward converter.

        :param d: Duty Cycle of the forward converter.
        :param V_s: Supply voltage of the forward converter.
        :param N1: The number of primary windings in the transformer.
        :param N2: The number of secondary windings in the transformer.
        :return: Average voltage output of the converter.
        """
        return N2 / N1 * d * V_s

    @staticmethod
    def V_s(d: DutyCycle, V_o: float, N1: float, N2: float) -> float:
        """Calculates the V_s (supply voltage) of the forward converter.

        :param d: Duty Cycle of the forward converter.
        :param V_o: Average voltage of the forward converter.
        :param N1: The number of primary windings in the transformer.
        :param N2: The number of secondary windings in the transformer.
        :return: The supply voltage of the converter.
        """
        return N1 / N2 * V_o / d

    @staticmethod
    def turns_ratio(d: DutyCycle, V_o: float, V_s: float) -> float:
        """Calculates the turns ratio of the forward converter.

        :param d: Duty Cycle of the forward converter.
        :param V_s: Supply voltage of the forward converter.
        :param V_o: Average voltage output of the converter.
        """
        return d * V_s / V_o


class IsolatedForwardConverter:
    def __init__(self):
        ...
