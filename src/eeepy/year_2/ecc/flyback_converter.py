#!/usr/bin/env python3
"""Functions and equations used to calculate flyback converter related values.
"""
from .forward_converters import DutyCycle


# Non-Isolated
class DutyCycleNonIsolated(DutyCycle):
    """Duty cycle equation for a non-isolated flyback converter.

    Formula: d = \frac {V_{o}}{V_{s} + V_{o}}
    """
    def __new__(cls, V_o: float, V_s: float):
        return super().__new__(cls, V_o / (V_s + V_o))

    def __init__(self, V_o: float, V_s: float):
        """Duty cycle equation for a non-isolated  flyback converter.

        :param V_o: Average voltage output of the converter.
        :param V_s: Supply voltage of the converter.
        """
        # pylint: disable=unused-argument
        super().__init__()

    @staticmethod
    def V_o(d: DutyCycle, V_s: float) -> float:
        """Calculates the V_o (average voltage) of a non-isolated flyback
        converter.

        :param d: Duty Cycle of the converter.
        :param V_s: Supply voltage of the converter.
        :return: The average voltage output of the converter.
        """
        return d * V_s / (1 - d)

    @staticmethod
    def V_s(d: DutyCycle, V_o: float) -> float:
        """Calculates the V_s (supply voltage) of a non-isolated flyback
        converter.

        :param d: Duty Cycle of the converter.
        :param V_o: Average voltage of the converter.
        :return: The supply voltage of the converter.
        """
        return V_o * (1 - d) / d


# Isolated
class DutyCycleIsolated(DutyCycle):
    """Duty cycle equation for an isolated flyback converter.

    Formula: d = \frac {N_{1}}{N_{2}} \frac {V_{o}}{V_{s} + V_{o}}
    """
    def __new__(cls, V_o: float, V_s: float,
                N1: float, N2: float):
        return super().__new__(cls, (V_o * N1) / ((V_s - V_o) * N2))

    def __init__(self, V_o: float, V_s: float,
                 N1: float, N2: float):
        """Duty cycle equation for an isolated flyback converter.

        :param V_o: Average voltage output of the converter.
        :param V_s: Supply voltage of the converter.
        :param N1: The number of primary windings in the transformer.
        :param N2: The number of secondary windings in the transformer.
        """
        # pylint: disable=unused-argument
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
        """Calculates the V_o (average voltage) of the flyback converter.

        :param d: Duty Cycle of the forward converter.
        :param V_s: Supply voltage of the forward converter.
        :param N1: The number of primary windings in the transformer.
        :param N2: The number of secondary windings in the transformer.
        :return: Average voltage output of the converter.
        """
        return N2 / N1 * d * V_s / (1 - d)

    @staticmethod
    def V_s(d: DutyCycle, V_o: float, N1: float, N2: float) -> float:
        """Calculates the V_s (supply voltage) of the flyback converter.

        :param d: Duty Cycle of the converter.
        :param V_o: Average voltage of the converter.
        :param N1: The number of primary windings in the transformer.
        :param N2: The number of secondary windings in the transformer.
        :return: The supply voltage of the converter.
        """
        return N1 / N2 * V_o * (1 - d) / d

    @staticmethod
    def turns_ratio(d: DutyCycle, V_o: float, V_s: float) -> float:
        """Calculates the turns ratio of the flyback converter.

        :param d: Duty Cycle of the forward converter.
        :param V_s: Supply voltage of the forward converter.
        :param V_o: Average voltage output of the converter.
        """
        return d * V_s / V_o / (1 - d)
