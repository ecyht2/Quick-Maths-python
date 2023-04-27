#!/usr/bin/env python3
"""Functions and equations relating to magnetic materials and characteristics.
"""


class EddyCurrentLoss:
    """Eddy current loss.

    P_e = K_e B_max^2 f^2 t^2 V
    """
    @staticmethod
    def power_loss(K_e: float, B_max: float, f: float,
                   t: float, V: float) -> float:
        """Eddy Current Loss."""
        return K_e * B_max**2 * f**2 * t**2 * V

    @staticmethod
    def eddy_current_constant(P_e: float, B_max: float, f: float,
                              t: float, V: float) -> float:
        """Eddy Current Loss Constant."""
        return P_e / (B_max**2 * f**2 * t**2 * V)


class HysterisisLoss:
    """Hysterisis Loss.

    P_h = K_h f V B_max^η
    """
    @staticmethod
    def power_loss(K_h: float, f: float, V: float,
                   B_max: float, eta: float) -> float:
        """Hysterisis Loss."""
        return K_h * f * V * B_max**eta

    @staticmethod
    def hysterisis_constant(P_h: float, f: float, V: float,
                            B_max: float, eta: float) -> float:
        """Hysterisis loss material dependant constant."""
        return P_h / (f * V * B_max**eta)
