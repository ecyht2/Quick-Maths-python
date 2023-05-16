#!/usr/bin/env python3
"""Functions and equations relating to magnetic materials and characteristics.
"""
from typing import Optional, Tuple


class EddyCurrentLoss:
    """Eddy current loss.

    Formula: :math:`P_e = K_e B_{\\max}^2 f^2 t^2 V`

    :math:`P_e` is the eddy current loss.

    :math:`K_e` is the material and lamination thickness dependant constant.

    :math:`B_{\\max}` is the maximum flux density.

    :math:`f` is the frequency of the exitation field.

    :math:`t` is the thickness of the core.

    :math:`V` is the volume of the magnetic material.
    """
    @staticmethod
    def power_loss(K_e: float, B_max: float, f: float,
                   t: float, V: float) -> float:
        """Calculates the Eddy Current Loss.

        :param K_e: The eddy current constant value.
        :param B_max: The maximum flux density.
        :param f: The frequency of the exitation field.
        :param t: The thickness of the core.
        :param V: The volume of the core.
        :returns: The eddy current loss.
        """
        return K_e * B_max**2 * f**2 * t**2 * V

    @staticmethod
    def eddy_current_constant(P_e: float, B_max: float, f: float,
                              t: float, V: float) -> float:
        """Calculates the Eddy Current Loss Constant.

        :param P_e: The eddy current loss.
        :param B_max: The maximum flux density.
        :param f: The frequency of the exitation field.
        :param t: The thickness of the core.
        :param V: The volume of the core.
        :returns: The eddy current constant value.
        """
        return P_e / (B_max**2 * f**2 * t**2 * V)

    @staticmethod
    def new_loss(Pold: float,
                 K_e: Optional[Tuple[float, float]] = None,
                 B_max: Optional[Tuple[float, float]] = None,
                 f: Optional[Tuple[float, float]] = None,
                 t: Optional[Tuple[float, float]] = None,
                 V: Optional[Tuple[float, float]] = None) -> float:
        """Calculates the new eddy current loss when the a value changed.

        Each of the values are represented as a tuple of 2 floats. The first
        value of the tuple is the old value and the second value is the value
        it changed to.

        :param Pold: The old eddy current loss.
        :param K_e: The old and new eddy current constant values.
        :param B_max: The old and new maximum flux densities.
        :param f: The old and new frequency.
        :param t: The old and new thickness of the core.
        :param V: The old and new volume of the core.
        :returns: The new eddy current loss.

        Example
        -------
        Find the new loss when the frequency changed to 50 to 60 Hz and the
        maximum flux density changed from 1 to 0.5 T.

        >>> from eeepy.year_2.ecc.magnetic_materials import EddyCurrentLoss
        >>> f = (50.0, 60.0)  # Frequency changes
        >>> B = (1.0, 0.5)    # Flux density changes
        >>> EddyCurrentLoss.new_loss(69, f=f, B_max=B)
        24.84
        """
        ratio = 1
        if K_e is not None:
            ratio *= K_e[1] / K_e[0]
        if B_max is not None:
            ratio *= B_max[1]**2 / B_max[0]**2
        if f is not None:
            ratio *= f[1]**2 / f[0]**2
        if t is not None:
            ratio *= t[1]**2 / t[0]**2
        if V is not None:
            ratio *= V[1] / V[0]
        return Pold * ratio


class HysterisisLoss:
    """Hysterisis Loss.

    Formula: :math:`P_h = K_h f V B_{\\max}^{\\eta}`

    :math:`P_h` is the hysterisis loss.

    :math:`K_h` is the material dependant constant.

    :math:`f` is the frequency of the exitation field.

    :math:`V` is the volume of the magnetic material.

    :math:`B_{max}` is the maximum flux density.

    :math:`\\eta` is the Steinmetz's exponent constant usually between 1.5 and
    2.5 depending on the material.
    """
    @staticmethod
    def power_loss(K_h: float, f: float, V: float,
                   B_max: float, eta: float) -> float:
        """Calculates the Hysterisis Loss.

        :param K_h: The hysterisis constant values.
        :param f: The frequency of the exitation field.
        :param V: The volume of the core.
        :param B_max: The maximum flux density.
        :param eta: The Steinmetz's exponent constant.
        :returns: The hysterisis loss.
        """
        return K_h * f * V * B_max**eta

    @staticmethod
    def hysterisis_constant(P_h: float, f: float, V: float,
                            B_max: float, eta: float) -> float:
        """Calculates the Hysterisis loss material dependant constant.

        :param P_h: The hysterisis loss.
        :param f: The frequency of the exitation field.
        :param V: The volume of the core.
        :param B_max: The maximum flux density.
        :param eta: The Steinmetz's exponent constant.
        :returns: The hysterisis constant values.
        """
        return P_h / (f * V * B_max**eta)

    @staticmethod
    def new_loss(Pold: float,
                 K_h: Optional[Tuple[float, float]] = None,
                 f: Optional[Tuple[float, float]] = None,
                 V: Optional[Tuple[float, float]] = None,
                 B_max: Optional[Tuple[float, float]] = None,
                 eta: Optional[Tuple[float, float]] = None) -> float:
        """Calculates the new eddy current loss when the a value changed.

        Each of the values are represented as a tuple of 2 floats. The first
        value of the tuple is the old value and the second value is the value
        it changed to.

        .. note::
            B_max and eta must be specified together.

        :param Pold: The old hysterisis loss.
        :param K_h: The old and new hysterisis constant values.
        :param f: The old and new frequency.
        :param V: The old and new volume of the core.
        :param B_max: The old and new maximum flux densitys.
        :param eta: The old and new Steinmetz's exponent constants.
        :returns: The new hysterisis loss.

        Example
        -------
        Find the new loss when the frequency changed to 50 to 60 Hz and the
        maximum flux density changed from 1 to 0.5 T at 2.5 Steinmetz's
        exponent.

        >>> from eeepy.year_2.ecc.magnetic_materials import HysterisisLoss
        >>> f = (50.0, 60.0)  # Frequency changes
        >>> B = (1.0, 0.5)    # Flux density changes
        >>> eta = (2.5, 2.5)    # Steinmet's exponent
        >>> HysterisisLoss.new_loss(69, f=f, B_max=B, eta=eta)
        14.637110370561533
        """
        ratio = 1
        if K_h is not None:
            ratio *= K_h[1] / K_h[0]
        if (B_max is not None) and (eta is not None):
            ratio *= B_max[1]**eta[1] / B_max[0]**eta[0]
        elif any(((B_max is not None) and (eta is None),
                 (B_max is None) and (eta is not None))):
            raise ValueError("Both B_max and eta must be specified together.")
        if f is not None:
            ratio *= f[1] / f[0]
        if V is not None:
            ratio *= V[1] / V[0]
        return Pold * ratio
