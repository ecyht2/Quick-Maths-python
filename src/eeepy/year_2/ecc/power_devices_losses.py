#!/usr/bin/env python3
"""Formulas used to calculate power losses of in a power electronics
converter circuit."""
import math


class ConductionLoss(float):
    """Conduction losses in power electronics."""
    @classmethod
    def diode(cls, V_FWD: float, I_D: float,
              R_s: float, I_D_RMS: float):
        """Calculates the conduction loss of a power diode.

        Formula: P_{loss} = F_{FWD}  Ī_{D} + R_{S} + (I_{D,RMS})^2

        :param: V_FWD: The forward voltage of the diode.
        :param I_D: The mean current experice by the diode.
        :param R_s: The series resistance of the dided.
        :param I_D_RMS: The RMS current of the diode.
        """
        return cls(V_FWD * I_D + R_s * I_D_RMS**2)

    @classmethod
    def mosfet(cls, R_DS, I_Q):
        """Calculates the conduction loss of a power MOSFET.

        Formula: P_{loss} = R_{DS(ON)} + (I_{Q,RMS})^2

        :param R_DS: The R_DS(ON) of the MOSFET.
        :param I_Q: The RMS current of the MOSFET.
        """
        return cls(R_DS * I_Q**2)


class SwitchingLoss(float):
    """Switching losses in forward, flyback and boos converters."""
    @classmethod
    def energy_turn_on(cls, I_t: float, V_s: float, t_TR_ON: float):
        """Calculates the energy loss during the turning on of the MOSFET.

        Formula: E_{ON} = \frac {1}{2} \\check{I}_{L} V_{S} t_{TR,ON}

        :param I_t: The minimum current of the inductor.
        :param V_s: The supply voltage.
        :param t_TR_ON: The turn on duration of the MOSFET.
        """
        return cls(0.5 * I_t * V_s * t_TR_ON)

    @classmethod
    def power_turn_on(cls, I_t: float, V_s: float, t_TR_ON: float,
                      f_sw: float):
        """Calculates the power loss during the turning on of the MOSFET.

        Formula: Power lost = \frac {1}{2} \\check{I}_{L} V_{S} t_{TR,ON}
        f_{sw}

        :param I_t: The minimum current of the inductor.
        :param V_s: The supply voltage.
        :param t_TR_ON: The turn on duration of the MOSFET.
        :param f_sw: The frequency of the converter.
        """
        return cls(0.5 * I_t * V_s * t_TR_ON * f_sw)

    @classmethod
    def power_turn_on_e(cls, E_ON: float, f_sw: float):
        """Calculates the power loss during the turning on of the MOSFET.

        Formula: Power lost = E_{ON} f_{sw}

        :param E_ON: The energy loss during the turning on of the MOSFET.
        :param f_sw: The frequency of the converter.
        """
        return cls(E_ON * f_sw)

    @classmethod
    def energy_turn_off(cls, I_p: float, V_s: float, t_TR_OFF: float):
        """Calculates the energy loss during the turning off of the MOSFET.

        Formula: E_{ON} = \frac {1}{2} \\hat{I}_{L} V_{S} t_{TR,OFF}

        :param I_p: The maximum current of the inductor.
        :param V_s: The supply voltage.
        :param t_TR_OFF: The turn off duration of the MOSFET.
        """
        return cls(0.5 * I_p * V_s * t_TR_OFF)

    @classmethod
    def power_turn_off(cls, I_p: float, V_s: float, t_TR_OFF: float,
                       f_sw: float):
        """Calculates the power loss during the turning off of the MOSFET.

        Formula: Power loss = \frac {1}{2} \\hat{I}_{L} V_{S} t_{TR,OFF}
        f_{sw}

        :param I_p: The maximum current of the inductor.
        :param V_s: The supply voltage.
        :param t_TR_OFF: The turn off duration of the MOSFET.
        :param f_sw: The frequency of the converter.
        """
        return cls(0.5 * I_p * V_s * t_TR_OFF * f_sw)

    @classmethod
    def power_turn_off_e(cls, E_OFF: float, f_sw: float):
        """Calculates the power loss during the turning off of the MOSFET.

        Formula: Power lost = E_{OFF} f_{sw}

        :param E_OFF: The energy loss during the turning off of the MOSFET.
        :param f_sw: The frequency of the converter.
        """
        return cls(E_OFF * f_sw)


def I_RMS(i1: float, i2: float, d: float):
    """Calculates the RMS value of any converter component.

    Formula: \\sqrt {\frac{i_{1}^{2} + i_{1} * i_{2} + i_{2}^{2}}{3}

    Parameters
    ----------
    i1: float
        The first current.
    i2: float
        The second current.
    d: float
        The ratio of how long the component is on for in each cycle
        e.g. The diode is only on when the transistor is of hence the d values
        needs to be (1 - dutycycle) where dutycycle is the dutycycle of the
        circuit.
    """
    ans = (i1**2 + i1 * i2 + i2**2) * d
    ans /= 3
    return math.sqrt(ans)
