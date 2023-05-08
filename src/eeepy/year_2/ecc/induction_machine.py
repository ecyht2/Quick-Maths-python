#!/usr/bin/env python3
"""Functions and equations related to induction machines."""


class Slip:
    """The slip of an induction machine.

    s = (ɷ_s - ɷ_r) / ɷ_s
    ----------------------
    s is the slip.
    ɷ_s is the synchronous speed in rad/s.
    ɷ_r is the rotor (mechanical) speed in rad/s.
    """
    @staticmethod
    def slip(w_s: float, w_r: float) -> float:
        """Calculates the slip.

        :param ɷ_s: The synchronous speed in rad/s.
        :param ɷ_r: The rotor (mechanical) speed in rad/s.
        :returns: The slip.
        """
        return (w_s - w_r) / w_s

    @staticmethod
    def rotor_speed(s: float, w_s: float) -> float:
        """Calculates the rotor speed in rad/s.

        :param s: The slip.
        :param ɷ_s: The synchronous speed in rad/s.
        :returns: The rotor (mechanical) speed in rad/s.
        """
        return (1 - s) * w_s

    @staticmethod
    def synchronous_speed(s: float, w_r: float) -> float:
        """Calculates the synchronous speed in rad/s.

        :param s: The slip.
        :param ɷ_r: The rotor (mechanical) speed in rad/s.
        :returns: The synchronous speed in rad/s.
        """
        return w_r / (1 - s)


class CopperLoss:
    """Copper loss in an induction machine.

    P_SCL = 3 I_1^2 R_1
    -------------------
    P_SCL The copper loss at the stator.
    I_1 is the current at stator.
    R_1 is the resistance at the stator.

    P_RCL = 3 I_2^2 R_2
    -------------------
    P_RCL is the copper loss at the rotor.
    I_2 is the current at rotor.
    R_2 is the resistance at the rotor.
    """
    @staticmethod
    def power_scl(I_1: float, R_1: float) -> float:
        """Copper loss at the stator.

        :param I_1: The current at stator.
        :param R_1: The resistance at the stator.
        :returns: The copper loss at the stator.
        """
        return 3 * I_1**2 * R_1

    @staticmethod
    def power_rcl(I_2: float, R_2: float) -> float:
        """Copper loss at the rotor.

        :param I_2: The Current at rotor.
        :param R_2: The Resistance at the rotor.
        :returns: The copper loss at the rotor.
        """
        return 3 * I_2**2 * R_2


class AirGapPower:
    """Power at the air gap of an induction machine.

    P_AG = P_in - P_SCL - P_core
    ----------------------------
    P_AG is the air gap power.
    P_in is the input power.
    P_SCL is the stator copper loss.
    P_core is the core loss.

    P_AG = 3 I_2**2 * R_2 / s
    -------------------------
    P_AG is the air gap power.
    I_2 is the current at rotor.
    R_2 is the resistance at the rotor.
    s is the slip of the induction machine.
    """
    @staticmethod
    def from_input_power(P_in: float, P_SCL: float, P_core: float) -> float:
        """Calculates the air gap power from input power.

        :param P_in: The input power.
        :param P_SCL: The stator copper loss.
        :param P_core: The core loss.
        :returns: The air gap power.
        """
        return P_in - P_SCL - P_core

    @staticmethod
    def from_rotor(I_2: float, R_2: float, s: float) -> float:
        """Calculates the air gap power from rotor parameters.

        :param I_2: The current at rotor.
        :param R_2: The resistance at the rotor.
        :param s: The slip of the induction machine.
        :returns: The air gap power.
        """
        return 3 * I_2**2 * R_2 / s


class PowerConverted:
    """The power converted in the induction machine.

    P_conv = P_AG - P_RCL
    ---------------------
    P_conv is the power converted.
    P_AG is the air gap power.
    P_RCL is the rotor copper loss.

    P_conv = P_AG (1 - s)
    ---------------------
    P_conv is the power converted.
    P_AG is the air gap power.
    s is the slip.
    """
    @staticmethod
    def from_copper_loss(P_AG: float, P_RCL: float) -> float:
        """Calculates the converted power from rotor copper loss.

        :param P_AG: The air gap power.
        :param P_RCL: The rotor copper loss.
        :returns: The power converted.
        """
        return P_AG - P_RCL

    @staticmethod
    def from_slip(P_AG: float, s: float) -> float:
        """Calculates the converted power from slip.

        :param P_AG: The air gap power.
        :param s: The slip.
        :returns: The power converted.
        """
        return P_AG * (1 - s)


class InducedToruqe:
    """The induced torque by the induction machine.

    T_ind = P_conv / ɷ_m
    --------------------
    T_ind is the induced torque.
    P_conv is the power converted.
    ɷ_m is the mechanical speed in rad/s.

    T_ind = P_AG / w_s
    --------------------
    T_ind is the induced torque.
    P_AG is the air gap power.
    ɷ_s is the synchronous speed in rad/s.
    """
    @staticmethod
    def from_power_converted(P_conv: float, w_m: float) -> float:
        """Calculates the induced torque using power converted.

        :param P_conv: The power converted.
        :param ɷ_m: The mechanical speed in rad/s.
        :returns: The induced torque.
        """
        return P_conv / w_m

    @staticmethod
    def from_air_gap_power(P_AG: float, w_s: float) -> float:
        """Calculates the induced torque using air gap power.

        :param P_AG: The air gap power.
        :param ɷ_s: The synchronous speed in rad/s.
        :returns: The induced torque.
        """
        return P_AG / w_s


class OutputPower:
    """Output power of an induction machine.

    P_out = P_conv - P_misc - P_fw
    ------------------------------
    P_out is the output power.
    P_conv is the power converted.
    P_misc is the stray power losses.
    P_fw is the friction and windage losses.

    P_out = T_load ɷ_m
    ------------------
    P_out is the output power.
    T_load is the load torque.
    ɷ_m is the mechanical speed in rad/s.
    """
    @staticmethod
    def output_power(P_conv: float, P_misc: float, P_fw: float) -> float:
        """Calculates the output power.

        :param P_conv: The power converted.
        :param P_misc: The stray power losses.
        :param P_fw: The friction and windage losses.
        :returns: The output power.
        """
        return P_conv - P_misc - P_fw

    @staticmethod
    def load_torque(P_out: float, w_m: float) -> float:
        """Calculates the load torque.

        :param P_out: The output power.
        :param ɷ_m: The mechanical speed in radians/s.
        :returns: The load torque.
        """
        return P_out / w_m
