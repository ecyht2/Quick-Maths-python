"""Functions and Equations taught in EEEE1028 Power and Energy
module."""
import cmath
from math import acos, cos, pi, radians, sin, sqrt, tan
from typing import Union

from eeepy.utils.helper import product
from eeepy.utils.units import angular_frequency


# Circuit Analysis
def current(Q: float, t: float) -> float:
    """
    Calculates the current from charge over time
    """
    return Q / t


def voltage(W: float, Q: float) -> float:
    """
    Calculates the voltage with energy (W) and charge (Q)
    """
    return W / Q


def electrical_power(V: float = 0, I: float = 0, R: float = 0) -> float:
    """
    Calculates the of an electrical component
    """
    P = 0
    if abs(V) == 0 and abs(I) > 0 and abs(R) > 0:
        P = I**2 * R
    elif abs(I) == 0 and abs(V) > 0 and abs(R) > 0:
        P = V**2 / R
    elif abs(R) == 0 and abs(V) > 0 and abs(I) > 0:
        P = V * I
    else:
        raise ValueError("At least 2 of the parameters must be provided")
    return P


def resistance(rho: float, l: float, A: float) -> float:
    """
    Calculates the resistance of a given material
    """
    return rho * l / A


def series_res(resistors: list or tuple, *res) -> float:
    """
    Calculates the resistance in series
    """
    if isinstance(resistors, Union[list, tuple]):
        R = sum(resistors)
    else:
        R = resistors + sum(res)
    return R


def par_res(resistors: list or tuple, *res) -> float:
    """Calculates the resistance in parallel."""
    if isinstance(resistors, Union[list, tuple]):
        flipped = []
        for i in resistors:
            flipped.append(1 / i)
        total = sum(flipped)**-1
    else:
        flipped = [1 / resistors]
        for i in res:
            flipped.append(1 / i)
        total = sum(flipped)**-1
    return total


# Ohm's Law
def ohms_V(I: float, R: float) -> float:
    """Calculates the voltage using Ohms's Law."""
    return I * R


def ohms_I(V: float, R: float) -> float:
    """Calculates the current using Ohms's Law."""
    return V / R


def ohms_R(V: float, I: float) -> float:
    """Calculates the resistance using Ohms's Law."""
    return V / I


# Divider Rules
def voltage_divider(V_total: float, R_out: float, R_rest: float,
                    *res) -> float:
    """Calculates the voltage output using voltage divider rule."""
    v_out = 0

    v_out = V_total * R_out
    if isinstance(R_rest, Union[list, tuple]):
        v_out = v_out / (R_out + sum(R_rest))
    else:
        v_out = v_out / (R_out + R_rest + sum(res))

    return v_out


def current_divider_rt(I_total: float, R_out: float, R_t: float):
    """Calculates the voltage output using voltage divider rule using Rt."""
    I_out = 0

    I_out = R_t * I_total / R_out
    return I_out


def current_divider_r_split(I_total: float, R_out: float, R_rest: float, *res):
    """Calculates the current output using current divider rule."""
    I_out = 0

    if isinstance(R_rest, Union[list, tuple]):
        I_out = I_total * product(R_rest) / (R_out + sum(R_rest))
    else:
        I_out = I_total * R_rest * product(res) / (R_out + R_rest + sum(res))

    return I_out


# Capacitor and Inductor
# Capacitor
def capacitance(Q: float, V: float) -> float:
    """
    Calculate a capacitance of a capacitor
    """
    return Q / V


def capacitor_current(C: float = 0, dv: float = 0, dQ: float = 0) -> float:
    """
    Calculates the current flow of a capacitor
    """
    I = 0
    if C > 0 and dv > 0:
        I = C * dv
    elif dQ > 0:
        I = dQ
    else:
        raise ValueError("No (C and dv) or dQ given")
    return I


def capacitor_power(V: float, C: float = 0, Q: float = 0) -> float:
    """
    Calculates the energy supplied for a given capacitor
    """
    P = 0
    if C > 0:
        P = 0.5 * C * V**2
    elif Q > 0:
        P = 0.5 * Q * V
    else:
        raise ValueError("No C or Q given")
    return P


def capacitor_charge(C: float, V: float) -> float:
    """
    Calculates the charge stored in a capacitor
    """
    return C * V


# Alias Functions
capacitance_series = par_res
capacitance_parallel = series_res
capacitor_voltage_distribution = voltage_divider


def RC_time_constant(R: float, C: float):
    """
    Calculates the time constant (τ) of an RC circuit
    """
    return R * C


# Inductor
def inductor_voltage(L: float, di: float) -> float:
    """
    Calculates the emf induced by and inductor
    """
    return L * di


def inductor_power(L: float, I: float) -> float:
    """
    Calculates the energy supplied for a given capacitor
    """
    return 0.5 * L * I**2


inductance_series = series_res
inductance_parallel = par_res


def RL_time_constant(R: float, L: float):
    """
    Calculates the time constant (τ) of an RC circuit
    """
    return L / R


# Network Theorem
# Source Transformation
def source_transform_V_to_I(Vs: float, Rs: float) -> float:
    """
    Converts a voltage source in series into current source in parallel
    """
    return ohms_I(Vs, Rs)


def source_transform_I_to_V(Is: float, Rs: float) -> float:
    """
    Converts a current source in parallel into voltage source in series
    """
    return ohms_V(Is, Rs)


# Delta-Star Transformation
def delta_to_star(Rbranch1: float, Rbranch2: float, Ropposite: float) -> float:
    """
    Converts a resistor in delta into star equivalent resistance
    """
    return Rbranch1 * Rbranch2 / (Rbranch1 + Rbranch2 + Ropposite)


def star_to_delta(Rnode1: float, Rnode2: float, Ropposite: float) -> float:
    """
    Converts a resistor in star into delta  equivalent resistance
    """
    return Rnode1 + Rnode2 + (Rnode1 * Rnode2) / Ropposite


# Impedence
def rct_cap(C: float, f: float = 0, omega: float = 0, T: float = 0) -> float:
    """
    Calculates the reactance of a capacitor
    """
    X = 1 / C
    if omega > 0:
        X *= 1 / omega
    elif f > 0:
        X *= 1 / (2 * pi * f)
    elif T > 0:
        X *= T / (2 * pi)
    else:
        raise ValueError("No omega or f or T given")
    return X


def rct_cap_rev(X: float, f: float = 0, omega: float = 0,
                T: float = 0) -> float:
    """
    Calculates the capacitance of a capacitor from the reactance
    """
    C = 1 / X
    if f > 0:
        C *= 1 / (2 * pi * f)
    elif omega > 0:
        C *= 1 / (omega)
    elif T > 0:
        C *= T / (2 * pi)
    return C


def rct_ind(L: float, f: float = 0, omega: float = 0, T: float = 0) -> float:
    """
    Calculates the reactance of an inductor
    """
    X = L
    if omega > 0:
        X *= omega
    elif f > 0:
        X *= 2 * pi * f
    elif T > 0:
        X *= 2 * pi / T
    else:
        raise ValueError("No omega or f or T given")
    return X


def rct_ind_rev(X: float, f: float = 0, omega: float = 0,
                T: float = 0) -> float:
    """
    Calculates the inductance of an inductor from the reactance
    """
    L = X
    if f > 0:
        L /= 2 * pi * f
    elif omega > 0:
        L /= omega
    elif T > 0:
        L /= 2 * pi / T
    return L


def imp_cap(C: float, f: float = 0, omega: float = 0, T: float = 0) -> complex:
    """Calculates the impedance of a capacitor."""
    Z = rct_cap(C, f, omega, T) / 1j
    return Z


def imp_ind(L: float, f: float = 0, omega: float = 0, T: float = 0) -> complex:
    """Calculates the impedance of an inductor."""
    Z = rct_ind(L, f, omega, T) * 1j
    return Z


def imp_cap_rev(Z: complex, f: float = 0,
                omega: float = 0, T: float = 0) -> float:
    """Calculates the capacitance from a purely capacitive impedance."""
    X = Z * 1j
    return rct_cap_rev(X.real, f, omega, T)


def imp_ind_rev(Z: complex, f: float = 0,
                omega: float = 0, T: float = 0) -> float:
    """Calculates the capacitance from a purely inductive impedance."""
    X = Z / 1j
    return rct_ind_rev(X.real, f, omega, T)


# AC Power
def active_power(rms_val: tuple = (0, 0), S: float = 0,
                 phi: float = 0, PF: float = 0,
                 radian: bool = True) -> float:
    """
    Calculates the active power of an AC circuit
    """
    Irms = rms_val[0]
    Vrms = rms_val[1]
    if S > 0:
        pass
    elif Irms > 0 and Vrms > 0:
        S = Irms * Vrms
    else:
        raise ValueError("No (Vrms and Irms) or S given")
    P = S
    if abs(phi) > 0:
        if not radian:
            phi = radians(phi)
        P *= cos(phi)
    elif PF > 0:
        P *= PF
    else:
        raise ValueError("No phi or PF given")
    return P


def apparent_power(Vrms: float = 0, Irms: float = 0,
                   P: float = 0, Q: float = 0) -> float:
    """
    Calculates the apparent power of an AC circuit
    """
    S = 0
    if Vrms > 0 and Irms > 0:
        S = Vrms * Irms
    elif P > 0 and Q > 0:
        S = sqrt(P**2 + Q**2)
    else:
        raise ValueError("No (Vrms and Irms) or (P and Q) given")
    return S


def power_factor(P: float = 0, S: float = 0,
                 theta: float = 0, radian: bool = True) -> float:
    """
    Calculates the power factor of an AC circuit
    """
    if not radian:
        theta = radians(theta)
    PF = 0
    if abs(P) > 0 and abs(S) > 0:
        PF = abs(P) / abs(S)
    elif theta != 0:
        PF = cos(theta)
    return PF


def reactive_power(phi: float,
                   Vrms: float = 0, Irms: float = 0, S: float = 0,
                   radian: bool = True) -> float:
    """
    Calculates the reactive power of an AC circuit
    """
    if S > 0:
        pass
    elif Irms > 0 and Vrms > 0:
        S = Irms * Vrms
    else:
        raise ValueError("No (Vrms and Irms) or S given")
    Q = S
    if not radian:
        phi = radians(phi)
    Q *= sin(phi)
    return Q


def complex_power(phi: float,
                  Vrms: float = 0, Irms: float = 0, S: float = 0,
                  radian: bool = True) -> complex:
    """
    Calculates the total complex power of an AC circuit
    """
    if S > 0:
        pass
    elif Irms > 0 and Vrms > 0:
        S = Irms * Vrms
    else:
        raise ValueError("No (Vrms and Irms) or S given")
    Sstar = S
    if not radian:
        phi = radians(phi)
    Sstar *= cos(phi) + sin(phi) * 1j
    return Sstar


def power_factor_correction(P: float, theta: tuple,
                            Vrms: float, f: float = 0,
                            omega: float = 0) -> float:
    """
    Calculates the capacitance needed to increase the power factor of a circuit
    without altering the voltage or current to the original load
    """
    if omega > 0:
        pass
    elif f > 0:
        omega = 2 * pi * f
    else:
        raise ValueError("No f or omega or T given")
    Qc = abs(P) * (tan(theta[0]) - tan(theta[1]))
    C = Qc / (omega * abs(Vrms)**2)
    return C


# Transient Analysis
def general_approach(final: float, initial: float, tau: float) -> str:
    """
    The transient state of and RC or RL circuit using general approach
    """
    if final == 0:
        eq = f"{initial - final}e^(-{1/tau}t)"
    else:
        eq = f"{final} + {initial - final}e^(-{1/tau}t)"
    return eq


# Three phase system
# Balanced System
# Current and Voltage
def line_voltage_star_polar(Vp: float, phase: float = None,
                            radian: bool = True) -> float:
    """Calculates the line voltage of a balanced three phase system in star\
configuration.
    """
    VL = sqrt(3) * Vp
    if phase is not None:
        if radian:
            phaseL = phase + pi / 6
        else:
            phaseL = phase + 30
        ret = (VL, phaseL)
    else:
        ret = VL
    return ret


def line_voltage_star_polar_rev(VL: float, phase: float = None,
                                radian: bool = True) -> float:
    """Calculates the phase voltage of a balanced three phase system in star\
configuration.
    """
    Vp = VL / sqrt(3)
    if phase is not None:
        if radian:
            phaseP = phase - pi / 6
        else:
            phaseP = phase - 30
        ret = (Vp, phaseP)
    else:
        ret = Vp
    return ret


def line_current_delta_polar(Ip: float, phase: float = None,
                             radian: bool = True) -> float:
    """Calculates the line current of a balanced three phase system in delta\
configuration.
    """
    IL = sqrt(3) * Ip
    if phase is not None:
        if radian:
            phaseL = phase - pi / 6
        else:
            phaseL = phase - 30
        ret = (IL, phaseL)
    else:
        ret = IL
    return ret


def line_current_delta_polar_rev(IL: float, phase: float = None,
                                 radian: bool = True) -> float:
    """Calculates the phase current of a balanced three phase system in delta\
configuration.
    """
    Ip = IL / sqrt(3)
    if phase is not None:
        if radian:
            phaseP = phase + pi / 6
        else:
            phaseP = phase + 30
        ret = (Ip, phaseP)
    else:
        ret = Ip
    return ret


def line_voltage_star(Vp: complex) -> complex:
    """Calculates the line voltage of a balanced three phase system in star\
configuration.
    """
    VL = cmath.rect(sqrt(3), pi / 6) * Vp
    return VL


def line_voltage_star_rev(VL: complex) -> complex:
    """Calculates the phase voltage of a balanced three phase system in star\
configuration.
    """
    Vp = VL / cmath.rect(sqrt(3), pi / 6)
    return Vp


def line_current_delta(Ip: complex) -> complex:
    """Calculates the line current of a balanced three phase system in delta\
configuration.
    """
    IL = cmath.rect(sqrt(3), - pi / 6) * Ip
    return IL


def line_current_delta_rev(IL: complex) -> complex:
    """Calculates the phase current of a balanced three phase system in delta\
configuration.
    """
    Ip = IL / cmath.rect(sqrt(3), - pi / 6)
    return Ip


# Power
def balanced_three_phase_total_power(phi: float,
                                     Vph: float = 0, Iph: float = 0,
                                     VL: float = 0, IL: float = 0) -> float:
    """
    Calculates the total power of a balanced three phase system
    """
    if Vph > 0 and Iph > 0:
        P = 3 * Vph * Iph * cos(phi)
    elif VL > 0 and IL > 0:
        P = sqrt(3) * Vph * Iph * cos(phi)
    else:
        raise ValueError("No (Vph and Iph) or (VL and IL) given")
    return P


def balanced_three_phase_total_apparent_power(Vph: float = 0, Iph: float = 0,
                                              VL: float = 0,
                                              IL: float = 0) -> float:
    """
    Calculates the total apparent power of a balanced three phase system
    """
    if Vph > 0 and Iph > 0:
        S = 3 * Vph * Iph
    elif VL > 0 and IL > 0:
        S = sqrt(3) * Vph * Iph
    else:
        raise ValueError("No (Vph and Iph) or (VL and IL) given")
    return S


# Transformers
# Transformers parameters
def turns_ratio(N1: int = 0, N2: int = 0,
                K: float = 0, a: float = 0) -> float:
    """
    Calculates the turns ratio of a transformer
    """
    if N1 > 0 and N2 > 0:
        a = N1 / N2
    elif a > 0:
        pass
    elif K > 0:
        a = 1 / K
    else:
        raise ValueError("No (N1 and N2) or a or K given")
    return a


def voltage_transformation_ratio(N1: int = 0, N2: int = 0,
                                 a: float = 0, K: float = 0) -> float:
    """
    Calculates the voltage transformation ratio of a transformer
    """
    if N1 > 0 and N2 > 0:
        K = N2 / N1
    elif a > 0:
        K = 1 / a
    elif K > 0:
        pass
    else:
        raise ValueError("No (N1 and N2) or a or K given")
    return K


def transformer_type(N1: int = 0, N2: int = 0,
                     a: float = 0, K: float = 0):
    """
    Determines the type of transformer a transformer is
    """
    K = voltage_transformation_ratio(N1, N2, a, K)
    if K > 1:
        tType = "Step-up Transformer"
    elif K < 1:
        tType = "Step-down Transformer"
    else:
        tType = "Isolation Transformer"
    return tType


# Ideal
def primary_voltage(V2: float, N1: int = 0, N2: int = 0,
                    a: float = 0, K: float = 0) -> float:
    """Calculates the voltage induced on the primary coil of an ideal or\
almost ideal transformer.
    """
    a = turns_ratio(N1, N2, K, a)
    return V2 * a


def primary_current_ideal(I2: float, N1: int = 0, N2: int = 0,
                          a: float = 0, K: float = 0) -> float:
    """
    Calculates the current induced on the primary coil of an ideal transformer
    """
    a = turns_ratio(N1, N2, K, a)
    return I2 / a


def primary_impedence(Z2: complex, N1: int = 0, N2: int = 0,
                      a: float = 0, K: float = 0) -> float:
    """
    Calculates the secondary impedence reflected on the primary
    """
    a = turns_ratio(N1, N2, K, a)
    return a**2 * Z2


def secondary_voltage(V1: float, N1: int = 0, N2: int = 0,
                      a: float = 0, K: float = 0) -> float:
    """Calculates the voltage induced on the secondary coil of an ideal or\
almost ideal transformer.
    """
    a = turns_ratio(N1, N2, K, a)
    return V1 / a


def secondary_current_ideal(I1: float, N1: int = 0, N2: int = 0,
                            a: float = 0, K: float = 0) -> float:
    """Calculates the current induced on the secondary coil of an ideal\
transformer.
    """
    a = turns_ratio(N1, N2, K, a)
    return I1 * a


def secondary_impedence(Z1: complex, N1: int = 0, N2: int = 0,
                        a: float = 0, K: float = 0) -> float:
    """
    Calculates the primary impedence reflected on the secondary
    """
    a = turns_ratio(N1, N2, K, a)
    return Z1 / a**2


# Almost ideal
def primary_current_almost_ideal(I2: float, Im: float = 0,
                                 N1: int = 0, N2: int = 0,
                                 a: float = 0, K: float = 0) -> float:
    """Calculates the current induced on the primary coil of an almost ideal\
transformer.
    """
    # pylint: disable = too-many-arguments
    a = turns_ratio(N1, N2, K, a)
    return Im + (I2 / a)


def secondary_current_almost_ideal(I1: float, Im: float = 0,
                                   N1: int = 0, N2: int = 0,
                                   a: float = 0, K: float = 0) -> float:
    """Calculates the current induced on the secondary coil of an almost ideal\
transformer.
    """
    # pylint: disable = too-many-arguments
    a = turns_ratio(N1, N2, K, a)
    return (I1 - Im) * a


# Efficiency and voltage regulation
def voltage_regulation(V2no: float, V2full: float) -> float:
    """
    Calculates the voltage regulation of a transformer
    """
    return (V2no - V2full) / V2full


def transformer_efficiency_x_load(Pout: float, Pcore: float,
                                  Pcopper: float, x: float = 1,
                                  formatted: bool = False) -> float:
    """
    Calculates the effieciency of a transformer at x*100% load
    """
    if x > 1 or x < 0:
        raise ValueError("x must be 0 <= x <= 1")
    ret = (x * Pout) / (x * Pout + Pcore + x**2 * Pcopper)
    if formatted:
        ret = str(ret * 100) + "%"
    return ret


# SC and OC test
def OC_test(P: float, I: float, V: float,
            f: float = 0, omega: float = 0, T: float = 0,
            reactance: bool = True) -> tuple:
    """
    Calculates the Rc and (Xm or Lm) of a transformer

    Parameters
    ----------
    P
        Power of the transformer
    I
        Magnitude of current
    V
        Magnitude of voltage
    f (Not needed if omega or T given)
        Frequency of circuit
    omega (Not needed if T or f given)
        Angular Frequency of circuit
    T (Not needed if omega or f given)
        Period of circuit
    reactance (True by default)
        If set to True the reactance will be given otherwise the inductance\
will be given

    Returns
    -------
    Tuple
        (Rc, Xm) if reactance is True
        (Rc, Lm) if reactance is False
    """
    # pylint: disable = too-many-arguments
    if omega > 0:
        pass
    elif f > 0 or T > 0:
        omega = angular_frequency(f, T)
    else:
        raise ValueError("No f or omega or T given")
    # Calculating PF and Θ
    PF = power_factor(P, V * I)
    theta = acos(PF)

    # Calculating Losses
    try:
        Rc = V / (I * PF)
    except ZeroDivisionError:
        Rc = 0
    try:
        Lm = V / (omega * I * sin(theta))
    except ZeroDivisionError:
        Lm = 0

    # Setting return
    if reactance:
        ret = (Rc, Lm * omega)
    else:
        ret = (Rc, Lm)
    return ret


def SC_test(P: float, I: float, V: float,
            f: float = 0, omega: float = 0, T: float = 0,
            reactance: bool = True) -> tuple:
    """
    Calculates the R and (X or L) of a transformer

    Parameters
    ----------
    P
        Power of the transformer
    I
        Magnitude of current
    V
        Magnitude of voltage
    f (Not needed if omega or T given)
        Frequency of circuit
    omega (Not needed if T or f given)
        Angular Frequency of circuit
    T (Not needed if omega or f given)
        Period of circuit
    reactance (True by default)
        If set to True the reactance will be given otherwise the inductance\
will be given

    Returns
    -------
    Tuple
        (R, X) if reactance is True
        (R, L) if reactance is False
    """
    # pylint: disable = too-many-arguments
    if omega > 0:
        pass
    elif f > 0 or T > 0:
        omega = angular_frequency(f, T)
    else:
        raise ValueError("No f or omega or T given")
    # Calculating PF and Θ
    PF = power_factor(P, V * I)
    theta = acos(PF)

    # Calculating Losses
    try:
        R = V * PF / I
    except ZeroDivisionError:
        R = 0
    try:
        L = V * sin(theta) / (omega * I)
    except ZeroDivisionError:
        L = 0

    # Setting return
    if reactance:
        ret = (R, L * omega)
    else:
        ret = (R, L)
    return ret
