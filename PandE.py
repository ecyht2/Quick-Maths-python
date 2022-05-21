import math
from Constants import *
from Maths import product
from helper import *

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
    if type(resistors) == list or type(resistors) == tuple:
        R = sum(resistors)
    else:
        R = resistors + sum(res)
    return R
def par_res(resistors: list or tuple, *res) -> float:
    """
    Calculates the resistance in parallel
    """
    if type(resistors) == list:
        flipped = []
        for i in resistors:
            flipped.append(1/i)
        total = sum(flipped)**-1
    else:
        flipped = [1/resistors]
        for i in res:
            flipped.append(1/i)
        total = sum(flipped)**-1
    return total

# Ohm's Law
def ohms_V(I: float, R: float) -> float:
    return I*R
def ohms_I(V: float, R: float) -> float:
    return V/R
def ohms_R(V: float, I: float) -> float:
    return V/I

# Divider Rules
def voltage_divider(V_total, R_out, R_rest, *res):
    v_out = 0

    v_out = V_total * R_out
    if type(R_rest) == list:
        v_out = v_out/(R_out + sum(R_rest))
    else:
        v_out = v_out/(R_out + R_rest + sum(res))

    return v_out
def current_divider_rt(I_total, R_out, R_t):
    I_out = 0

    I_out = R_t * I_total / R_out
    return I_out
def current_divider_r_split(I_total, R_out, R_rest, *res):
    I_out = 0

    if type(R_rest) == list:
        I_out = I_total*product(R_rest)/(R_out + sum(R_rest))
    else:
        I_out = I_total*R_rest*product(res)/(R_out + R_rest + sum(res))

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
    return Rside1 * Rside2 / (Rside1 + Rside2 + Ropposite)
def star_to_delta(Rnode1: float, Rnode2: float, Ropposite: float) -> float:
    """
    Converts a resistor in star into delta  equivalent resistance
    """
    return Rnode1 + Rnode2 + (Rnode1 * Rnode2) / Ropposite

# Impedence
def rct_cap(C: float, f: float = 0, omega: float = 0, T: float = 0):
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
def rct_cap_rev(X: float, f: float = 0, omega: float = 0, T: float = 0):
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
def rct_ind(L: float, f: float = 0, omega: float = 0, T: float = 0):
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
def rct_ind_rev(X: float, f: float = 0, omega: float = 0, T: float = 0):
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
    Z = rct_cap(C, f, omega, T) / 1j
    return Z
def imp_ind(L: float, f: float = 0, omega: float = 0, T: float = 0) -> complex:
    Z = rct_ind(L, f, omega, T) * 1j
    return Z
def imp_cap_rev(Z: complex, f: float = 0, omega: float = 0, T: float = 0) -> float:
    X = Z * 1j
    return rct_cap_rev(X.real, f, omega, T)
def imp_ind_rev(Z: complex, f: float = 0, omega: float = 0, T: float = 0) -> float:
    X = Z / 1j
    return rct_ind_rev(X.real, f, omega, T)

# AC Power
def active_power(Vrms: float = 0, Irms: float = 0, S: float = 0,
                 phi: float = 0, PF: float = 0,
                 radian: bool = True) -> float:
    """
    Calculates the active power of an AC circuit
    """
    if S > 0:
        S = S
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
    if P > 0 and S > 0:
        PF = P / S
    elif not theta == 0:
        PF = cos(theta)
    return PF
def reactive_power(phi: float,
                   Vrms: float = 0, Irms: float = 0, S: float = 0,
                   radian: bool = True) -> float:
    """
    Calculates the reactive power of an AC circuit
    """
    if S > 0:
        S = S
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
        S = S
    elif Irms > 0 and Vrms > 0:
        S = Irms * Vrms
    else:
        raise ValueError("No (Vrms and Irms) or S given")
    Sstar = S
    if not radian:
        phi = radians(phi)
    Sstar *= cos(phi) + sin(phi) * 1j
    return Sstar
def power_factor_correction(P: float, thetaNew: float, thetaOld: float, Vrms: float,
                            f: float = 0, omega: float = 0, T: float = 0) -> float:
    """
    Calculates the capacitance needed to increase the power factor of a circuit
    without altering the voltage or current to the original load
    """
    if omega > 0:
        omega = omega
    elif f > 0:
        omega = 2 * pi * f
    elif T > 0:
        omega = 2 * pi / T
    else:
        raise ValueError("No f or omega or T given")
    Qc = P * (tan(thetaOld) - tan(thetaNew))
    C = Qc / (omega * Vrms**2)
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
