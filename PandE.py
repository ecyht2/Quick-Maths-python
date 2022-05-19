import math
from Constants import *
from Maths import product

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

# RMS
def RMS_sinusoidal(value: float) -> float:
    """
    Converts a sinusoidal sigal peak value into RMS
    """
    return value / 2**0.5
def RMS_sinusoidal_reverse(value: float) -> float:
    """
    Converts an RMS value into sinusoidal sigal peak value
    """
    return value * 2**0.5

# Impedence
def imp_cap(capacitor, freq):
    return 1/(2*pi*freq*capacitor*1j)
def imp_ind(ind, freq):
    return 2*pi*freq*ind*1j
def imp_cap_rev(imp, freq):
    return 1/(2*pi*freq*capacitor*1j)
def imp_ind(ind, freq):
    return 2*pi*freq*ind*1j

# Horse Power
def hp_to_W(hp):
    W = 746*hp
    return W
def W_to_hp(W):
    hp = W/746
    return hp

#RPM
def rpm_to_rad(rpm):
    rad = rpm * 2*pi/60
    return rad
def rad_to_rpm(rad):
    rpm = rad * 60/(2*pi)
    return rpm
