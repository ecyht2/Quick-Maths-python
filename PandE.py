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
        P = I**2 / R
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

def par_res(resistors, *res):
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
