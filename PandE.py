import math
from Constants import *

# Circuit Analysis
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
def ohms_V(I, R):
    return I*R
def ohms_I(V, R):
    return V/R
def ohms_R(V, I):
    return V/I

# Voltage Divider
def voltage_divider(V_total, R_out, R_rest, *res):
    v_out = 0

    v_out = V_total * R_out
    if type(R_rest) == list:
        v_out = v_out/(R_out + sum(R_rest))
    else:
        v_out = v_out/(R_out + R_rest + sum(res))

    return v_out

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
