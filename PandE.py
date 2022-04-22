import math
from Constants import *

def par_res(resistors, *res):
    if type(resistors) == list:
        sum(resisstors**-1)**-1

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
