from math import log10
from statistics import mean,median,mode
from statistics import quantiles, variance, stdev
from math import sqrt, atan, sin, cos

# db Functions
# Power
def db_power(input, output):
    return 10*log10(output/input)
def db_power_reverse(db, initial):
    db = db / 10
    return 10**db * initial
def db_power_find_initial(db, output):
    db = db / 10
    return output / 10**db

# Volts/Current
def db_volts(input, output):
    return 20*log10(output/input)
def db_volts_reverse(db, initial):
    db = number / 20
    return 10**db * initial
def db_volts_find_initial(db, output):
    db = db / 20
    return output / 10**db

# Data Presentation
def stat_range(data):
    minimum = min(data)
    maximum = max(data)
    return maximum - minimum
def IQR(data):
        quartiles = quantile(data)
        return quartiles[2] - quartiles[0]

# Vectors
def coulomb_law(q1, q2, r):
    K = 9e9
    return (K * abs(q1) * abs(q2))/(r**2)
def E_field(q, r):
    K = 9e9
    return (K * abs(q))/(r**2)
