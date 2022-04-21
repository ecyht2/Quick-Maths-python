from math import log10
from statistics import mean,median,mode
from statistics import quantiles, variance, stdev
from math import sqrt, atan, sin, cos
from Constants import *

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
class myStats:
    mean = 0
    mode = 0
    median = 0
    iqr = 0
    upper_quartile = 0
    lowwer_quartile = 0
    variance = 0
    std = 0
    range = 0
    data = []
    def __init__(self, data):
        self.data = data

        self.mean = mean(data)
        self.mode = mode(data)
        self.median = median(data)

        self.iqr = IQR(data)
        quartiles = quantiles(data)
        self.upper_quartile = quantiles[2]
        self.lower_quartile = quantiles[0]

        self.variance = variance(data)
        self.std = stdev(data)

        self.range = stat_range(data)
    def __get_mean():
        self.mean = mean(data)

class myStatsGrouped(myStats):
    freq_data = {}
    def __init__(self, data):
        self.freq_data = data
        for value in data.keys():
            for frequency in range(data[value]):
                self.data.append(value)

        self.mean = mean(self.data)
        self.mode = mode(self.data)
        self.median = median(self.data)

        self.iqr = IQR(self.data)
        quartiles = quantiles(data)
        self.upper_quartile = quantiles[2]
        self.lower_quartile = quantiles[0]

        self.variance = variance(self.data)
        self.std = stdev(self.data)

        self.range = stat_range(self.data)

def stat_range(data):
    minimum = min(data)
    maximum = max(data)
    return maximum - minimum
def IQR(data):
        quartiles = quantiles(data)
        return quartiles[2] - quartiles[0]

# Vectors
def coulomb_law(q1, q2, r):
    K = 9e9
    return (K * abs(q1) * abs(q2))/(r**2)
def E_field(q, r):
    K = 9e9
    return (K * abs(q))/(r**2)
def mag(a = 0, b = 0, c = 0, l = 0):
    distance = sqrt(a**2 + b**2 + c**2)
    if type(l) == list or type(l) == tuple:
        distance = sqrt(sum(i**2 for i in l))
    return distance

class Vector(list):
    def __init__(self, vector):
        vector = list(vector)
        while len(vector) != 3:
            vector.append(0)
        for i in range(len(vector)):
            self.append(vector[i])
        self.mag = mag(l=vector)

    def cross(self, vector):
        x = self[1]*vector[2] - self[2]*vector[1]
        y = -(self[0]*vector[2] - self[2]*vector[0])
        z = self[0]*vector[1] - self[1]*vector[0]
        return (x, y ,z)

    def dot(self, vector):
        product = vector[0]*self[0] + vector[1]*self[1] + vector[2]*self[2]
        return product

    def add(self, vector):
        return_vector = []
        for i in range(3):
            return_vector.append(self[i] + Vector[i])

        return return_vector

    def subtract(self, vector):
        return_vector = []
        for i in range(3):
            return_vector.append(self[i] - Vector[i])

        return return_vector
