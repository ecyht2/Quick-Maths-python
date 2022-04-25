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
        self.upper_quartile = quartiles[2]
        self.lower_quartile = quartiles[0]

        self.variance = variance(data)
        self.std = stdev(data)

        self.range = stat_range(data)

    def __str__(self):
        info = {
            "mean": self.mean,
            "mode": self.mode,
            "median": self.median,
            "iqr": self.iqr,
            "upper quartile": self.upper_quartile,
            "lower quartile": self.lower_quartile,
            "variance": self.variance,
            "std": self.std,
            "range": self.range,
        }
        return str(info)

    def __repr__(self):
        info = {
            "mean": self.mean,
            "mode": self.mode,
            "median": self.median,
            "iqr": self.iqr,
            "upper quartile": self.upper_quartile,
            "lower quartile": self.lower_quartile,
            "variance": self.variance,
            "std": self.std,
            "range": self.range,
        }
        return str(info)

class myStatsGrouped(myStats):
    freq_data = {}
    def __init__(self, data):
        self.freq_data = data
        self.data = []
        for value in data.keys():
            for frequency in range(data[value]):
                self.data.append(value)

        self.mean = mean(self.data)
        self.mode = mode(self.data)
        self.median = median(self.data)

        self.iqr = IQR(self.data)
        quartiles = quantiles(data)
        self.upper_quartile = quartiles[2]
        self.lower_quartile = quartiles[0]

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
def coulomb_law(q1, q2, r, vector = None, k = k):
    magnitude = (k * abs(q1) * abs(q2))/(r**2)
    if not vector == None:
        vector = Vector(vector)
        force = vector.unit_vector() * magnitude
        return force
    else:
        return magnitude
def E_field(q = 0, r = 0, sigma = 0, epsilon0 = epsilon0, k = k, vector = None):
    eField = 0
    if sigma == 0:
        eField = (k * abs(q))/(r**2)
    else:
        eField = sigma / (2*epsilon0)

    if not vector == None:
        vector = Vector(vector)
        return_value = vector.unit_vector() * eField
    else:
        return_value = eField

    return return_value
def EMF(v = 0, B = 0, L = 0, E = 0):
    if E == 0:
        return v*B*L
    else:
        return E*L
def mag(a = 0, b = 0, c = 0, vector = None) -> float:
    """
    Finds the magnitude of any vector of (a, b, c) or vector

    Parameters
    ----------
    a
        x̂ value of the vector
    b
        ŷ value of the vector
    c
        ẑ value of the vector
    vector
        A vector of type tuple, list or Vector

    Returns
    -------
    float
        The magnitude of the Vector
    """
    if not vector == None:
        if type(vector) == list or type(vector) == tuple or type(vector) == Vector:
            raise TypeError("Invalid vector passed in. It must be an array_like object")
        if len(vector) > 3:
            raise ValueError("Invalid vector, a vector must include at most 3 values")
        distance = sqrt(sum(i**2 for i in vector))
    else:
        distance = sqrt(a**2 + b**2 + c**2)
    return distance

class Vector(list):
    """
    Create a Vector

    Parameters
    ----------
    vector
        A list or tuple with the values of (x, y, z), 0 will be added for the z value if only x and y are given
    """
    mag = 0
    def __init__(self, vector):
        # Checking if parameter is valid
        # Type
        if not (type(vector) == list or type(vector) == tuple or type(vector) == Vector):
            raise TypeError("Vector() argument must be an array_like (list, tuple, Vector) object")
        # Size
        if len(vector) > 3:
            raise ValueError("Invalid vector to be converted")
        # Indexes
        else:
            for i in vector:
                if not (type(i) == int or type(i) == float):
                    raise ValueError("Invalid vector to be converted")

        vector = list(vector)
        while len(vector) != 3:
            vector.append(0)
        for i in range(len(vector)):
            self.append(vector[i])
        self.mag = mag(vector=vector)

    def cross(self, vector):
        # Checking vector type
        if not type(vector) == Vector:
            raise TypeError("Only a type Vector can be crossed with a Vector")

        # a x b = absin(thetha)c
        x = self[1]*vector[2] - self[2]*vector[1]
        y = -(self[0]*vector[2] - self[2]*vector[0])
        z = self[0]*vector[1] - self[1]*vector[0]
        return Vector([x, y ,z])

    def dot(self, vector):
        # Checking vector type
        if not type(vector) == Vector:
            raise TypeError("Only a type Vector can be dot with a Vector")

        # a.b = |a||b|cos(thetha)
        product = vector[0]*self[0] + vector[1]*self[1] + vector[2]*self[2]
        return product

    def __add__(self, vector):
        # Checking vector type
        if not type(vector) == Vector:
            raise TypeError("Only a type Vector can be added with a Vector")

        vector = Vector(vector)
        return_vector = []
        for i in range(3):
            return_vector.append(self[i] + vector[i])

        return Vector(return_vector)

    def __sub__(self, vector):
        # Checking vector type
        if not type(vector) == Vector:
            raise TypeError("Only a type Vector can be subtracted with a Vector")

        vector = Vector(vector)
        return_vector = []
        for i in range(3):
            return_vector.append(self[i] - vector[i])

        return Vector(return_vector)

    def __mul__(self, scalar):
        # Checking vector type
        if not (type(scalar) == int or type(scalar) == float):
            raise TypeError("Invalid scalar value for scalar multiplication")

        return_vector = []
        for i in range(3):
            return_vector.append(self[i] * scalar)

        return Vector(return_vector)

    def unit_vector(self):
        unit_vector = []

        # Finding values of x, y  and z
        # (x̂ + ŷ + ẑ)/(sqrt(x̂ + ŷ + ẑ))
        for i in self:
            unit_vector = self * (1 / self.mag)

        return Vector(unit_vector)
