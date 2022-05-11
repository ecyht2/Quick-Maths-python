from math import log10, exp, factorial
from statistics import mean,median,mode
from statistics import quantiles, variance, stdev
from math import sqrt, atan, sin, cos
from Constants import *
import numpy as np

# Helper
def product(iterable, start = 0):
    """

    """
    result = 1
    for i in range(len(iterable) + start):
        result *= iterable[i]

    return result

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
    # Measure of Location
    mean = 0
    mode = 0
    median = 0
    # Measure of Spread
    iqr = 0
    upperQuartile = 0
    lowwerQuartile = 0
    variance = 0
    std = 0
    range = 0
    # Other attributes
    max = 0
    min = 0
    data = []
    def __init__(self, data):
        self.data = data

        self.mean = mean(data)
        self.mode = mode(data)
        self.median = median(data)

        self.iqr = IQR(data)
        quartiles = quantiles(data)
        self.upperQuartile = quartiles[2]
        self.lowerQuartile = quartiles[0]

        self.variance = variance(data)
        self.std = stdev(data)

        self.range = stat_range(data)

        self.max = max(data)
        self.min = min(data)

    def __str__(self):
        info = {
            "mean": self.mean,
            "mode": self.mode,
            "median": self.median,
            "iqr": self.iqr,
            "upper quartile": self.upperQuartile,
            "lower quartile": self.lowerQuartile,
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
            "upper quartile": self.upperQuartile,
            "lower quartile": self.lowerQuartile,
            "variance": self.variance,
            "std": self.std,
            "range": self.range,
        }
        return str(info)

class myStatsGrouped(myStats):
    freqData = {}
    def __init__(self, data: dict):
        self.freqData = data
        self.data = []
        for value in data.keys():
            for frequency in range(data[value]):
                self.data.append(value)

        self.mean = mean(self.data)
        self.mode = mode(self.data)
        self.median = median(self.data)

        self.iqr = IQR(self.data)
        quartiles = quantiles(data)
        self.upperQuartile = quartiles[2]
        self.lowerQuartile = quartiles[0]

        self.variance = variance(self.data)
        self.std = stdev(self.data)

        self.range = stat_range(self.data)

        self.max = max(data)
        self.min = min(data)

class myStatsGroupedRange():
    freqData = {}
    # Measure of Location
    mean = 0
    median = 0
    medianClass = ""
    mode = 0
    # Measure of Spread
    variance = 0
    std = 0
    def __init__(self, data: dict):
        self.freqData = data

        self.boundaries = dict()
        for keys, values in data.items():
            if "-" in keys:
                tempData = keys.split("-")
                self.boundaries[keys] = [int(tempData[0]), int(tempData[1])]

        self.midpoints = dict()
        for boundary in self.boundaries.values():
            self.midpoints[keys] = mean(boundary)

        self.observations = 0
        for freq in data.values():
            self.observations += freq

        # Measure of Location
        self.mean = self.__find_mean()
        self.median = self.__find_median()
        self.mode = max(data, key=data.get)

        # Measure of Spread
        self.variance = self.__find_variance()
        self.std = sqrt(self.variance)

    def __find_mean(self) -> float:
        """
        Finds the mean for the data of grouped data with ranges
        """
        products = []
        for f, m in zip(self.midpoints.values(), self.freqData.values()):
            products.append(f*m)
        mean = sum(products)/self.observations
        return mean

    def __find_median(self) -> float:
        """
        Finds the mean for the data of grouped data with ranges
        """
        # Finding the median class
        sumFreq = 0
        for freqClass, freq in self.freqData.items():
            sumFreq += freq
            if sumFreq >= self.observations/2:
                self.medianClass = freqClass
                break
            previousClass = freqClass

        # Finding the median
        boundary = self.boundaries[self.medianClass]
        medianClass = self.freqData[self.medianClass]
        freqPrevious = self.freqData[previousClass]
        median = boundary[0] + (self.observations/2 - freqPrevious)*(boundary[1] - boundary[0])/medianClass
        return median

    def __find_variance(self) -> float:
        """
        Finds the variance for the data of grouped data with ranges
        """
        variance = 1/(self.observations - 1)
        # Mean Square
        meansquare = (self.mean*self.observations)**2 / self.observations

        # Mean of x square
        products = []
        for m, f in zip(self.midpoints.values(), self.freqData.values()):
            products.append(f*(m**2))
        meanxsquare = sum(products)

        variance *= (meanxsquare - meansquare)
        return variance

    def __str__(self):
        info = {
            "mean": self.mean,
            "mode": self.mode,
            "median": self.median,
            "median class": self.medianClass,
            "variance": self.variance,
            "std": self.std,
        }
        return str(info)

    def __repr__(self):
        info = {
            "mean": self.mean,
            "mode": self.mode,
            "median": self.median,
            "median class": self.medianClass,
            "variance": self.variance,
            "std": self.std,
        }
        return str(info)

def stat_range(data):
    minimum = min(data)
    maximum = max(data)
    return maximum - minimum
def IQR(data):
        quartiles = quantiles(data)
        return quartiles[2] - quartiles[0]

# Vectors
# Classes
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
        if not (type(vector) == list or type(vector) == tuple or type(vector) == type(self)):
            raise TypeError("A vector must be an array_like (list, tuple, Vector) object")
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
        """
        Returns the cross product of the vectors (self)x(vector)

        Parameters
        ----------
        vector
            The vector to perform the cross product with

        Returns
        -------
        Vector
            The cross product
        """
        # Checking vector type
        if not type(vector) == type(self):
            raise TypeError("Only a type Vector can be crossed with a Vector")

        # a x b = absin(thetha)c
        x = self[1]*vector[2] - self[2]*vector[1]
        y = -(self[0]*vector[2] - self[2]*vector[0])
        z = self[0]*vector[1] - self[1]*vector[0]
        return_vector = [x, y, z]
        return Vector(return_vector)

    def dot(self, vector):
        """
        Returns the dot product of the vectors (self).(vector)

        Parameters
        ----------
        vector
            The vector to perform the dot product with

        Returns
        -------
        float
            The dot product
        """
        # Checking vector type
        if not type(vector) == type(self):
            raise TypeError("Only a type Vector can be dot with a Vector")

        # a.b = |a||b|cos(thetha)
        product = vector[0]*self[0] + vector[1]*self[1] + vector[2]*self[2]
        return product

    def __add__(self, vector):
        """
        Returns self + vector
        """
        # Checking vector type
        if not type(vector) == type(self):
            raise TypeError("Only a type Vector can be added with a Vector")

        vector = Vector(vector)
        return_vector = []
        for i in range(3):
            return_vector.append(self[i] + vector[i])

        return Vector(return_vector)

    def __sub__(self, vector):
        """
        Returns self - vector
        """
        # Checking vector type
        if not type(vector) == type(self):
            raise TypeError("Only a type Vector can be subtracted with a Vector")

        vector = Vector(vector)
        return_vector = []
        for i in range(3):
            return_vector.append(self[i] - vector[i])

        return Vector(return_vector)

    def __mul__(self, scalar):
        """
        Returns self * scalar
        """
        # Checking vector type
        if not (type(scalar) == int or type(scalar) == float):
            raise TypeError("Invalid scalar value for scalar multiplication")

        return_vector = []
        for i in range(3):
            return_vector.append(self[i] * scalar)

        return Vector(return_vector)

    def unit_vector(self):
        """
        Returns the unit vector of self

        Returns
        -------
        Vector
            The unit vector
        """
        unit_vector = []

        # Finding values of x, y  and z
        # (x̂ + ŷ + ẑ)/(sqrt(x̂ + ŷ + ẑ))
        for i in self:
            unit_vector = self * (1 / self.mag)

        return Vector(unit_vector)

class Charge(Vector):
    """
    Create a charge object

    Parameters
    ----------
    q
        The charge of the charge in coulomb (C)
    vector
        The vector location of the charge
    """
    q = 0
    def __init__(self, q, vector):
        self.q = q
        # Checking if parameter is valid
        # Type
        if not (type(vector) == list or type(vector) == tuple or type(vector) == type(self)):
            raise TypeError("A vector must be an array_like (list, tuple, Vector) object")
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

    def E_Field(self, location: Vector):
        """
        Calculates the Electric Field produced by the charge at a point

        Parameters
        ----------
        location
            The location vector

        Returns
        -------
        Vector
            The Electric Field Vector
        """
        # Converting to required format
        location = Vector(location)
        self_vector = list(self)
        self_vector = Vector(self_vector)

        # Determining force direction
        if self.q > 0:
            diff = location - self_vector
        else:
            diff = self_vector - location
        r = diff.mag

        # Returning result
        return E_field(self.q, r = r, vector=diff)

    def Electric_Force(self, charge):
        """
        Calculates the Electric Force Caused by this charge to another

        Parameters
        ----------
        charge
            The charge which the force is act apon
        """
        # Checking if parameters is a compatible type
        if not type(charge) == type(self):
            raise TypeError("Charge given must be of type Charge")

        # Converting to required format
        charge_vector = Vector(list(charge))
        self_vector = Vector(list(self))

        # Determine force direction
        if (self.q > 0 and charge.q > 0) or (self.q < 0 and charge.q < 0):
            diff = charge_vector - self_vector
        else:
            diff = self_vector - charge_vector
        r = diff.mag

        # Returning result
        return coulomb_law(self.q, charge.q, r=r, vector=diff)

# Electric Field
def coulomb_law(q1, q2, r, vector = None, k = kCoulomb):
    magnitude = (k * abs(q1) * abs(q2))/(r**2)
    if not vector == None:
        vector = Vector(vector)
        force = vector.unit_vector() * magnitude
        return force
    else:
        return magnitude
def electric_force(q: float, E: float) -> float:
    F = 0
    F = abs(q)*abs(E)
    return F
def E_field(q = 0, r = 0, sigma = 0, epsilon0 = epsilon0, k = kCoulomb, vector = None):
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

# Magnetic Field
def B_field(H: float = 0, I: float = 0, r: float = 0, mu0: float = mu0, mur: float = 1) -> float:
    """
    Calculates the magnetic field
    """
    B = 0
    if I == 0 and H == 0:
        raise ValueError("No I or H was given")
    elif H == 0:
        B = B_field_I(r, I, mu0, mur)
    elif I == 0:
        B = B_field_H(r, H, mu0, mur)

    return B
def B_field_H(H: float, mu0: float = mu0, mur: float = 1) -> float:
    """
    Calculates the magnetic field given H
    """
    B = mu0*mur*H
    return B
def B_field_I(r: float, I: float, mu0: float = mu0, mur: float = 1) -> float:
    """
    Calculates the magnetic field given the current
    """
    H = I/(2*pi*r)
    B = mu0*mur*H
    return B
def magnetic_force(q, v, B):
    """
    Calculate the magnetic force on charge q caused by B
    """
    F = abs(q)*v*abs(B)
    return F
def EMF(v = 0, B = 0, L = 0, E = 0):
    if E == 0:
        return v*B*L
    else:
        return E*L

# Electromagnetic Waves
def energy_density_E(E: float, epsilon0: float = epsilon0) -> float:
    """
    Calculate the energy density for the electric field of the EM wave
    """
    uE = 0.5 * epsilon0*E**2
    return uE
def energy_density_B(B: float, mu0: float = mu0) -> float:
    """
    Calculate the energy density for the magnetic field of the EM wave
    """
    uB = 0.5 * B**2/mu0
    return uB
def energy_density(B: float = 0, E: float = 0,
                   mu0: float = mu0, epsilon0: float = epsilon0,
                   uE: float = 0, uB: float = 0) -> float:
    """
    Calculate the energy density of the EM wave
    """
    u = 0

    if E > 0:
        uE = energy_density_E(E, epsilon0)
    if B > 0:
        uB = energy_density_B(B, mu0)

    if uE == 0 and uB > 0:
        u = 2*uB
    elif uB == 0 and uE > 0:
        u = 2*uE
    elif uB > 0 and uE > 0:
        u = uB + uE
    else:
        raise ValueError("Atleast one of these: B, E, uE or uB must be passed in")

    return u
def EM_E_field(B: float, c: float = c) -> float:
    """
    Calculates the Electric Field of an EM wave
    """
    E = c*B
    return E
def EM_B_field(E: float, c: float = c,
               epsilon0: float = epsilon0, mu0: float = mu0) -> float:
    """
    Calculates the Magnetic Field of an EM wave
    """
    B = epsilon0 * mu0 * c * E
    # B = sqrt(mu0*epsilon0) * E
    return B
def poynting_vector(E: Vector, B: Vector, mu0: float = mu0) -> float:
    """
    Calculates the poynting vector of and EM wave
    """
    S = 1/mu0 * E.cross(B)
    return S

# Vector Maths
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
        if not (type(vector) == list or type(vector) == tuple or type(vector) == Vector):
            raise TypeError("Invalid vector passed in. It must be an array_like (list, tuple, Vector) object")
        if len(vector) > 3:
            raise ValueError("Invalid vector, a vector must include at most 3 values")
        distance = sqrt(sum(i**2 for i in vector))
    else:
        distance = sqrt(a**2 + b**2 + c**2)
    return distance
