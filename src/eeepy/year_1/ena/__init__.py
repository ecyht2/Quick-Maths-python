#!/usr/bin/env python3
"""Functions and Equations taught in EEEE1030 Engineering Mathematics
module."""
from math import pi, sqrt
from statistics import mean, median, mode, quantiles, stdev, variance
from typing import Iterable, Union

from eeepy.utils.constants import c, epsilon0, kCoulomb, mu0


# Data Presentation
class MyStats:
    """A class that calculates a statistics of a given data.

    :param data: The data.
    """
    def __init__(self, data: Iterable):
        """A class that calculates a statistics of a given data.

        :param data: The data.
        """
        self.data = tuple(data)

        quartiles = quantiles(data)
        self._upper_quartile = quartiles[2]
        self._lower_quartile = quartiles[0]

    @property
    def mean(self) -> float:
        """Return the mean of the data."""
        return mean(self.data)

    @property
    def median(self) -> float:
        """Return the median of the data."""
        return median(self.data)

    @property
    def mode(self) -> float:
        """Return the mode of the data."""
        return mode(self.data)

    @property
    def iqr(self) -> float:
        """Return the IRO of the data."""
        return IQR(self.data)

    @property
    def upper_quartile(self) -> float:
        """Return the upper quartile of the data."""
        return self._upper_quartile

    @property
    def lower_quartile(self) -> float:
        """Return the lower quartile of the data."""
        return self._lower_quartile

    @property
    def var(self) -> float:
        """Return the variance of the data."""
        return variance(self.data)

    @property
    def std(self) -> float:
        """Return the standard deviation of the data."""
        return stdev(self.data)

    @property
    def stat_range(self) -> float:
        """Return the range of the data."""
        return stat_range(self.data)

    @property
    def stat_min(self) -> float:
        """Return the minimum of the data."""
        return min(self.data)

    @property
    def stat_max(self) -> float:
        """Return the maximum of the data."""
        return max(self.data)

    def __str__(self):
        """Returns str(self)."""
        info = {
            "mean": self.mean,
            "mode": self.mode,
            "median": self.median,
            "iqr": self.iqr,
            "upper quartile": self._upper_quartile,
            "lower quartile": self._lower_quartile,
            "variance": self.var,
            "std": self.std,
            "range": self.stat_range,
        }
        return str(info)

    def __repr__(self):
        """Returns repr(self)."""
        info = {
            "mean": self.mean,
            "mode": self.mode,
            "median": self.median,
            "iqr": self.iqr,
            "upper quartile": self.upper_quartile,
            "lower quartile": self.lower_quartile,
            "variance": self.var,
            "std": self.std,
            "range": self.stat_range,
        }
        return str(info)


class MyStatsGrouped(MyStats):
    """A class that calculates a statistics of a given data with
    frequencies.

    :param data: The data.
    """
    def __init__(self, data: dict):
        """A class that calculates a statistics of a given data with
        frequencies.

        :param data: The data.
        """
        self.freq_data = data
        self.data = []
        for value in data.keys():
            for _frequency in range(data[value]):
                self.data.append(value)

        super().__init__(self.data)


class MyStatsGroupedRange():
    """A class that calculates a statistics of a given a grouped data.

    :param data: The data.
    """
    def __init__(self, data: dict):
        """A class that calculates a statistics of a given a grouped data.

        :param data: The data.
        """
        self.freq_data = data

        self.boundaries = {}
        for keys, _values in data.items():
            if "-" in keys:
                tempData = keys.split("-")
                self.boundaries[keys] = [int(tempData[0]), int(tempData[1])]

        self.midpoints = {}
        for keys, boundary in self.boundaries.items():
            self.midpoints[keys] = mean(boundary)

        self.observations = 0
        for freq in data.values():
            self.observations += freq

        # Finding the median class
        sum_freq = 0
        for freq_class, freq in self.freq_data.items():
            sum_freq += freq
            if sum_freq >= self.observations / 2:
                self.median_class = freq_class
                break
            self.previous_class = freq_class

    @property
    def mean(self) -> float:
        """Return the mean of the data."""
        return self.__find_mean()

    @property
    def median(self) -> float:
        """Return the median of the data."""
        return self.__find_median()

    @property
    def mode(self) -> float:
        """Return the mode of the data."""
        return max(self.freq_data, key=self.freq_data.get)

    @property
    def var(self) -> float:
        """Return the variance of the data."""
        return self.__find_variance()

    @property
    def std(self) -> float:
        """Return the standard deviation of the data."""
        return sqrt(self.var)

    def __find_mean(self) -> float:
        """Finds the mean for the data of grouped data with ranges.

        :returns: The mean.
        """
        products = []
        for f, m in zip(self.midpoints.values(), self.freq_data.values()):
            products.append(f * m)
        data_mean = sum(products) / self.observations
        return data_mean

    def __find_median(self) -> float:
        """Finds the mean for the data of grouped data with ranges.

        :returns: The median.
        """
        # Finding the median
        boundary = self.boundaries[self.median_class]
        self.median_class = self.freq_data[self.median_class]
        freq_previous = self.freq_data[self.previous_class]
        med = boundary[0] + (self.observations / 2 - freq_previous)\
            * (boundary[1] - boundary[0]) / self.median_class
        return med

    def __find_variance(self) -> float:
        """Finds the variance for the data of grouped data with ranges.

        :returns: The variance.
        """
        var = 1 / (self.observations - 1)
        # Mean Square
        meansquare = (self.mean * self.observations)**2 / self.observations

        # Mean of x square
        products = []
        for m, f in zip(self.midpoints.values(), self.freq_data.values()):
            products.append(f * (m**2))
        meanxsquare = sum(products)

        var *= (meanxsquare - meansquare)
        return var

    def __str__(self):
        """Returns str(self)."""
        info = {
            "mean": self.mean,
            "mode": self.mode,
            "median": self.median,
            "median class": self.median_class,
            "variance": self.var,
            "std": self.std,
        }
        return str(info)

    def __repr__(self):
        """Returns repr(self)."""
        info = {
            "mean": self.mean,
            "mode": self.mode,
            "median": self.median,
            "median class": self.median_class,
            "variance": self.var,
            "std": self.std,
        }
        return str(info)


def stat_range(data: Iterable) -> float:
    """The range of a given data.

    :return: The range.
    """
    minimum = min(data)
    maximum = max(data)
    return maximum - minimum


def IQR(data: Iterable) -> float:
    """The interquartile range of a given data.

    :return: The interquartile range (IQR).
    """
    quartiles = quantiles(data)
    return quartiles[2] - quartiles[0]


# Vectors
# Classes
class Vector():
    """Create a Vector.

    Parameters
    ----------
    vector
        A list or tuple with the values of (x, y, z),
        0 will be added for the z value if only x and y are given
    """
    mag = 0

    def __init__(self, vector: Union[list, tuple]):
        # Checking if parameter is valid
        # Type
        if not isinstance(type(vector), Union[list, tuple, Vector]):
            raise TypeError("A vector must be an array_like"
                            "(list, tuple, Vector) object")
        # Size
        if len(vector) > 3:
            raise ValueError("Invalid vector to be converted")
        # Indexes
        for i in vector:
            if not issubclass(type(i), Union[float, int]):
                raise ValueError("Invalid vector to be converted")

        vector = list(vector)
        while len(vector) != 3:
            vector.append(0)
        self.data = vector
        self.mag = mag(vector=vector)

    def cross(self, vector):
        """Returns the cross product of the vectors (self)x(vector).

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
        if not issubclass(type(vector), Union[list, tuple, Vector]):
            raise TypeError("Only a type Vector can be crossed with a Vector")
        if len(vector) > 3:
            raise ValueError("Invalid vector")

        # a x b = absin(thetha)c
        x = self[1] * vector[2] - self[2] * vector[1]
        y = -(self[0] * vector[2] - self[2] * vector[0])
        z = self[0] * vector[1] - self[1] * vector[0]
        return_vector = [x, y, z]
        return Vector(return_vector)

    def dot(self, vector):
        """Returns the dot product of the vectors (self).(vector).

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
        if not issubclass(type(vector), Union[list, tuple, Vector]):
            raise TypeError("Only a type Vector can be dot with a Vector")
        if len(vector) > 3:
            raise ValueError("Invalid vector")

        # a.b = |a||b|cos(thetha)
        product = vector[0] * self[0]\
            + vector[1] * self[1]\
            + vector[2] * self[2]
        return product

    def __getitem__(self, key: int) -> float:
        """x.__getitem__(y) <==> x[y]."""
        return self.data[key]

    def __repr__(self):
        """Return repr(self)."""
        return str(self.data)

    def __str__(self):
        """Return str(self)."""
        return str(self.data)

    def __len__(self):
        """Return len(self)."""
        return len(self.data)

    def __add__(self, vector):
        """Returns self + vector."""
        # Checking vector type
        if not issubclass(type(vector), Union[list, tuple, Vector]):
            return NotImplemented
        if len(vector) > 3:
            raise ValueError("Invalid vector")

        vector = Vector(vector)
        return_vector = []
        for i in range(3):
            return_vector.append(self[i] + vector[i])

        return Vector(return_vector)

    def __radd__(self, vector):
        """Returns vector + self."""
        return_vector = self + vector

        return Vector(return_vector)

    def __sub__(self, vector):
        """Returns self - vector."""
        # Checking vector type
        if not issubclass(type(vector), Union[list, tuple, Vector]):
            return NotImplemented
        if len(vector) > 3:
            raise ValueError("Invalid vector")

        vector = Vector(vector)
        return_vector = []
        for i in range(3):
            return_vector.append(self[i] - vector[i])

        return Vector(return_vector)

    def __rsub__(self, vector):
        """Returns vector - self."""
        if not issubclass(type(vector), Union[list, tuple, Vector]):
            return NotImplemented
        if len(vector) > 3:
            raise ValueError("Invalid vector")

        vector = Vector(vector)
        return_vector = []
        for i in range(3):
            return_vector.append(vector[i] - self[i])

        return Vector(return_vector)

    def __mul__(self, scalar: float):
        """Performs Scalar multiplication of Vector."""
        # Checking vector type
        if not isinstance(scalar, Union[int, float]):
            return NotImplemented

        return_vector = []
        for i in range(3):
            return_vector.append(self[i] * scalar)

        return Vector(return_vector)

    def __rmul__(self, scalar: float):
        """Performs Scalar multiplication of Vector."""
        # Checking vector type
        if not isinstance(scalar, Union[int, float]):
            return NotImplemented

        return_vector = []
        for i in range(3):
            return_vector.append(self[i] * scalar)

        return Vector(return_vector)

    def __truediv__(self, scalar: float):
        """Performs Scalar division of Vector."""
        # Checking vector type
        if not isinstance(scalar, Union[int, float]):
            return NotImplemented

        return_vector = []
        for i in range(3):
            return_vector.append(self[i] / scalar)

        return Vector(return_vector)

    def __abs__(self):
        """Returns the Magnitude of the Vector."""
        return self.mag

    def unit_vector(self):
        """Returns the unit vector of self.

        Returns
        -------
        Vector
            The unit vector
        """
        unit_vector = []

        # Finding values of x, y  and z
        # (x̂ + ŷ + ẑ)/(sqrt(x̂ + ŷ + ẑ))
        for _ in self:
            unit_vector = self * (1 / self.mag)

        return Vector(unit_vector)


class Charge(Vector):
    """Create a Charge object.

    Parameters
    ----------
    q
        The charge of the charge in coulomb (C)
    vector
        The vector location of the charge
    """
    q = 0

    def __init__(self, q: float, vector: Union[list, tuple]):
        self.q = q
        super().__init__(vector)

    def e_field(self, location: Vector) -> Vector:
        """Calculates the Electric Field produced by the charge at a point.

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
        return e_field(self.q, r=r, vector=diff)

    def electric_force(self, charge: float) -> Vector:
        """Calculates the Electric Force Caused by this charge to another.

        Parameters
        ----------
        charge
            The charge which the force is act apon
        """
        # Checking if parameters is a compatible type
        if not isinstance(charge, type(self)):
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
def coulomb_law(q1: float, q2: float, r: float,
                vector: Vector = None, k: float = kCoulomb) -> float or Vector:
    """Calculates the magnitude force induced on a charge.
    If a vector is given, the vector of the force is returned instead.

    Parameters
    ----------
    q1
        The charge of on of the charge
    q2
        The charge of on of the charge
    r
        The distance between the charges
    vector (optional)
        The vector of the charge in which the force is acted on
    k (optional)
        Coulomb's Constant

    Returns
    -------
    float if no vector is given
        The magnitude of the force
    Vector
        The force on the charge
    """
    magnitude = (k * abs(q1) * abs(q2)) / (r**2)
    if vector is not None:
        vector = Vector(vector)
        force = vector.unit_vector() * magnitude
        return force

    return magnitude


def electric_force(q: float, E: float) -> float:
    """Calculates the electric force caused by an electric field on a charge.

    Parameters
    ----------
    q
        The charge in which the force is acted on
    E
        The electric field

    Returns
    -------
    float
        The electric force
    """
    F = 0
    F = abs(q) * abs(E)
    return F


def e_field(q: float = 0, r: float = 0,
            sigma: float = 0, k: float = kCoulomb,
            vector: float = None) -> float or Vector:
    """Calculates the electric field of a charge.

    Parameters
    ----------
    q (Not needed if sigma and epsilon0 is provided)
        The charge
    r (Not needed if sigma and epsilon0 is provided)
        The distance from the charge
    sigma (Not needed if q and r is provided)
        The sigma value
    k (optional)
        Coulomb's constant
    vector (optional)
        The vector of the charge

    Returns
    -------
    float (if no vector is given)
        The magnitude of the electric field
    Vector
        The Electric Field of the charge
    """
    eField = 0
    if sigma == 0:
        eField = (k * abs(q)) / (r**2)
    else:
        eField = sigma / (2 * epsilon0)

    if vector is not None:
        vector = Vector(vector)
        return_value = vector.unit_vector() * eField
    else:
        return_value = eField

    return return_value


# Magnetic Field
def B_field(H: float = 0, current: float = 0, r: float = 0,
            mur: float = 1) -> float:
    """Calculates the magnetic field.

    Parameters
    ----------
    H (Not needed if I and r is given)
        Magnetic field intensity
    current (Not needed if H is given)
        The current of the wire
    r (Not needed if H is given)
        The distance from the wire
    mur
        Relative permeability of the material

    Returns
    -------
    float
        The magnetic field
    """
    B = 0
    if current == 0 and H == 0:
        raise ValueError("No I or H was given")
    if H == 0:
        B = B_field_I(r, current, mur)
    elif current == 0:
        B = B_field_H(H, mur)

    return B


def B_field_H(H: float, mur: float = 1) -> float:
    """Calculates the magnetic field given H.

    Parameters
    ----------
    H (Not needed if I and r is given)
        Magnetic field intensity
    mur
        Relative permeability of the material

    Returns
    -------
    float
        The magnetic field
    """
    B = mu0 * mur * H
    return B


def B_field_I(r: float, current: float, mur: float = 1) -> float:
    """Calculates the magnetic field given the current.

    Parameters
    ----------
    r
        The distance from the wire
    current
        The current of the wire
    mur
        Relative permeability of the material

    Returns
    -------
    float
        The magnetic field
    """
    H = current / (2 * pi * r)
    B = mu0 * mur * H
    return B


def magnetic_force(q: float, v: float, B: float) -> float:
    """Calculate the magnetic force on charge q caused by B.

    Parameters
    ----------
    q
        The charge
    v
        The speed the charge is moving at
    B
        The magnetic field the the charge is moving in

    Returns
    -------
    float
        The magnetic force
    """
    F = abs(q) * v * abs(B)
    return F


def EMF(v: float = 0, B: float = 0, L: float = 0, E: float = 0) -> float:
    """Calculates the EMF induced by a moving rod in a magnetic field.

    Parameters
    ----------
    v
        The speed the rod is moving at
    B
        The magnetic field
    L
        The length of the rod
    E
        The electric field

    Returns
    -------
    float
        The EMF induced
    """
    if E == 0:
        return v * B * L
    return E * L


# Electromagnetic Waves
def energy_density_E(E: float) -> float:
    """Calculate the energy density for the electric field of the EM wave.

    Parameters
    ----------
    E
        The electric field

    Returns
    -------
    float
        The energy density
    """
    uE = 0.5 * epsilon0 * E**2
    return uE


def energy_density_B(B: float) -> float:
    """Calculate the energy density for the magnetic field of the EM wave.

    Parameters
    ----------
    B
        The magnetic field

    Returns
    -------
    float
        The energy density
    """
    uB = 0.5 * B**2 / mu0
    return uB


def energy_density(B: float = 0, E: float = 0,
                   uE: float = 0, uB: float = 0) -> float:
    """Calculate the energy density of the EM wave.

    Parameters
    ----------
    B (Not needed if uE or uB or E is given)
        The magnetic field
    E (Not needed if uE or uB or B is given)
        The electric field
    uE (Not needed if B or E or uB is given)
        Electric field energy density
    uB (Not needed if B or E or uE is given)
        Magnetic field energy density

    Returns
    -------
    float
        The energy density
    """
    u = 0

    if E > 0:
        uE = energy_density_E(E)
    if B > 0:
        uB = energy_density_B(B)

    if uE == 0 and uB > 0:
        u = 2 * uB
    elif uB == 0 and uE > 0:
        u = 2 * uE
    elif uB > 0 and uE > 0:
        u = uB + uE
    else:
        raise ValueError("Atleast one of these: B, E, uE or uB"
                         "must be passed in")

    return u


def EM_E_field(B: float) -> float:
    """Calculates the Electric Field of an EM wave.

    Parameters
    ----------
    B
        Magnetic field of the EM wave

    Returns
    -------
    float
        The electric field
    """
    E = c * B
    return E


def EM_B_field(E: float) -> float:
    """Calculates the Magnetic Field of an EM wave.

    Parameters
    ----------
    E
        Electric field of the EM wave

    Returns
    -------
    float
        The magnetic field
    """
    B = epsilon0 * mu0 * c * E
    # B = sqrt(mu0*epsilon0) * E
    return B


def poynting_vector(E: Vector, B: Vector) -> float:
    """Calculates the poynting vector of and EM wave.

    Parameters
    ----------
    E
        Electric field of the EM wave
    B
        Magnetic field of the EM wave

    Returns
    -------
    float
        The poynting vector of the EM wave
    """
    S = 1 / mu0 * E.cross(B)
    return S


# Vector Maths
def mag(a: float = 0, b: float = 0, c_vec: float = 0,
        vector: Vector = None) -> float:
    """Finds the magnitude of any vector of (a, b, c) or vector.

    Parameters
    ----------
    a
        x̂ value of the vector
    b
        ŷ value of the vector
    c_vec
        ẑ value of the vector
    vector (Not needed if a, b and c is given)
        A vector of type tuple, list or Vector

    Returns
    -------
    float
        The magnitude of the Vector
    """
    if vector is not None:
        if not issubclass(type(vector), Union[list, tuple, Vector]):
            raise TypeError("Invalid vector passed in. It must be an"
                            "array_like (list, tuple, Vector) object")
        if len(vector) > 3:
            raise ValueError("Invalid vector, a vector must"
                             " include at most 3 values")
        distance = sqrt(sum(i**2 for i in vector))
    else:
        distance = sqrt(a**2 + b**2 + c_vec**2)
    return distance
