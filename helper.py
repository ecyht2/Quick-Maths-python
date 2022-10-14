#!/usr/bin/env python3
from math import exp, factorial, sqrt
from math import log10, log1p, log2
from statistics import mean, median, mode
from statistics import quantiles, variance, stdev
from math import acos, atan, asin
from math import tan, sin, cos
from math import radians, degrees
import cmath

# Imports
# Math Functions
exp = exp,
factorial = factorial,
sqrt = sqrt
log10 = log10
log1p = log1p
log2 = log2
# Stats Function
# Distance
mean = mean
median = median
mode = mode
# Spread
quantiles = quantiles
variance = variance
stdev = stdev
# Trigonometric
# Reverse Functions
acos = acos
atan = atan
asin = asin
# Normal Functions
tan = tan
sin = sin
cos = cos


class Complex():
    """Create a complex number from a real part and an imaginary part."""
    real: float = 0
    imag: float = 0
    r: float = 0
    phi: float = 0
    cartesian: complex = 0 + 0j
    polar: tuple[float, float] = (0, 0)

    def __init__(self, Z: complex):
        # Settings
        if not hasattr(self, "isRadians"):
            self.isRadians = True
        else:
            self.isRadians = self.isRadians
        if not hasattr(self, "isPolar"):
            self.isPolar = False
        else:
            self.isPolar = self.isPolar
        # Forms
        self.cartesian = complex(Z)
        self.polar = cmath.polar(Z)
        # Attributes
        # Cartesian Attributes
        self.real = Z.real
        self.imag = Z.imag
        # Polar Attributes
        self.r = abs(Z)
        if self.isRadians:
            self.phi = cmath.phase(Z)
        else:
            self.phi = degrees(cmath.phase(Z))
            self.polar = (self.polar[0], self.polar[1])

    @classmethod
    def from_polar(cls, r: float, phi: float, radian: bool = True):
        if not radians:
            phi = radians(phi)
        cls.isRadians = radian
        cls.isPolar = True
        Z = cmath.rect(r, phi)
        return cls(Z)

    @classmethod
    def from_polar_tuple(cls, Z: tuple[float, float], radian: bool = True):
        if not radians:
            Z = cmath.rect(Z[0], radians(Z[1]))
        else:
            Z = cmath.rect(Z[0], Z[1])
        cls.isRadians = radian
        cls.isPolar = True
        return cls(Z)

    @classmethod
    def from_cartesian(cls, Z: complex, radian: bool = True):
        cls.isRadians = radian
        cls.isPolar = False
        return cls(Z)

    def __get_class(self, value, polar, radian):
        """Returns a complex class of value with in the
        form specified by polar and radian."""
        # Cartesian
        if not polar:
            ret = Complex.from_cartesian(value, radian)
        # Polar
        else:
            polarForm = cmath.polar(value)
            # Converting to degrees
            if self.isRadians:
                phi = polarForm[1]
            else:
                phi = degrees(polarForm[1])
            ret = Complex.from_polar(polarForm[0], phi, radian)
        return ret

    def conjugate(self, polar: bool = None, radian: bool = None):
        """Return the complex conjugate of its argument.

        Parameters
        ----------
        polar
        radian

        Returns
        -------
        """
        # Manual overide polar
        if polar is not None:
            polar = polar
        else:
            polar = self.isPolar
        # Manual overide radians
        if radian is not None:
            radian = radian
        else:
            radian = self.isRadians

        conjugate = self.cartesian.conjugate()
        # Getting return class
        ret = self.__get_class(conjugate, polar, radian)
        return ret

    def __add__(self, value: complex):
        """Returns self + value."""
        # Chacking validity
        if not hasattr(value, "imag") or not hasattr(value, "real")\
                or not hasattr(value, "conjugate"):
            return NotImplemented
        sum = self.real + value.real + (self.imag + value.imag) * 1j
        # sum = self.cartesian + value

        # Getting return class
        ret = self.__get_class(sum, self.isPolar, self.isRadians)
        return ret

    def __radd__(self, value: complex):
        """Returns value + self."""
        # Chacking validity
        if not hasattr(value, "imag") or not hasattr(value, "real")\
                or not hasattr(value, "conjugate"):
            return NotImplemented
        sum = self.real + value.real + (self.imag + value.imag) * 1j
        # sum = value + self.cartesian

        # Getting return class
        ret = self.__get_class(sum, self.isPolar, self.isRadians)
        return ret

    def __sub__(self, value: complex):
        """Returns self - value."""
        # Chacking validity
        if not hasattr(value, "imag") or not hasattr(value, "real")\
                or not hasattr(value, "conjugate"):
            return NotImplemented
        diff = self.real - value.real + (self.imag - value.imag) * 1j
        # diff = self.cartesian - value

        # Getting return class
        ret = self.__get_class(diff, self.isPolar, self.isRadians)
        return ret

    def __rsub__(self, value: complex):
        """Returns value - self."""
        # Chacking validity
        if not hasattr(value, "imag") or not hasattr(value, "real")\
                or not hasattr(value, "conjugate"):
            return NotImplemented
        diff = value.real - self.real + (value.imag - self.imag) * 1j
        # diff = value - self.cartesian

        # Getting return class
        ret = self.__get_class(diff, self.isPolar, self.isRadians)
        return ret

    def __mul__(self, value: complex):
        """Returns self * value."""
        # Chacking validity
        if not hasattr(value, "imag") or not hasattr(value, "real")\
                or not hasattr(value, "conjugate"):
            return NotImplemented
        prod = self.real * value.real + self.imag * 1j * value.real + \
            self.real * value.imag * 1j + self.imag * 1j * value.imag * 1j
        # prod = self.cartesian * value

        # Getting return class
        ret = self.__get_class(prod, self.isPolar, self.isRadians)
        return ret

    def __rmul__(self, value: complex):
        """Returns value * self."""
        # Chacking validity
        if not hasattr(value, "imag") or not hasattr(value, "real")\
                or not hasattr(value, "conjugate"):
            return NotImplemented
        prod = self.real * value.real + self.imag * 1j * value.real + \
            self.real * value.imag * 1j + self.imag * 1j * value.imag * 1j
        # prod = value * self.cartesian

        # Getting return class
        ret = self.__get_class(prod, self.isPolar, self.isRadians)
        return ret

    def __truediv__(self, value: complex):
        """Returns self / value"."""
        # Chacking validity
        if not hasattr(value, "imag") or not hasattr(value, "real")\
                or not hasattr(value, "conjugate"):
            return NotImplemented
        conjugate = value.conjugate()
        denominator = (value * conjugate).real
        numerator = self.cartesian * conjugate
        quotient = numerator.real / denominator + numerator.imag / denominator
        # quotient = self.cartesian / value

        # Getting return class
        ret = self.__get_class(quotient, self.isPolar, self.isRadians)
        return ret

    def __rtruediv__(self, value: complex):
        """Returns value / self."""
        # Chacking validity
        if not hasattr(value, "imag") or not hasattr(value, "real")\
                or not hasattr(value, "conjugate"):
            return NotImplemented
        conjugate = self.conjugate(False, True)
        denominator = (self * conjugate).real
        numerator = value * conjugate
        quotient = numerator.real / denominator + numerator.imag / denominator
        # quotient = value / self.cartesian

        # Getting return class
        ret = self.__get_class(quotient, self.isPolar, self.isRadians)
        return ret

    def __abs__(self):
        """abs(self)."""
        return self.r

    def __str__(self):
        """Returns str(self)."""
        # Deciding what to return
        if self.isPolar:
            string = self.polar
        else:
            string = self.cartesian
        return str(string)

    def __repr__(self):
        """Returns repr(self)."""
        if self.isPolar:
            string = self.polar
        else:
            string = self.cartesian
        return str(string)


# Helper
def product(iterable: list | tuple, start: int = 0) -> float:
    """Returns the product of all the element of an iterable object."""
    result = 1
    for i in range(len(iterable) + start):
        result *= iterable[i]
    return result
