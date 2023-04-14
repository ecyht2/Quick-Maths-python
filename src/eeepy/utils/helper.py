#!/usr/bin/env python3
"""Useful classes and functions to perform basic tasks."""
import cmath
from math import degrees, radians
from typing import Union


class Complex():
    """Create a complex number from a real part and an imaginary part."""
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

    def polar(self) -> tuple:
        """Returns the polar form of the Complex number."""
        return cmath.polar(self.cartesian)

    @classmethod
    def from_polar(cls, r: float, phi: float, radian: bool = True):
        """Creates a Complex class from polar form."""
        if not radians:
            phi = radians(phi)
        cls.isRadians = radian
        cls.isPolar = True
        Z = cmath.rect(r, phi)
        return cls(Z)

    @classmethod
    def from_polar_tuple(cls, Z: tuple, radian: bool = True):
        """Creates a Complex class from polar form in a tuple.

        This classmethod treats the first item in the tuple as the radius from
        the origin and the second item as the angle from the x-axis. If the
        parameter radian is True the angle is in radians, otherwise it is in
        degrees.
        """
        if not radians:
            Z = cmath.rect(Z[0], radians(Z[1]))
        else:
            Z = cmath.rect(Z[0], Z[1])
        cls.isRadians = radian
        cls.isPolar = True
        return cls(Z)

    @classmethod
    def from_cartesian(cls, Z: complex, radian: bool = True):
        """Creates a Complex class from cartesian form."""
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
        if polar is None:
            polar = self.isPolar
        # Manual overide radians
        if radian is None:
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
        total = self.real + value.real + (self.imag + value.imag) * 1j

        # Getting return class
        ret = self.__get_class(total, self.isPolar, self.isRadians)
        return ret

    def __radd__(self, value: complex):
        """Returns value + self."""
        # Chacking validity
        if not hasattr(value, "imag") or not hasattr(value, "real")\
                or not hasattr(value, "conjugate"):
            return NotImplemented
        total = self.real + value.real + (self.imag + value.imag) * 1j

        # Getting return class
        ret = self.__get_class(total, self.isPolar, self.isRadians)
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
            string = self.polar()
        else:
            string = self.cartesian
        return str(string)

    def __repr__(self):
        """Returns repr(self)."""
        if self.isPolar:
            string = self.polar()
        else:
            string = self.cartesian
        return str(string)


# Helper
def product(iterable: Union[list, tuple], start: int = 0) -> float:
    """Returns the product of all the element of an iterable object."""
    # list | tuple
    result = 1
    for i in range(len(iterable) + start):
        result *= iterable[i]
    return result
