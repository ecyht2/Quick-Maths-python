#!/usr/bin/env python3
"""Useful classes and functions to perform basic tasks."""
import cmath
import itertools
from collections.abc import Iterable
from typing import TypeVar

import numpy as np

_T = TypeVar("_T", bound=complex)


class Complex(complex):
    """Create a complex number from a real part and an imaginary part."""

    def __new__(cls, z: complex):
        return super().__new__(cls, z)

    def __init__(self, z: complex):
        self._is_polar = False
        self._is_radians = True
        # Polar Attributes
        self._r = abs(z)
        self._phi = cmath.phase(z)

    @property
    def is_polar(self) -> bool:
        """Check if the complex number is in polar form."""
        return self._is_polar

    @property
    def is_cartesian(self) -> bool:
        """Check if the complex number is in cartesian form."""
        return not self._is_polar

    @property
    def r(self) -> float:
        """The magnitude of the complex number."""
        return self._r

    @property
    def phi(self) -> float:
        """The phase/angle of the complex number."""
        return self._phi

    @property
    def is_radians(self) -> bool:
        """Check if the angle is in radians or not."""
        return self._is_radians

    def cartesian(self) -> tuple[float, float]:
        """Gets the complex number in cartesian form.

        :returns: A tuple where the first element is the real part and the second
            element is the imaginary part of the complex number.
        """
        return (self.real, self.imag)

    def polar(self) -> tuple[float, float]:
        """Gets the complex number in polar form.

        :returns: A tuple where the first element is the magnitude and the second
            element is the phase/angle of the complex number.
        """
        return (self.r, self.phi)

    @classmethod
    def from_polar(cls, r: float, phi: float, radian: bool = True):
        """Creates a new complex class from polar form.

        :param r: The magnitude of the complex number.
        :param phi: The phase/angle of the complex number in radians if `radian` is set
            to True and degrees if set to False (defaults to radians).
        :param radian: The specified angle in `phi` is in radians or degrees.
        """
        if not radian:
            phi = np.radians(phi)
        z = cmath.rect(r, phi)
        z = cls(z)
        z._is_polar = True
        z._is_radians = radian
        return z

    @classmethod
    def from_polar_tuple(cls, Z: tuple, radian: bool = True):
        """Creates a Complex class from polar form in a tuple.

        This classmethod treats the first item in the tuple as the radius from
        the origin and the second item as the angle from the x-axis. If the
        parameter radian is True the angle is in radians, otherwise it is in
        degrees.
        """
        z = cmath.rect(Z[0], np.radians(Z[1])) if not radian else cmath.rect(Z[0], Z[1])
        z = cls(z)
        z._is_radians = radian
        z._is_polar = True
        return z

    @classmethod
    def from_cartesian(cls, Z: complex, radian: bool = True):
        """Creates a Complex class from cartesian form."""
        z = cls(Z)
        z._is_radians = radian
        z._is_polar = False
        return z

    def conjugate(self, polar: bool | None = None, radian: bool | None = None):
        """Return the complex conjugate of its argument.

        :param polar:
        :param radian:
        :returns:
        """
        if (polar is None and self.is_polar) or polar:
            return self.from_polar_tuple(
                cmath.polar(super().conjugate()),
                self.is_radians if radian is None else radian,
            )

        return self.__class__(super().conjugate())

    def __repr__(self):
        """Returns repr(self)."""
        # Deciding what to return
        if self.is_polar:
            return str(self.polar())

        return super().__repr__()


# Helper
def product(iterable: Iterable[_T], /, start: int = 0) -> _T:
    """Returns the product of all the element of an iterable object.

    :param iterable: An iterable object.
    :param start: The nth element to start the product.
    :returns: The product of the iterable.
    """
    iterator = itertools.islice(iterable, start, None)
    result = next(iterator)
    for i in iterator:
        result *= i  # type: ignore
    return result


def quadratic(a: float, b: float, c: float) -> tuple[Complex, Complex]:
    """Find the roots of a quadratic equation.

    The quardatic formula:

    $x = \\frac{-b \\pm \\sqrt{b^{2} - 4 a c}}{2 a}$

    The values of the parameters follows the equation:

    $a x^{2} + b x + c$

    :param a: The first value of the quadratic formula.
    :param b: The second value of the quadratic formula.
    :param c: The third value of the quadratic formula.
    :returns: The roots of the equation.
    """
    z1 = -b + cmath.sqrt(b**2 - 4 * a * c)
    z1 /= 2 * a
    z2 = -b - cmath.sqrt(b**2 - 4 * a * c)
    z2 /= 2 * a
    return (Complex(z1), Complex(z2))
