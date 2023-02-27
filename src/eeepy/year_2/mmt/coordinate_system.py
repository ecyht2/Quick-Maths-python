#!/usr/bin/env python3
"""Functions and Equations that helps with calculating values related to
 coordinate system."""
import math


class Point:
    """A Point in a 2D plane."""
    def __init__(self, x: float, y: float):
        """Creates a Point in a 2D plane.

        :param x: The x coordinate of the point.
        :param y: The y coordinate of the point.
        """
        self._x = float(x)
        self._y = float(y)
        self._r = (x**2 + y**2)**0.5
        self._phi = math.atan2(y / x)

    @property
    def x(self) -> float:
        """The x coordinate of the 2D point."""
        return self._x

    @property
    def y(self) -> float:
        return self._y

    @property
    def r(self) -> float:
        return self._r

    @property
    def phi(self) -> float:
        return self._phi

    @classmethod
    def from_polar(cls, r: float, phi: float):
        """Creates a Point in a 2D plane using polar coordinates.

        :param r: The radius from the origin.
        :param phi: The angle from the positive x axis.
        """
        x = r * math.cos(phi)
        y = r * math.sin(phi)
        point = cls(x, y)
        return point

    def cartesian(self) -> tuple[float, float]:
        """The Point in cartesian form.

        :return: (x, y)
        """
        return (self.x, self.y)

    def polar(self) -> tuple[float, float]:
        """The Point in polar form.

        :return: (r, phi)
        """
        return (self.r, self.phi)

    def __repr__(self) -> str:
        """The representation of a Point."""
        return str({
            "cartesian": self.cartesian(),
            "polar": self.polar(),
        })

    def __abs__(self) -> float:
        """Returns the r of the Point."""
        return self.r


class Point3D:
    def __init__(self):
        ...
