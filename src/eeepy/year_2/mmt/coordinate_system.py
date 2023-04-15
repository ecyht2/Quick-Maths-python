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
        self._phi = math.atan2(y, x)

    @property
    def x(self) -> float:
        """The x coordinate of the point."""
        return self._x

    @property
    def y(self) -> float:
        """The y coordinate of the point."""
        return self._y

    @property
    def r(self) -> float:
        """The distance from origin."""
        return self._r

    @property
    def phi(self) -> float:
        """The azimuthal angle of the point."""
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

    def cartesian(self) -> tuple:
        """The Point in cartesian form.

        :return: (x, y)
        """
        return (self.x, self.y)

    def polar(self) -> tuple:
        """The Point in polar form.

        :return: (r, ɸ)
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
    """A Point in a 3D plane."""
    def __init__(self, x: float, y: float, z: float):
        """Creates a Point in a 2D plane.

        :param x: The x coordinate of the point.
        :param y: The y coordinate of the point.
        :param z: The z coordinate of the point.
        """
        self._x = float(x)
        self._y = float(y)
        self._z = float(z)

    @property
    def x(self):
        """The x coordinate of the point."""
        return self._x

    @property
    def y(self):
        """The y coordinate of the point."""
        return self._y

    @property
    def z(self):
        """The z coordinate of the point."""
        return self._z

    @property
    def r_spherical(self):
        """The distance from origin."""
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    @property
    def theta(self):
        """The zenith angle of the point."""
        return math.acos(self.z / self.r_spherical)

    @property
    def phi(self):
        """The azimuthal angle of the point."""
        return math.atan2(self.y, self.x)

    @property
    def r_cylindrical(self):
        """The radial distance of the point."""
        return math.sqrt(self.x**2 + self.y**2)

    def cartesian(self) -> tuple:
        """The Point in cartesian form.

        :return: (x, y, z)
        """
        return (self.x, self.y, self.z)

    def spherical(self) -> tuple:
        """The Point in spherical form.

        :return: (r, ɸ, θ)
        """
        return (self.r_spherical, self.phi, self.theta)

    def cylindrical(self) -> tuple:
        """The Point in cylindrical form.

        :return: (r, ɸ, z)
        """
        return (self.r_cylindrical, self.phi, self.z)

    @classmethod
    def from_spherical(cls, r: float, phi: float, theta: float):
        """Creates a Point in a 3D plane from spherical coordinates.

        :param r: The distance from origin.
        :param phi: The azimuthal angle of the point.
        :param theta: The zenith angle of the point.
        """
        x = r * math.sin(theta) * math.cos(phi)
        y = r * math.sin(theta) * math.sin(phi)
        z = r * math.cos(theta)
        return cls(x, y, z)

    @classmethod
    def from_cylindrical(cls, r: float, phi: float, z: float):
        """Creates a Point in a 3D plane from cylindrical coordinates.

        :param r: The radial distance of the point.
        :param phi: The azimuth angle of the point.
        :param z: The z coordinate of the point.
        """
        x = r * math.cos(phi)
        y = r * math.sin(phi)
        return cls(x, y, z)
