#!/usr/bin/env python3
"""This file contains test for coordinate_system module of eeepy.year_2.mmt"""
import math
from eeepy.year_2.mmt.coordinate_system import Point, Point3D


class TestPoint:
    """Tests for Point class."""
    def test_constructor(self):
        """Tests for Point class constructor."""
        point = Point(2, 2)

        # Cartesian Values
        assert point.x == 2
        assert point.y == 2

        # Polar Values
        assert point.r == math.sqrt(8)
        assert point.phi == math.pi / 4

    def test_from_polar(self):
        """Tests for Point class from_polar classmethod."""
        point = Point.from_polar(2, math.pi / 4)

        # Polar Values
        assert point.r == 2
        assert point.phi == math.pi / 4

        # Cartesian Values
        assert point.x == math.sqrt(2)
        assert point.y == 1.414213562373095

    def test_methods(self):
        """Tests for Point class methods."""
        point = Point(1, -math.sqrt(3))

        assert point.cartesian() == (1.0, -math.sqrt(3))
        assert point.polar() == (1.9999999999999998, - math.pi / 3)
        assert abs(point) == 1.9999999999999998
        assert repr(point) == str({
            "cartesian": point.cartesian(),
            "polar": point.polar(),
        })


class TestPoint3D:
    """Tests for Point3D class."""
    def test_constructor(self):
        """Tests for Point3D class constructor."""
        x = 69
        y = 420
        z = 69420
        point = Point3D(x, y, z)

        # Cartesian Values
        assert point.x == x
        assert point.y == y
        assert point.z == z

        # Cylindrical Values
        assert point.r_cylindrical == math.sqrt(x**2 + y**2)
        assert point.phi == math.atan2(y, x)

        # Spherical Values
        assert point.r_spherical == math.sqrt(x**2 + y**2 + z**2)
        assert point.theta == math.acos(z / point.r_spherical)

    def test_from_cylindrical(self):
        """Tests for Point3D class from_cylindrical classmethod."""
        r = 69
        phi = math.pi / 3
        z = 420
        point = Point3D.from_cylindrical(r, phi, z)

        # Cylindrical Values
        assert point.r_cylindrical == r
        assert point.phi == phi
        assert point.z == z

        # Cartesian Values
        x = r * math.cos(phi)
        y = r * math.sin(phi)
        assert point.x == x
        assert point.y == y

        # Spherical Values
        r = math.sqrt(x**2 + y**2 + z**2)
        assert point.r_spherical == r
        assert point.theta == math.acos(z / r)

    def test_from_spherical(self):
        """Tests for Point3D class from_spherical classmethod."""
        r = 69
        phi = math.pi / 3
        theta = math.pi / 4
        point = Point3D.from_spherical(r, phi, theta)

        # Spherical Values
        assert point.r_spherical == r
        assert point.phi == phi
        assert point.theta == theta

        # Cartesian Values
        x = r * math.sin(theta) * math.cos(phi)
        y = r * math.sin(theta) * math.sin(phi)
        z = r * math.cos(theta)
        assert point.x == x
        assert point.y == y
        assert point.z == z

        # Cylindrical Values
        assert point.r_cylindrical == math.sqrt(x**2 + y**2)

    def test_methods(self):
        """Tests for Point3D class methods."""
        x = 69
        y = 420
        z = 69420
        point = Point3D(x, y, z)

        r_cyl = math.sqrt(x**2 + y**2)
        phi = math.atan2(y, x)
        r_sph = math.sqrt(x**2 + y**2 + z**2)
        theta = math.acos(z / point.r_spherical)

        assert point.cartesian() == (x, y, z)
        assert point.cylindrical() == (r_cyl, phi, z)
        assert point.spherical() == (r_sph, phi, theta)
