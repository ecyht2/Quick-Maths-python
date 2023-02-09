#!/usr/bin/env python3
"""Functions and Equations that helps with calculating values related to
 coordinate system."""


class Point:
    def __init__(self, x: float, y: float):
        self.x = float(x)
        self.y = float(y)

    def from_polar(cls, r: float, phi: float):
        return cls(0, 0, 0)


class Point3D:
    def __init__(self):
        ...
