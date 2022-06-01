#!/usr/bin/env python3
from math import exp, factorial, sqrt
from math import log10, log1p, log2
from statistics import mean,median,mode
from statistics import quantiles, variance, stdev
from math import acos, atan, asin
from math import tan, sin, cos
from math import radians, degrees
import cmath

# Helper
def product(iterable: list | tuple, start: int = 0) -> float:
    """Returns the product of all the element of an iterable object."""
    result = 1
    for i in range(len(iterable) + start):
        result *= iterable[i]
    return result
