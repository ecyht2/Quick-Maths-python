#!/usr/bin/env python3
from math import exp, factorial
from math import log10, log1p, log2
from statistics import mean,median,mode
from statistics import quantiles, variance, stdev
from math import sqrt, atan, sin, cos, exp

# Helper
def product(iterable, start = 0):
    """

    """
    result = 1
    for i in range(len(iterable) + start):
        result *= iterable[i]

    return result
