#!/usr/bin/env python3
"""Functions and Equations relating to solving integration numerically."""
from typing import Callable


def riemann_sum(fx: Callable[[float], float], n: int,
                h: float = 1.0, initial: float = 0) -> float:
    """Calculates the area under the curve using the Riemann sum
    approximation.

    :param fx: The curve function.
    :param n: The number of sample.
    :param h: The distance between each sample.
    :param initial: The initial x value to start with.
    """
    fx_h = []
    for i in range(n):
        fx_h.append(fx(i * h + initial) * h)
    return sum(fx_h)


def trapezoid_rule(fx: Callable[[float], float], n: int,
                   h: float = 1.0, initial: float = 0) -> float:
    """Calculates the area under the curve using the Trapezoid Rule
    approximation.

    :param fx: The curve function.
    :param n: The number of sample.
    :param h: The distance between each sample.
    :param initial: The initial x value to start with.
    """
    fx_h = []
    for i in range(n):
        y1 = fx(i * h + initial)
        y2 = fx((i + 1) * h + initial)
        fx_h.append(0.5 * (y1 + y2) * h)
    return sum(fx_h)


def simpson_13(fx: Callable[[float], float], n: int,
               h: float = 1.0, initial: float = 0) -> float:
    """Calculates the area under the curve using the Simpson's 1/3 Rule
    approximation.

    :param fx: The curve function.
    :param n: The number of sample.
    :param h: The distance between each sample.
    :param initial: The initial x value to start with.
    """
    f0 = fx(initial)
    fx_odd = []
    fx_even = []
    fn = fx(n * h + initial)

    for i in range(1, n):
        val = fx(i * h + initial)

        # Appending Values
        # Divisible by 3
        if i % 2 == 0:
            fx_even.append(val)
        # Non Divisible by 3
        else:
            fx_odd.append(val)

    summation = f0 + 4 * sum(fx_odd) + 2 * sum(fx_even) + fn
    return h / 3 * summation


def simpson_38(fx: Callable[[float], float], n: int,
               h: float = 1.0, initial: float = 0) -> float:
    """Calculates the area under the curve using the Simpson's 3/8 Rule
    approximation.

    :param fx: The curve function.
    :param n: The number of sample.
    :param h: The distance between each sample.
    :param initial: The initial x value to start with.
    """
    f0 = fx(initial)
    fx_div_3 = []
    fx_non_div_3 = []
    fn = fx(n * h + initial)

    for i in range(1, n):
        val = fx(i * h + initial)

        # Appending Values
        # Divisible by 3
        if i % 3 == 0:
            fx_div_3.append(val)
        # Non Divisible by 3
        else:
            fx_non_div_3.append(val)

    summation = f0 + 2 * sum(fx_div_3) + 3 * sum(fx_non_div_3) + fn
    return 3 * h / 8 * summation
