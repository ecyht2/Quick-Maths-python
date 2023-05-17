#!/usr/bin/env python3
"""Functions and equations related to solving Ordinary Differential
 Equation(ODE) numerically."""
from typing import Callable, Tuple, Union, List

from ..coordinate_system import Point


def eulers_forward(ode: Callable[[float, float], float], h: float, n: int,
                   initial: Union[Tuple[float, float], Point]) -> List[float]:
    """Solves the ODE using forward difference (FDM) Euler's method.

    :param ode: The odidary differential equation. The first parameter of the
        ode function is the x and second is the y value.
    :param h: The step size.
    :param n: The number of iterations.
    :param initial: The initial condition where the first value of the tuple is
        the x and the second is the y.
    :returns: The list of solutions to the ODE at each step.
    """
    result = []

    result.append(initial[1])
    x0 = initial[0]
    y0 = initial[1]

    for _ in range(n):
        y = y0 + h * ode(x0, y0)
        result.append(y)
        y0 = y
        x0 += h

    return result


def heun_method(ode: Callable[[float, float], float], h: float, n: int,
                initial: Union[Tuple[float, float], Point]) -> List[float]:
    """Solves the ODE using Heun method.

    :param ode: The odidary differential equation. The first parameter of the
        ode function is the x and second is the y value.
    :param h: The step size.
    :param n: The number of interations.
    :param initial: The initial condition where the first value of the tuple is
        the x and the second is the y.
    :returns: The list of solutions to the ODE at each step.
    """
    result = []

    result.append(initial[1])
    x0 = initial[0]
    y0 = initial[1]

    for _ in range(n):
        K1 = ode(x0, y0)
        K2 = ode(x0 + h, y0 + h * K1)
        y = y0 + h / 2 * (K1 + K2)
        result.append(y)
        y0 = y
        x0 += h

    return result


def rk4_method(ode: Callable[[float, float], float], h: float, n: int,
               initial: Union[Tuple[float, float], Point]) -> List[float]:
    """Solves the ODE using RK4 method.

    :param ode: The odidary differential equation. The first parameter of the
        ode function is the x and second is the y value.
    :param h: The step size.
    :param n: The number of iterations.
    :param initial: The initial condition where the first value of the tuple is
        the x and the second is the y.
    :returns: The list of solutions to the ODE at each step.
    """
    result = []

    result.append(initial[1])
    x0 = initial[0]
    y0 = initial[1]

    for _ in range(n):
        K1 = ode(x0, y0)
        K2 = ode(x0 + h / 2, y0 + K1 * h / 2)
        K3 = ode(x0 + h / 2, y0 + K2 * h / 2)
        K4 = ode(x0 + h, y0 + K3 * h)
        y = y0 + h / 6 * (K1 + 2 * K2 + 2 * K3 + K4)
        result.append(y)
        y0 = y
        x0 += h

    return result
