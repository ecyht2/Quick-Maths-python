#!/usr/bin/env python3
"""This file contains tests for numerical_integration module of
 eeepy.year_2.mmt."""
from eeepy.year_2.mmt import numerical_integration


def test_riemann():
    """Test for Reimann Sum Approximation."""
    def fx(x: float) -> float:
        return x**3 / 20 - x**2 / 1.7 + 2 * x + 2

    res = numerical_integration.riemann_sum(fx, 6)
    assert res == 20.89705882352941

    res = numerical_integration.riemann_sum(fx, 12, 0.5)
    assert res == 21.406617647058827

    res = numerical_integration.riemann_sum(fx, 12, 0.5, 1)
    assert res == 22.727941176470587


def func(x: float) -> float:
    """Shared fx function.

    f(x) = 1 / (x + 1)
    """
    return 1 / (x + 1)


def test_trepezoid():
    """Test for Trapezoid Rule Approximation."""
    res = numerical_integration.trapezoid_rule(func, 6, 1 / 6)
    assert round(res, 4) == 0.6949


def test_simpson_13():
    """Test for Simpson's 1/3 Rule Approximation."""
    res = numerical_integration.simpson_13(func, 6, 1 / 6)
    assert round(res, 6) == 0.693170


def test_simpson_38():
    """Test for Simpson's 3/8 Rule Approximation."""
    res = numerical_integration.simpson_38(func, 6, 1 / 6)
    assert round(res, 6) == 0.693195
