#!/usr/bin/env python3
"""Tests for eeepy.year_2.mmt.numerical_integration.ode module."""
import numpy as np
from eeepy.year_2.mmt.numerical_integration.ode import (eulers_forward,
                                                        heun_method,
                                                        rk4_method)


def ode(x: float, y: float) -> float:
    """Shared ODE equation used for tests."""
    return - np.sin(x) / y


def test_eulers_forward():
    """Test for eulers_forward calulation."""
    desired = [1.00, 1.00, 0.96, 0.88, 0.75, 0.56]

    result = eulers_forward(ode, 0.2, 5, (0, 1.0))
    result = list(map(lambda x: round(x, 2), result))

    assert result == desired


def test_heun_method():
    """Test for eulers_forward calulation."""
    desired = [1.00, 0.980, 0.918, 0.808, 0.631, 0.309]

    result = heun_method(ode, 0.2, 5, (0, 1.0))
    result = list(map(lambda x: round(x, 3), result))

    assert result == desired


def test_rk4_method():
    """Test for rk4_method calulation."""
    desired = [1.000, 0.980, 0.918, 0.807, 0.627, 0.283]

    result = rk4_method(ode, 0.2, 5, (0, 1.0))
    result = list(map(lambda x: round(x, 3), result))

    assert result == desired
