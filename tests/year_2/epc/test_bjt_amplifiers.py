#!/usr/bin/env python3
"""Tests for the bjt_amplifiers module of year_2 epc."""
from eeepy.year_2.epc import bjt_amplifiers


def test_g_m():
    g_m = bjt_amplifiers.calculate_g_m(0.0009486540969998821)
    assert g_m == 0.03648669603845701


def test_r_o():
    r_o = bjt_amplifiers.calculate_r_o(100, 0.0009486540969998821)
    assert r_o == 105412.49999999991


def test_r_pi():
    r_pi = bjt_amplifiers.calculate_r_pi(9.486540969998822e-06)
    assert r_pi == 2740.7249999999976


def test_ssparam_dataclass():
    ss_param = bjt_amplifiers.SSParameters(69, 420e-6, 69420e-9)
    assert ss_param.g_m() == 0.016153846153846154
    assert ss_param.r_pi() == 374.53183520599254
    assert ss_param.r_o() == 164285.7142857143
