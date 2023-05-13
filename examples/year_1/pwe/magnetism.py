#!/usr/bin/env python3
"""Example usage of magnetic circuits functions.

This example is based on Q1c of 2021-2022 final exams.
"""
import math

from eeepy.year_1.pwe.magnetism import FluxIntensity, FluxDensity, flux


if __name__ == '__main__':
    I = 1.8
    a = 3.5e-2
    b = 4.5e-2
    mu_r = 1e3
    N = 300

    # i
    r = (a + b) / 2
    print(f"{r=}")
    l = 2 * math.pi * r
    print(f"{l=}")

    # ii
    A = (b - a)**2
    print(f"{A=}")

    # iii
    H = FluxIntensity.H_field(N, I, l)
    print(f"{H=}")

    # iv
    B = FluxDensity.B_field(H, mu_r)
    print(f"{B=}")

    # v
    phi = flux(B, A)
    print(f"{phi=}")

    # vi
    B = 1.5
    H = FluxDensity.flux_intensity(B, mu_r)
    I = FluxIntensity.current(H, N, l)
    print(f"{I=}")
