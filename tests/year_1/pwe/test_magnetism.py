#!/usr/bin/env python3
"""Tests for eeepy.year_1.pwe.magnetism module."""
import math

from eeepy.year_1.pwe import magnetism
from eeepy.utils.constants import mu0


class TestFluxIntensity:
    """Tests for flux intensity equation."""
    r = 4.5e-2
    N = 500
    I = 2

    def test_h_field_circular(self):
        """Test for H_field_circular staticmethod."""
        H = magnetism.FluxIntensity.H_field_circular(self.N, self.I, self.r)
        assert H == self.N * self.I / (2 * math.pi * self.r)

    def test_h_field(self):
        """Test for H_field staticmethod."""
        l = 2 * math.pi * self.r
        H = magnetism.FluxIntensity.H_field(self.N, self.I, l)
        assert H == self.N * self.I / l

    def test_current(self):
        """Test for current staticmethod."""
        l = 2 * math.pi * self.r
        H = 5
        I = magnetism.FluxIntensity.current(H, self.N, l)
        assert I == H * l / self.N


class TestFluxDensity:
    """Tests for flux density equation."""
    mu_r = 1000

    def test_b_field(self):
        """Test for B_field staticmethod."""
        H = 1
        B = magnetism.FluxDensity.B_field(H, self.mu_r)
        assert B == self.mu_r * H * mu0

    def test_b_field_from_current(self):
        """Test for B_field_from_current staticmethod."""
        N = 2
        I = 2
        l = 2
        B = magnetism.FluxDensity.B_field_from_current(N, I, l, self.mu_r)
        assert B == self.mu_r * mu0 * N * I / l

    def test_flux_intensity(self):
        """Test for flux_intensity staticmethod."""
        B = 69
        H = magnetism.FluxDensity.flux_intensity(B, self.mu_r)
        assert H == B / (self.mu_r * mu0)


def test_magnetic_force_wire():
    """Test for magnetic_for_wire function."""
    I = 1
    B = 1.5
    l = 5e-2
    F = magnetism.magnetic_force_wire(I, B, l)
    assert F == abs(I * B * l)


def test_flux():
    """Test for flux function."""
    B = 1
    A = 6.9
    assert magnetism.flux(B, A) == B * A


def test_flux_linkage():
    """Test for flux linkage function."""
    phi = 3
    N = 100
    assert magnetism.flux_linkage(N, phi) == N * phi


def test_inductance():
    """Test for inductance function."""
    I = 1
    N = 100
    phi = 6.9
    linkage = magnetism.flux_linkage(N, phi)

    assert magnetism.inductance(I, linkage) == linkage / I
    assert magnetism.inductance(I, N=N, phi=phi) == N * phi / I


def test_flux_linkage_two_coil():
    """Test for flux_linkage_two_coil function."""
    L = 2
    I1 = 1
    I2 = 3
    M = 2

    assert magnetism.flux_linkage_two_coil(L, I1, I2, 2) ==\
        L * I1 + M * I2
    assert magnetism.flux_linkage_two_coil(L, I1, I2, 2, False) ==\
        L * I1 - M * I2


def test_volatage_two_coil():
    """Test for voltage_two_coil function."""
    R = 3
    I1 = 1
    L = 1e-6
    dI1 = 0.1
    I2 = 3
    M = 1

    V = R * I1

    assert magnetism.voltage_two_coil(R, I1, L, dI1, I2, M) ==\
        V + (L * dI1 + M * I2)
    assert magnetism.voltage_two_coil(R, I1, L, dI1, I2, M, False) ==\
        V + (L * I1 - M * I2)


def voltage_two_coil_phasor():
    """Test for voltage_two_coil_phasor function."""
    I = [1, 2]
    R = 3
    L = 1e-6
    f = 50
    M = 1

    omega = 2 * math.pi * f

    Z = R + magnetism.imp_ind(L, omega=omega)
    V = I[0] * Z

    assert magnetism.voltage_two_coil_phasor(R, L, I, M, omega, f) ==\
        V + (M * 1j * omega * I[1])
    assert magnetism.voltage_two_coil_phasor(R, L, I, M, omega, f, False) ==\
        V - (M * 1j * omega * I[1])
