#!/usr/bin/env python3
"""Tests for eeepy.year_2.ecc.dc_machine module."""
import math

import numpy as np
from eeepy.year_2.ecc import dc_machine


class TestPoleFlux:
    """Tests for PoleFlux equation."""
    B = 1
    R = 0.05
    p = 4
    l = 0.12
    lambda_f = 68
    lambda_p = 90

    def test_flux(self):
        """Test for flux calculation."""
        flux = dc_machine.PoleFlux.flux(
            self.B,
            self.R,
            self.p / 2,
            self.l,
            self.lambda_f,
            self.lambda_p
        )

        assert float(np.format_float_scientific(flux, 3)) == 7.121e-3

    def test_flux_poles(self):
        """Test for flux_poles calculation."""
        flux = dc_machine.PoleFlux.flux_poles(
            self.B,
            self.R,
            self.p,
            self.l,
            self.lambda_f,
            self.lambda_p
        )

        assert float(np.format_float_scientific(flux, 3)) == 7.121e-3


class TestInducedTorque:
    """Tests for InducedTorque equation."""
    N = 32
    Z = 12
    I = 40
    flux = 7.121e-3
    T = 17.41

    def test_torque(self):
        """Test for torque calculation."""
        T_d = dc_machine.InducedTorque.torque(self.N, self.Z,
                                              self.I, self.flux)
        assert round(T_d, 2) == 17.41

    def test_torque_electrical_loading(self):
        """Test for torque_electrical_loading calculation."""
        T_d = dc_machine.InducedTorque.torque_electrical_loading(
            0.1,
            1,
            1,
            60,
            90
        )
        assert T_d, 2 == 2 * 0.1 * 1 * 1 * 60 / 90

    def test_armature_current(self):
        """Test for armature_current calculation."""
        I = dc_machine.InducedTorque.armature_current(self.T, self.N,
                                                      self.Z, self.flux)
        assert round(I, 0) == self.I


class TestEMFInduced:
    """Tests for EMFInduced equation."""
    N = 32
    Z = 12
    flux = 7.121e-3
    f = 172.31 / 2 / math.pi
    emf = 74.99

    def test_emf_induced(self):
        """Test for emf_induced calculation."""
        emf = dc_machine.EMFInduced.emf_induced(self.N, self.Z,
                                                self.flux, self.f)
        assert round(emf, 2) == self.emf

    def test_mechanical_frequency(self):
        """Test for mechanical_frequency calculation."""
        f_m = dc_machine.EMFInduced.mechanical_frequency(self.emf, self.N,
                                                         self.Z, self.flux)
        assert round(f_m, 2) == round(self.f, 2)
