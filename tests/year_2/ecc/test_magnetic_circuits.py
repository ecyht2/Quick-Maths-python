#!/usr/bin/env python3
"""Tests for eeepy.year_2.ecc.magnetic_circuits module."""
import math

from eeepy.year_2.ecc import magnetic_circuits
from eeepy.utils.constants import mu0


class TestAmperesLaw:
    """Tests for AmperesLaw equation."""
    def test_mmf(self):
        """Test for MMF calculation."""
        N = 69
        I = 4.20

        mmf = magnetic_circuits.AmperesLaw.mmf(N, I)
        assert mmf == N * I

    def test_number_of_coils(self):
        """Test for number of coils calculation."""
        mmf = 10
        I = 2

        N = magnetic_circuits.AmperesLaw.number_of_coils(mmf, I)
        assert N == math.ceil(mmf / I)

    def test_current(self):
        """Test for current calculation."""
        mmf = 10
        N = 69

        I = magnetic_circuits.AmperesLaw.current(mmf, N)
        assert I == mmf / N


class TestHopkinsonLaw:
    """Tests for HopkinsonLaw equation."""
    phi = 0.0012
    mmf = 100 * 4
    r = 3.37 * 105

    def test_reluctance(self):
        """Test for reluctance calculation."""
        R = magnetic_circuits.HopkinsonLaw.reluctance(self.mmf, self.phi)
        assert R == self.mmf / self.phi

    def test_total_flux(self):
        """Test for total flux calculation."""
        flux = magnetic_circuits.HopkinsonLaw.total_flux(self.mmf, self.r)
        assert flux == self.mmf / self.r

    def test_mmf(self):
        """Test for mmf calculation."""
        mmf = magnetic_circuits.HopkinsonLaw.mmf(self.phi, self.r)
        assert mmf == self.phi * self.r


class TestReluctanceEquation:
    """Tests for ReluctanceEquation calculations."""
    l = 0.65
    mu = mu0 * 1000
    A = 2.5e-3

    def test_reluctance(self):
        """Test for reluctance calculation."""
        reluctance = magnetic_circuits.ReluctanceEquation\
                                      .reluctance(self.l, self.mu, self.A)
        assert round(reluctance, -2) == 2.069e5

    def test_reluctance_relative_permeability(self):
        """Test for reluctance using relative permeability calculation."""
        reluctance = magnetic_circuits.ReluctanceEquation\
            .reluctance_relative_permeability(self.l, 1000, self.A)
        assert round(reluctance, -2) == 2.069e5

    def test_length(self):
        """Test for length calculation."""
        l = magnetic_circuits.ReluctanceEquation\
            .length(2.069e5, self.mu, self.A)
        assert round(l, 2) == self.l
