#!/usr/bin/env python3
"""Tests for eeepy.year_2.ecc.magnetic_materials module."""
from eeepy.year_2.ecc import magnetic_materials


class TestEddyCurrentLoss:
    """Tests for EddyCurrentLoss equation."""
    def test_power_loss(self):
        """Test for power_loss calculation."""
        P_e = magnetic_materials.EddyCurrentLoss.power_loss(0.1, 0.5, 50, 1,
                                                            1.104)
        assert P_e == round(69.0, 1)

    def test_eddy_current_constant(self):
        """Test for eddy_current_constant calculation."""
        K_e = magnetic_materials.EddyCurrentLoss.eddy_current_constant(
            69.0,
            0.5,
            50,
            1,
            1.104
        )

        assert round(K_e, 1) == 0.1

    def test_new_loss(self):
        """Test for new_loss calculation."""
        P = magnetic_materials.EddyCurrentLoss.new_loss(108, f=(60, 50))
        assert P == 75
        P = magnetic_materials.EddyCurrentLoss.new_loss(108, K_e=(0.1, 0.2))
        assert P == 216
        P = magnetic_materials.EddyCurrentLoss.new_loss(108, B_max=(0.6, 0.5))
        assert P == 75
        P = magnetic_materials.EddyCurrentLoss.new_loss(108, t=(60e-2, 50e-2))
        assert P == 75
        P = magnetic_materials.EddyCurrentLoss.new_loss(108, V=(3e-3, 6e-3))
        assert P == 216
        P = magnetic_materials.EddyCurrentLoss.new_loss(108,
                                                        K_e=(0.1, 0.2),
                                                        B_max=(0.5, 1.0),
                                                        f=(60, 50),
                                                        t=(1e-2, 5e-2),
                                                        V=(1, 0.8))
        assert round(P, 3) == 12000


class TestHysterisisLoss:
    """Tests for HysterisisLoss equation."""
    def test_power_loss(self):
        """Test for power_loss calculation."""
        P_h = magnetic_materials.HysterisisLoss.power_loss(0.1, 50, 1,
                                                           0.5, 1.0)
        assert P_h == 2.5

    def test_hysterisis_constant(self):
        """Test for hysterisis_constant calculation."""
        K_h = magnetic_materials.HysterisisLoss.hysterisis_constant(
            2.5,
            50,
            1,
            0.5,
            1.0
        )
        assert K_h == 0.1

    def test_new_loss(self):
        """Test for new_loss calculation."""
        P = magnetic_materials.HysterisisLoss.new_loss(240, f=(60, 50))
        assert P == 200
        P = magnetic_materials.HysterisisLoss.new_loss(240, K_h=(0.1, 0.2))
        assert P == 480
        P = magnetic_materials.HysterisisLoss.new_loss(240, V=(1e-3, 1.5e-3))
        assert P == 360
        P = magnetic_materials.HysterisisLoss.new_loss(240,
                                                       B_max=(0.1, 0.8),
                                                       eta=(1.0, 1.0))
        assert round(P) == 1920
        P = magnetic_materials.HysterisisLoss.new_loss(240, B_max=(0.1, 0.1),
                                                       eta=(2.0, 1.0))
        assert round(P, 3) == 2400
        P = magnetic_materials.HysterisisLoss.new_loss(240,
                                                       K_h=(0.1, 0.2),
                                                       f=(60, 50),
                                                       V=(1.5e-3, 2.0e-3),
                                                       B_max=(0.2, 0.1),
                                                       eta=(2.0, 1.0))
        assert round(P, 3) == 1333.333
