#!/usr/bin/env python3
"""Tests for the power_devices_losses module of year_2 ecc."""
from eeepy.year_2.ecc import power_devices_losses


class TestConductionLoss:
    """Tests for conduction losses equations."""
    def test_mosfet(self) -> None:
        """Testing MOSFET conduction loss calculations."""
        R_q: float = 120e-3
        i_q_rms: float = 5.860887304837042
        loss: float = power_devices_losses.ConductionLoss.mosfet(
            R_q, i_q_rms)
        assert loss == 4.122

    def test_diode(self) -> None:
        """Testing diode conduction loss calculations."""
        V_on_d: float = 0.7
        i_d: float = 8.25
        R_d: float = (1.7 - 0.7) / 20
        i_d_rms: float = 11.32843031197762
        loss: float = power_devices_losses.ConductionLoss.diode(
            V_on_d, i_d, R_d, i_d_rms)
        assert loss == 12.191666666666666


class TestSwitchingLoss:
    """Tests for switching losses equations."""
    V_s: int = 100
    f: float = 200e3
    t_on: float = 50e-9
    i_t: int = 5
    i_p: int = 12
    t_off: float = 100e-9

    def test_on_energy(self) -> None:
        """Testing turn on energy loss."""
        loss: float = power_devices_losses.SwitchingLoss.energy_turn_on(
            self.i_t, self.V_s, self.t_on
        )
        assert loss == 1.2499999999999999e-05

    def test_on_power_1(self) -> None:
        """Testing turn on power loss."""
        loss: float = power_devices_losses.SwitchingLoss.power_turn_on(
            self.i_t, self.V_s, self.t_on, self.f
        )
        assert loss == 2.5

    def test_on_power_2(self) -> None:
        """Testing turn on power loss using energy loss."""
        energy: float = power_devices_losses.SwitchingLoss.energy_turn_on(
            self.i_t, self.V_s, self.t_on
        )
        loss: float = power_devices_losses.SwitchingLoss.power_turn_on_e(
            energy, self.f
        )
        assert loss == 2.5

    def test_off_energy(self) -> None:
        """Testing turn off energy loss."""
        loss: float = power_devices_losses.SwitchingLoss.energy_turn_off(
            self.i_p, self.V_s, self.t_off
        )
        assert loss == 5.9999999999999995e-05

    def test_off_power_1(self) -> None:
        """Testing turn off power loss."""
        loss: float = power_devices_losses.SwitchingLoss.power_turn_off(
            self.i_p, self.V_s, self.t_off, self.f)
        assert loss == 11.999999999999998

    def test_off_power_2(self) -> None:
        """Testing turn off power loss using energy loss."""
        energy: float = power_devices_losses.SwitchingLoss.energy_turn_off(
            self.i_p, self.V_s, self.t_off
        )
        loss: float = power_devices_losses.SwitchingLoss.power_turn_off_e(
            energy, self.f)
        assert loss == 11.999999999999998


def test_rms() -> None:
    """Test for RMS current equation."""
    i1: int = 12
    i2: int = 5
    d = 0.45

    rms = power_devices_losses.I_RMS(i1, i2, d)
    assert rms == 5.860887304837042
