#!/usr/bin/env python
"""Example script that uses the power_devices_losses module.

The questions are from the power devices losses example question.
"""
import math
from eeepy.year_2.ecc import power_devices_losses


def run(*args, **kwargs):
    """Function that does all the calculations."""
    d = 0.45
    i_p = 20
    i_t = 10
    i_p1 = 12
    i_t1 = 5
    V_s = 100
    I_o = 15
    V_o = 0.5 * d * V_s
    f = 200e3

    # MOSFET
    # i_q = (i_t1 + i_p1) / 2 * d
    i_q_rms = math.sqrt(d * (i_p1**2 + i_p1 * i_t1 + i_t1**2) / 3)
    i_q_rms = power_devices_losses.I_RMS(i_p1, i_t1, d)

    # D
    i_d = (i_t + i_p) / 2 * (1 - d)
    i_d_rms = math.sqrt((1 - d) * (i_p**2 + i_p * i_t + i_t**2) / 3)

    # D3
    i_d3 = (i_t + i_p) / 2 * d
    i_d3_rms = math.sqrt(d * (i_p**2 + i_p * i_t + i_t**2) / 3)

    # Loss numbers
    # MOSFET
    R_q = 120e-3
    t_on = 50e-9
    t_off = 100e-9

    # diodes
    R_d = (1.7 - 0.7) / 20
    V_on_d = 0.7

    # Conduction losses
    cl_q = R_q * i_q_rms**2
    cl_q_2 = power_devices_losses.ConductionLoss.mosfet(R_q, i_q_rms)
    cl_d = i_d * V_on_d + R_d * i_d_rms**2
    cl_d_2 = power_devices_losses.ConductionLoss.diode(V_on_d, i_d, R_d,
                                                       i_d_rms)
    cl_d3 = i_d3 * V_on_d + R_d * i_d3_rms**2
    cl_d3_2 = power_devices_losses.ConductionLoss.diode(V_on_d, i_d3, R_d,
                                                        i_d3_rms)
    print("Conduction Loss:")
    print(f"Q: {cl_q}, {cl_q_2}")
    print(f"D: {cl_d}, {cl_d_2}")
    print(f"D3: {cl_d3}, {cl_d3_2}")

    # Switching losses
    # ON
    sl_on = 0.5 * i_t1 * t_on * V_s * f
    sl_on_2 = power_devices_losses.SwitchingLoss.power_turn_on(i_t1, V_s,
                                                               t_on, f)
    # OFF
    sl_off = 0.5 * i_p1 * t_off * V_s * f
    sl_off_2 = power_devices_losses.SwitchingLoss.power_turn_off(i_p1, V_s,
                                                                 t_off, f)
    print(f"On loss: {sl_on}, {sl_on_2}")
    print(f"Off loss: {sl_off}, {sl_off_2}")

    # Efficiency
    P = V_o * I_o
    P_loss = 2 * cl_q + cl_d + cl_d3 + sl_on + sl_off
    eff = P / (P + P_loss)
    print(f"Efficiency: {eff:%}")


if __name__ == '__main__':
    run()
