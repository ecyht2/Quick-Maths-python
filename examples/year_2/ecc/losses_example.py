#!/usr/bin/env python
"""Example script that uses the power_devices_losses module.

The questions are from the power devices losses example question.
"""
import math
from eeepy.year_2.ecc import power_devices_losses


def run():
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
    i = {}
    # i_q = (i_t1 + i_p1) / 2 * d
    i["q"] = {}
    i["q"]["rms"] = math.sqrt(d * (i_p1**2 + i_p1 * i_t1 + i_t1**2) / 3)
    i["q"]["rms"] = power_devices_losses.I_RMS(i_p1, i_t1, d)

    # D
    i["d"] = {}
    i["d"]["avg"] = (i_t + i_p) / 2 * (1 - d)
    i["d"]["rms"] = math.sqrt((1 - d) * (i_p**2 + i_p * i_t + i_t**2) / 3)

    # D3
    i["d3"] = {}
    i["d3"]["avg"] = (i_t + i_p) / 2 * d
    i["d3"]["rms"] = math.sqrt(d * (i_p**2 + i_p * i_t + i_t**2) / 3)

    # Loss numbers
    # MOSFET
    mosfet = {}
    mosfet["R_q"] = 120e-3
    mosfet["t_on"] = 50e-9
    mosfet["t_off"] = 100e-9

    # diodes
    diodes = {}
    diodes["R_d"] = (1.7 - 0.7) / 20
    diodes["V_on_d"] = 0.7

    # Conduction losses
    conduction_loss = {}
    conduction_loss["cl_q"] = mosfet["R_q"] * i["q"]["rms"]**2
    conduction_loss["cl_q_2"] = power_devices_losses.ConductionLoss.mosfet(
        mosfet["R_q"], i["q"]["rms"]
    )
    conduction_loss["cl_d"] = i["d"]["avg"] * diodes["V_on_d"] +\
        diodes["R_d"] * i["d"]["rms"]**2
    conduction_loss["cl_d_2"] = power_devices_losses.ConductionLoss.diode(
        diodes["V_on_d"], i["d"]["avg"], diodes["R_d"], i["d"]["rms"]
    )
    conduction_loss["cl_d3"] = i["d3"]["avg"] * diodes["V_on_d"] +\
        diodes["R_d"] * i["d3"]["rms"]**2
    conduction_loss["cl_d3_2"] = power_devices_losses.ConductionLoss.diode(
        diodes["V_on_d"], i["d3"]["avg"], diodes["R_d"], i["d3"]["rms"]
    )
    print("Conduction Loss:")
    print(f"Q: {conduction_loss['cl_q']}, {conduction_loss['cl_q_2']}")
    print(f"D: {conduction_loss['cl_d']}, {conduction_loss['cl_d_2']}")
    print(f"D3: {conduction_loss['cl_d3']}, {conduction_loss['cl_d3_2']}")

    # Switching losses
    switching_loss = {}
    # ON
    switching_loss["sl_on"] = 0.5 * i_t1 * mosfet["t_on"] * V_s * f
    switching_loss["sl_on_2"] =\
        power_devices_losses.SwitchingLoss.power_turn_on(
            i_t1, V_s, mosfet["t_on"], f
    )
    # OFF
    switching_loss["sl_off"] = 0.5 * i_p1 * mosfet["t_off"] * V_s * f
    switching_loss["sl_off_2"] =\
        power_devices_losses.SwitchingLoss.power_turn_off(
            i_p1, V_s, mosfet["t_off"], f
    )
    print(f"On loss: {switching_loss['sl_on']}, {switching_loss['sl_on_2']}")
    print(f"Off loss: {switching_loss['sl_off']}, "
          f"{switching_loss['sl_off_2']}")

    # Efficiency
    power = {}
    power["in"] = V_o * I_o
    power["loss"] = 2 * conduction_loss["cl_q"] + conduction_loss["cl_d"] +\
        conduction_loss["cl_d3"] + switching_loss["sl_on"] +\
        switching_loss["sl_off"]
    power["eff"] = power["in"] / (power["in"] + power["loss"])
    print(f"Efficiency: {power['eff']:%}")


if __name__ == '__main__':
    run()
