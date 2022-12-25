#!/usr/bin/env python
"""Noise Example Sheet using eeepy."""
from eeepy.year_1.InfoandSystem import (SNR, NF_cascade, db_power,
                                        db_power_reverse,
                                        effective_noise_temperature_NF,
                                        thermal_noise)

if __name__ == '__main__':
    # Q1
    q1 = thermal_noise(300, 1e6)
    print(q1)

    # Q2
    q2_noise = thermal_noise(600, 10e6)
    q2 = SNR(5e-9, q2_noise)
    print(q2)

    # Q3
    q3_gain = [6, 12, 20]
    q3_gain_linear = [db_power_reverse(i, 1) for i in q3_gain]
    q3_NF = [1.7, 2, 4]
    q3_b = NF_cascade(q3_gain, q3_NF, NoiseFigureIn=True,
                      NoiseFigureReturn=True)
    q3_b = NF_cascade(q3_gain_linear, q3_NF, NoiseFigureIn=False,
                      NoiseFigureReturn=True)
    q3_linear = NF_cascade(q3_gain_linear, q3_NF, NoiseFigureIn=False,
                           NoiseFigureReturn=False)
    q3_b = db_power(1, q3_linear)
    print(q3_b)
    q3_c = effective_noise_temperature_NF(17, q3_linear, kelvin=False,
                                          NoiseFigure=False)
    q3_c = effective_noise_temperature_NF(17, q3_b, kelvin=False,
                                          NoiseFigure=True)
    print(q3_c)

    # Q4
    q4_A = effective_noise_temperature_NF(30, 1.87, False, False)
    q4_B = effective_noise_temperature_NF(13.5, 1.92, False, False)
    print(q4_A, q4_B)
