#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from string import ascii_uppercase
import numpy as np
from string import digits
from Constants import *
from Maths import db_power, db_volts, db_power_reverse, product

def kMap(size, equation):
    pass

def letter_to_binary(equation: str) -> str:
    """
    Convert equations writen in letters into it's binary form

    Parameters
    ----------
    equation
        The equation to convert to binary

    Returns
    -------
    str
        The equation converted to binary
    """
    eq = ""
    for i in range(len(equation)):
        letter = equation[i]
        if letter in ascii_uppercase:
            try:
                if(equation[i+1] == '\''):
                    letter = "0"
                else:
                    letter = "1"
            except IndexError:
                letter = "1"
        else:
            if letter == '\'':
                letter = ""
        eq = "".join([eq, letter])
    return eq

convolution = np.convolve

def filter(b: list, a: list, x: list) -> list:
    """
    Filters the input data x using "FIR or IIR"

    Parameters
    ----------
    b
        Numerator coefficient
    a
        Denominator coefficient
    x
        An array that contains the intput data

    Returns
    -------
    list
        The output of the filter
    """
    pass

def moving_average_filter(x: list, n: int) -> list:
    """
    Find the mean average filter for data x

    Parameters
    ----------
    x
        An array that contains the input data
    n
        The number of points to average

    Returns
    -------
    list
        The result of the moving_average_filter
    """
    # Getting the impulse response needed for convolution
    h = []
    for i in range(n):
        h.append(1/n)

    # Getting the moving average result array
    result = convolution(x, h)
    while len(result) > len(x):
        result.pop()

    # Returning result
    return result

# Signal Plotting
def plot_digital_as_digital(signal, modulation: str, vMode: bool = True) -> None:
    """
    Plot a digital signal that is transmitted is as a digital signal

    Parameters
    ----------
    signal
        The signal being transmitted
        It can be in a str or an array_like
        Only 0 and 1 will be considered
    modulation
        The type of modulation used

    Returns
    -------
    None
    """
    # Setting up ID for each type of modulation
    modulationKey = {
        "nrz unipolar": 0,
        "nrz bipolar": 1,
        "rz unipolar": 2,
        "rz bipolar": 3,
        "manchester": 4,
    }
    # Checking if the modulation given is valid or not
    if modulation.lower() not in modulationKey:
        raise ValueError("Invalid Modulation Type")

    # Saving modulation ID
    modulationType = modulationKey[modulation.lower()]
    # Creating subplots
    fig, ax = plt.subplots()

    # Finding the axis values
    # Only extracting 0 and 1 (IDK if I should raise and error or not)
    signal_array = []
    for i in signal:
        if i in "01":
            signal_array.append(int(i))
    # x values
    xs = np.arange(0, len(signal_array) + 1, 0.5)
    # Required for RZ bipolar
    state = True
    # y values
    ys = [0]
    for i in signal_array:
        # Decides what to append
        if modulationType == 0:
            if i == 0:
                append_no = [0, 0]
            else:
                append_no = [1, 1]
        elif modulationType == 1:
            if i == 0:
                append_no = [-1, -1]
            else:
                append_no = [1, 1]
        elif modulationType == 2:
            if i == 0:
                append_no = [0, 0]
            else:
                append_no = [1, 0]
        elif modulationType == 3:
            if i == 0:
                append_no = [0, 0]
            else:
                if state == True:
                    append_no = [1, 0]
                else:
                    append_no = [-1, 0]
                state = not state
        elif modulationType == 4:
            if i == 0:
                append_no = [0, 1]
            else:
                append_no = [1, 0]

        # Appending numbers to y list
        for j in append_no:
            ys.append(j)
    ys.append(0)

    # Formatting function
    def format_fn(tick_val, tick_pos):
        if int(tick_val) == 1:
            return "V+"
        elif int(tick_val) == -1:
            return "V-"
        else:
            return '0'

    # Setting up axis
    # plt.axis([0, 9, -1.1, 1.1])
    ax.axis([0, len(xs)/2 - 1, min(ys) - 0.1, max(ys) + 0.1])
    # Setting up axis ticks
    # plt.yticks([1, 0, -1])
    ax.yaxis.set_ticks([max(ys), 0, min(ys)])
    ax.xaxis.set_ticks(np.arange(0, len(xs)/2, 1), ["" for i in range(int(len(xs)/2))])
    # Setting up grid
    # plt.grid(axis='x', color='b', linestyle='--')
    ax.grid(axis='x', color='b', linestyle='--', linewidth=1)

    # A FuncFormatter is created automatically.
    if vMode:
        ax.yaxis.set_major_formatter(format_fn)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.step(xs, ys)
    plt.show()

# Noise
def thermal_noise(T, B, k = kBolzman):
    return k*T*B

def C_to_kelvin(T):
    return T + 273
def kelvin_to_C(T):
    return T - 273

# SNR
def SNR(S, N, dB = True, power=True):
    returnValue = 0
    if dB:
        if power:
            returnValue = db_power(N, S)
        else:
            returnValue = db_volts(N, S)
    else:
        returnValue = S/N

    return returnValue

# NF
def NF(inputSNR, outputSNR, NoiseFigure: bool = True):
    returnValue = 0
    if NoiseFigure:
        returnValue = db_power(inputSNR, outputSNR)
    else:
        returnValue = inputSNR/outputSNR

    return returnValue
def noise_factor_internal(internal, inputNoise):
    return 1 + internal/inputNoise
def NF_cascade(gain: iter, NF: iter, NoiseFigureIn: bool = True, NoiseFigureReturn: bool = True):
    # Checking for conditions
    if not len(gain) == len(NF):
        raise ValueError("Size of gain and NF aren't equal")

    # Changing to noise figure into noise factor
    if NoiseFigureIn:
        for i in range(len(gain)):
            gain[i] = db_power_reverse(gain[i], 1)
    # Using absolute values
    for i in range(len(gain)):
        gain[i] = abs(gain[i])
        NF[i] = abs(NF[i])

    # Calculating NF
    totalNF = NF[0]
    for i in range(len(NF) - 1):
        totalNF += (NF[i+1] - 1) / product(gain[:i+1])

    returnNF = 0
    if NoiseFigureReturn:
        returnNF = db_power(1, totalNF)
    else:
        returnNF = totalNF

    return returnNF

def effective_noise_temperature_NF(T, NF, kelvin: bool = True, NoiseFigure: bool = True):
    if not kelvin:
        T = C_to_kelvin(T)
    if NoiseFigure:
        NF = db_power_reverse(NF, 1)

    return T * (NF - 1)
def effective_noise_temperature_gain(T, gain, kelvin: bool = True):
    if not kelvin:
        for i in range(len(T)):
            T[i] = C_to_kelvin(T[i])

    totalT = T[0]
    for i in range(len(T) - 1):
        totalT += T[i+1] / product(gain[:i+1])

    return totalT
