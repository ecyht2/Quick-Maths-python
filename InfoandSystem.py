#!/usr/bin/env python
from string import ascii_uppercase
from Constants import *
from Maths import db_power, db_volts, db_power_reverse
from Maths import exp, product, log1p

# Diodes
def bolzmann_diode_equation(IS: float, VD: float, VT: float = VT) -> float:
    """
    Calculates the current in a pn junction diode
    """
    ID = IS * (exp(VD/VT) - 1)
    return ID
def bolzmann_diode_equation_rev(ID: float, IS: float, VT: float = VT) -> float:
    """
    Calculates the voltage in a pn junction diode
    """
    VD = VT * (log1p(ID/IS + 1))
    return VD

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

def convolution(x: list, h: list) -> list:
    """
    Find the convolution of data stream x given impulse response h

    Parameters
    ----------
    x
        A list containing all the data of the data stream
    h
        A list containing all the impulse response

    Returns
    -------
    list
        The convoluted result
    """
    # Can be achieved via numpy.convolve(x, h)
    result = []
    row = []
    xSize = len(x)
    hSize = len(h)

    # Finding y per sample
    for rows in range(xSize):
        # Appending leading 0
        column = [0 for i in range(rows)]
        # Finding output of x[n]
        for i in h:
            column.append(i * x[rows])
        # Appending trailing 0
        while len(column) < xSize + hSize - 1:
            column.append(0)
        # Appending row to column
        row.append(column)

    # Summing y per sample
    for columns in range(xSize + hSize - 1):
        # Summing the columns
        total = 0
        for rows in range(xSize):
            total += row[rows][columns]
        # Appending to rsult
        result.append(total)

    # Returning result
    return result

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
