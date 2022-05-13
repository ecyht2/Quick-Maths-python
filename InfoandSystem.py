#!/usr/bin/env python
from string import ascii_uppercase
from string import digits
from Constants import *
from Maths import db_power, db_volts, db_power_reverse
from Maths import exp, product, log1p, log2
import numpy as np

# Diodes
def bolzmann_diode_equation(IS: float, VD: float, VT: float = VT) -> float:
    """
    Calculates the current in a pn junction diode
    """
    ID = IS * (exp(VD/VT) - 1)
    return ID
def bolzmann_diode_equation_reverse(ID: float, IS: float, VT: float = VT) -> float:
    """
    Calculates the voltage in a pn junction diode
    """
    VD = VT * (log1p(ID/IS + 1))
    return VD

# Op-Amp
def inverting_op_amp(R1: float, R2: float) -> dict:
    """
    Calculates the gain and input resistance of an Inverting Op-Amp
    """
    A = - R2/R1
    R = R1
    return {"Gain": A, "Input Resistance": R}
def non_inverting_op_amp(R1: float, R2: float) -> dict:
    """
    Calculates the gain and input resistance of a Non-Inverting Op-Amp
    """
    A = 1 + R2/R1
    R = "Infinity"
    return {"Gain": A, "Input Resistance": R}
def transimpedence_op_amp(Iin: float, Rf: float) -> float:
    """
    Calculates the output voltage of a transimpedence Op-Amp given an input current
    """
    return Iin * Rf

# Transistor
def transistor_beta(IC: float, IB: float) -> float:
    """
    Calculates the β (common-emitter current gain) of a transistor
    Should be between 50 and 200
    """
    beta = IC/IB
    return beta
def transistor_alpha(beta: float) -> float:
    """
    Calculates the α common-base current gain
    Should be slightly less than 1
    """
    alpha = beta / (1 + beta)
    return alpha
def transistor_alpha_IE(IC: float, IE: float) -> float:
    """
    Calculates the α common-base current gain using IC and IE
    Should be slightly less than 1
    """
    alpha = IC/IE
    return alpha
def transistor_mode(VE: float, VB: float, VC: float, transistorType: str = "NPN") -> str:
    """
    Determines the mode the transistor is operating in given VE, VB and VC
    """
    # Defining Modes
    modes = ["Active", "Saturation", "Cutoff", "Reverse"]

    # Logic
    if VC > VB and VB > VE:
        mode = 0
    elif VB > VE and VB > VC:
        mode = 1
    elif VE > VB and VC > VB:
        mode = 2
    elif VE > VB and VB > VC:
        mode = 3

    if transistorType.upper() == "PNP":
        mode = -(mode + 1)
    elif transistorType.upper() == "NPN":
        mode = mode
    else:
        raise ValueError("Invalid BJT type")
    return modes[mode]
def drain_current_mosfet(K: float, VGS: float, VThresh: float, mosfetType: str) -> float:
    """
    Calculates the drain current (ID) of a mosfet
    """
    ID = 0
    if mosfetType.upper() == "NMOS":
        ID = K * (VGS - VThresh)**2
    elif mosfetType.upper() == "PMOS":
        ID = K * (VGS + VThresh)
    else:
        raise ValueError("Invalid MOSFET type")
    return ID

# Number System
# Commit # 0110 1001  Nice
# BCD
def decimal_to_bcd(number: float) -> str:
    """
    Converts Decimal Number to BCD
    """
    if not (type(number) == int or type(number) == float):
        raise TypeError("Input number must be an integer")

    numberString = str(number)
    BCD = ""
    for place in numberString:
        if place == ".":
            BCD += ". "
            continue
        placeInt = int(place)
        BCD += bin(placeInt)
        BCD += " "

    BCD = BCD.rstrip()
    BCD = BCD.replace("0b", "")
    splitedBCD = BCD.split(" ")
    for i in range(len(splitedBCD)):
        if splitedBCD[i] == ".":
            continue
        while len(splitedBCD[i]) < 4:
            splitedBCD[i] = "0" + splitedBCD[i]
    BCD = " ".join(splitedBCD)

    return BCD
def bcd_to_decimal(number: str) -> float:
    """
    Converts BCD Number to Decimal
    """
    splittedNumber: list[str] = number.split(" ")
    number = ""
    for digit in splittedNumber:
        if digit == ".":
            number += "."
            continue
        digitInt = int(digit, base = 2)
        number += str(digitInt)
    return float(number)
# Gray Code
def binary_to_gray(number: str) -> str:
    """
    Converts binary number to gray code
    """
    grayCode = "1"
    for i in range(len(number) - 1):
        binary = int(number[i + 1])
        binaryPrev = int(number[i])
        grayCode += str(int(bool(binary) ^ bool(binaryPrev)))
    return grayCode
def gray_to_binary(number: str) -> str:
    """
    Converts gray code number to binary
    """
    binaryNumber = "1"
    for i in range(len(number) - 1):
        binary = int(number[i + 1])
        binaryPrev = int(binaryNumber[i])
        binaryNumber += str(int(bool(binary) ^ bool(binaryPrev)))
    return binaryNumber

# Combinational Logic Circuit
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

# A2D or ADC
def nyquist_shanon_sampling_frequency(fmax: float = 0, B: float = 0) -> float:
    """
    Find the minimum sampling frequency required to digitize an analogue signal
    """
    fs = 0
    if fmax > 0:
        fs = 2 * fmax
    elif B > 0:
        fs = 2 * B
    return fs
def quantization_level(b: int) -> int:
    """
    Find the number of quantization level for a b amount of bits ADC
    """
    return 2**b
def quantization_level_reverse(m: int) -> int:
    """
    Find the number of bits an ADC has given it has m number of quantization level
    """
    return log2(m)

# BW
def bandwidth(fmax: float, fmin: float) -> float:
    """
    Calculates the bandwidth of a signal
    """
    return fmax - fmin

# Information
def information_content(N: int, b: int = 0, P: float = 0) -> float:
    """
    Calculates the information content (H) of a digitalize signal
    """
    H = 0
    if P == 0:
        H = N * b
    elif b == 0:
        H = (N * P * -log2(P))
    else:
        raise ValueError("No b or P is given")
    return H

# AM
def AM_modulating_index(Vm: float = 0, Vc: float = 0,
                        Vmax: float = 0, Vmin: float = 0,
                        Pt: float = 0, Pc: float = 0) -> float:
    """
    Calculates the modulating index of an AM signal
    """
    m = 0
    if Vm > 0 and Vc > 0:
        m = Vm / Vc
    elif Vmax > 0 and Vmin > 0:
        m = (Vmax - Vmin) / (Vmax + Vmin)
    elif Pt > 0 and Pc > 0:
        m = 2 * ((Pt / Pc) - 1)
        m = m**0.5
    else:
        raise ValueError("No (Vm and Vc) or (Vmax and Vmin) or (Pt and Pc) given")
    return m
def AM_Vm(Vmax: float, Vmin: float) -> float:
    """
    Calculates the Voltage of the modulating signal (Vm)
    """
    return (Vmax - Vmin) / 2
def AM_Vc(Vmax: float, Vmin: float) -> float:
    """
    Calculates the Voltage of the carrier signal (Vc)
    """
    return (Vmax + Vmin) / 2
def AM_BW(fm: float) -> float:
    """
    Calculates the bandwidth of an AM signal
    """
    return 2 * fm
def AM_sidebands(fc: float, fm: float) -> tuple[float, float]:
    """
    Calculates the sidebands of an AM signal
    """
    return (fc - fm, fc + fm)
def AM_power_transmitted(Pc: float, m: float) -> float:
    """
    Calculates the power of the transmitted AM signal
    """
    return Pc * (1 + m**2 / 2)
def AM_power_carrier(Pt: float, m: float) -> float:
    """
    Calculates the power of the carrier of an AM signal
    """
    return Pt / (1 + m**2 / 2)
def AM_modulating_index_sum(m: list or tuple, *argc: tuple[float]) -> float:
    """
    Calculates the total modulating index of mutiple AM singal simultaneously
    """
    mType = type(m)
    mt: float = 0
    if mType == list or mType == tuple:
        for i in m:
            mt += i**2
    elif mType == int or mType == float:
        mt += m**2
        for i in argc:
            mt += i**2

    return mt**0.5
def AM_modulating_index_sum_voltage(Vc: float, Vm: list or tuple, *argc: tuple[float]) -> float:
    """
    Calculates the total modulating index of mutiple AM singal simultaneously given the voltages
    """
    VType = type(Vm)
    mt: float = 0
    if VType == list or VType == tuple:
        for i in Vm:
            mt += i**2
    elif VType == int or VType == float:
        mt += Vm**2
        for i in argc:
            mt += i**2

    return mt**0.5/Vc

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

def plot_error_bar(x: iter, y: iter,
                   xerr: iter, yerr: iter,
                   xLabel: str = "x", yLabel: str = "y",
                   showGraph: bool = True) -> None:
    """
    Plot a graph with an errorbar

    Parameters
    ----------
    x
        An array_like containing all the x values
    y
        An array_like containing all the y values
    xerr
        An array_like containing all error of the x values
    yerr
        An array_like containing all error of the y values
    xLabel
        The label for x-axis
    yLabel
        The label for y-axis
    showGraph
        To show the graph or not

    Returns
    -------
    None
    """
    fig = plt.figure()
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)

    plt.errorbar(x, y, xerr=xerr, yerr=yerr, capsize=2)

    if showGraph:
        plt.show()

# Digital Filter
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

# Noise
def thermal_noise(T: float, B: float, k: float = kBolzman) -> float:
    return k*T*B

def C_to_kelvin(T: float):
    return T + 273
def kelvin_to_C(T: float):
    return T - 273

# SNR
def SNR(S: float, N: float, dB: bool = True, power: bool = True):
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
