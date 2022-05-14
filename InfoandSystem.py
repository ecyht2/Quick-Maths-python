#!/usr/bin/env python
from string import ascii_uppercase
from string import digits
from Constants import *
from Maths import db_power, db_volts, db_power_reverse
from Maths import product, log1p, log2
from Maths import exp, factorial
import csv

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

# Angle Modulation
def FM_PM_deviation_sensitivity(delta: float, em: float) -> float:
    """
    Calculates the deviation sensitivity of the modulator of a FM or PM signal (kf/kp)
    """
    return delta / em
def FM_PM_modulating_index(delta: float, n: float) -> float:
    """
    Calculates the modulating index of a FM or PM signal
    """
    return delta / n
def bessel_function(x: float, v: int) -> float:
    """
    Calculates the J value of a given order v with the β value x
    Jv(x)
    """
    if v < 0:
        raise ValueError("v must be a positive integer number")
    # Initializing variables
    sValues:list = list()
    ds: float = 10
    s: int = 0
    while ds > 0.0001:
        # https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470054208.app3
        sValues.append((x / 2)**(2 * s + v)*(-1)**s / (factorial(s) * factorial(s+v)))
        # Updating change in the current value
        if len(sValues) > 1:
            ds = abs(sValues[-1] - sValues[-2])
        # Going to next s
        s += 1
    # Summing up the values to get J
    J = sum(sValues)
    return J

def FM_PM_J_values(m: float, cached: bool = True) -> tuple[float]:
    """
    Find the J values of a given modulating index
    """
    J = list()
    # Getting items from cached
    if m <= 5 and cached:
        m = round(m, 1)
        with open('bessel.csv', 'r') as csvfile:
            bessel = csv.reader(csvfile)
            for i in bessel:
                for j in range(len(i)):
                    try:
                        i[j] = float(i[j])
                    except ValueError:
                        pass
                if i[0] == m:
                    J = i[1:]
    # Calculating values manually
    else:
        Jvalue: float = 10
        # Getting the first 14 J values
        for v in range(15):
            Jvalue = bessel_function(m, v)
            J.append(round(Jvalue, 4))
        # Removing trailing values that are < 0.01
        while abs(J[-1]) < 0.01:
            J.pop()
    return tuple(J)
def FM_PM_bandwidth_bessel(fm: float, n: float = 0, m: float = 0) -> float:
    """
    Find the bandwith of a given FM or PM signal using Bessel's frequency spectrum
    """
    bandwidth: float = 0
    if n > 0:
        bandwidth = 2*(n * fm)
    elif m > 0:
        n = len(FM_PM_J_values(m)) - 1
        bandwidth = 2*(n * fm)
    else:
        raise ValueError("No n or m given")
    return bandwidth
def FM_PM_bandwidth_carson(fm: float, delta: float) -> float:
    """
    Find the bandwith of a given FM or PM signal using Carson's rule
    """
    return 2 * (fm + delta)
def FM_PM_power_transmitted(R: float, Vcrms: float = 0,
                            J: tuple[float] or list[float] = tuple()) -> float:
    """
    Find the transmitted power of a FM or PM signal
    """
    Pt: float = 0
    if Vcrms > 0:
        Pt = Vcrms**2 / R
    elif len(J) > 0:
        totalJ = J[0]**2
        for i in J[1:]:
            totalJ += 2 * i**2
        Pt = totalJ / R
    else:
        raise ValueError("No Vcrms or J given")
    return Pt

# PM Pulse Modulation
def mu_law(mu: int, Vin: float, Vmax: float) -> float:
    """
    Compands the signal using μ-law (American)
    """
    Vout = Vmax * log1p(1 + mu * Vin / Vmax) / log1p(1 + mu)
    return Vout
def shanon_hartleys_formula(BW: float, M: int) -> float:
    """
    Calculates the maximum data that can be sent in a given bandwidth
    """
    return 2 * BW * log2(M)
def shanon_limit(BW: float, SNR: float) -> float:
    """
    Calculates the maximum data that can be sent in a given bandwidth
    """
    return BW * log2(1 + SNR)

# Fourier Analysis
def dft_number_of_samples(T: float, ts: float) -> float:
    """
    Calculates the amount of samples
    """
    return T / ts
def dft_sampling_period(T: float, N: float) -> float:
    """
    Calculates the sampling period
    """
    return T / N
def dft_sample_duration(N: float, ts: float) -> float:
    """
    Calculates the duration of the sampled data
    """
    return N * ts
def dft_max_frequency(fs: float) -> float:
    """
    Calculates the max frequency that can be sampled
    """
    return fs /2
def dft_min_frequency(T: float) -> float:
    """
    Calculates the min frequency that can be sampled
    """
    return 1 / T

# Analogue Filter
def low_pass_filter_voltage_gain(Af: float, f: float, fc: float) -> float:
    """
    Calculates the voltage gain of a low pass filter
    """
    Av = Af / (1 + (f / fc)**2)**0.5
    return Av
def low_pass_filter_cutoff_frequency(R: float, C: float) -> float:
    """
    Calculates the cut-off frequency of a low pass filter
    """
    return 1 / (2 * pi * R * C)
def low_pass_filter_capacitor(R: float, fc: float) -> float:
    """
    Calculates the capacitor needed of a low pass filter
    """
    return 1 / (2 * pi * R * fc)
def second_order_low_pass_filter_cutoff_frequency(R: list[float], C: list[float]) -> float:
    """
    Calculates the cut-off frequency of a second order low pass filter
    """
    return 1 / (2 * pi * (product(R) * product(C))**0.5)

# Digital Filter
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
def thermal_noise(T: float, B: float, k: float = kBolzman) -> float:
    return k*T*B

def C_to_kelvin(T: float) -> float:
    return T + 273
def kelvin_to_C(T: float) -> float:
    return T - 273

# SNR
def SNR(S: float, N: float, dB: bool = True, power: bool = True) -> float:
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
def NF(inputSNR: float = 0, outputSNR: float = 0,
       Nai: float = 0, Ni: float = 0,
       NoiseFigure: bool = True) -> float:
    returnValue = 0
    if inputSNR > 0 and outputSNR > 0:
        returnValue = inputSNR/outputSNR
    elif Nai > 0 and Ni > 0:
        returnValue = 1 + Nai / Ni
    else:
        raise ValueError("No (inputSNR and outputSNR) or (Nai and Ni) given")
    if NoiseFigure:
        returnValue = db_power(1, returnValue)

    return returnValue
def NF_cascade(gain: list or tuple, NF: list or tuple,
               NoiseFigureIn: bool = True, NoiseFigureReturn: bool = True) -> float:
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

    # Sorting values
    gain.sort()
    NF.sort()

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

def effective_noise_temperature_NF(T: float, NF: float,
                                   kelvin: bool = True, NoiseFigure: bool = True) -> float:
    if not kelvin:
        T = C_to_kelvin(T)
    if NoiseFigure:
        NF = db_power_reverse(NF, 1)

    return T * (NF - 1)
def effective_noise_temperature_gain(T: list or tuple, gain: list or tuple,
                                     kelvin: bool = True) -> float:
    if not kelvin:
        for i in range(len(T)):
            T[i] = C_to_kelvin(T[i])

    totalT = T[0]
    for i in range(len(T) - 1):
        totalT += T[i+1] / product(gain[:i+1])

    return totalT
