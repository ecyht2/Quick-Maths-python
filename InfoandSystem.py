#!/usr/bin/env python
from string import ascii_uppercase

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
