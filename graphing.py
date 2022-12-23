#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np


# Signal Plotting
def plot_digital_as_digital(signal, modulation: str,
                            vMode: bool = True) -> None:
    """Plot a digital signal that is transmitted is as a digital signal.

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
                if state:
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
    ax.axis([0, len(xs) / 2 - 1, min(ys) - 0.1, max(ys) + 0.1])
    # Setting up axis ticks
    # plt.yticks([1, 0, -1])
    ax.yaxis.set_ticks([max(ys), 0, min(ys)])
    ax.xaxis.set_ticks(np.arange(0, len(xs) / 2, 1), [
        "" for i in range(int(len(xs) / 2))])
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
    """Plot a graph with an errorbar
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
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)

    plt.errorbar(x, y, xerr=xerr, yerr=yerr, capsize=2)

    if showGraph:
        plt.show()
