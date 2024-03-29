#!/usr/bin/env python
"""This modules provides functions to plot signals."""
import matplotlib.pyplot as plt
import numpy as np


# Signal Plotting
def plot_digital_as_digital(signal, modulation: str,
                            v_mode: bool = True) -> None:
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

    TODO
    ----
    Use Enums instead of numbers.
    """
    # Setting up ID for each type of modulation
    modulation_key = {
        "nrz unipolar": 0,
        "nrz bipolar": 1,
        "rz unipolar": 2,
        "rz bipolar": 3,
        "manchester": 4,
    }
    # Checking if the modulation given is valid or not
    if modulation.lower() not in modulation_key:
        raise ValueError("Invalid Modulation Type")

    # Saving modulation ID
    modulation_type = modulation_key[modulation.lower()]
    # Creating subplots
    _, ax = plt.subplots()

    # Finding the axis values
    # Only extracting 0 and 1 (IDK if I should raise and error or not)
    signal_array = []
    for i in signal:
        if i in "01":
            signal_array.append(int(i))
    # x values
    xs = np.arange(0, len(signal_array) + 1, 0.5)
    # Required for RZ bipolar
    ys = plot_ys(modulation_type, signal_array)

    # Formatting function
    def format_fn(tick_val, _):
        if int(tick_val) == 1:
            return "V+"
        if int(tick_val) == -1:
            return "V-"
        return '0'

    # Setting up axis
    ax.axis([0, len(xs) / 2 - 1, min(ys) - 0.1, max(ys) + 0.1])
    # Setting up axis ticks
    ax.yaxis.set_ticks([max(ys), 0, min(ys)])
    ax.xaxis.set_ticks(np.arange(0, len(xs) / 2, 1), [
        "" for i in range(int(len(xs) / 2))])
    # Setting up grid
    ax.grid(axis='x', color='b', linestyle='--', linewidth=1)

    # A FuncFormatter is created automatically.
    if v_mode:
        ax.yaxis.set_major_formatter(format_fn)
        ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.step(xs, ys)
    plt.show()


def nrz_unipolar(signal: list) -> None:
    """Return the y value of NRZ Unipolar."""
    if signal == 0:
        append_no = [0, 0]
    else:
        append_no = [1, 1]

    return append_no


def nrz_bipolar(signal: list) -> int:
    """Return the y value of NRZ Bipolar."""
    if signal == 0:
        append_no = [-1, -1]
    else:
        append_no = [1, 1]

    return append_no


def rz_unipolar(signal: list) -> int:
    """Return the y value of RZ Unipolar."""
    if signal == 0:
        append_no = [0, 0]
    else:
        append_no = [1, 0]

    return append_no


def rz_bipolar(signal: list, state: bool) -> tuple:
    """Return the y value of RZ Bipolar."""
    if signal == 0:
        append_no = [0, 0]
    else:
        if state:
            append_no = [1, 0]
        else:
            append_no = [-1, 0]
        state = not state

    return (append_no, state)


def manchester(signal: list) -> int:
    """Return the y value of Manchester."""
    if signal == 0:
        append_no = [0, 1]
    else:
        append_no = [1, 0]

    return append_no


def plot_ys(modulation_type: int, signal_array: list):
    """Plot the y values."""
    state = True
    # y values
    ys = [0]
    for i in signal_array:
        # Decides what to append
        if modulation_type == 0:
            append_no = nrz_unipolar(i)
        elif modulation_type == 1:
            append_no = nrz_bipolar(i)
        elif modulation_type == 2:
            append_no = rz_unipolar(i)
        elif modulation_type == 3:
            append_no, state = rz_bipolar(i, state)
        elif modulation_type == 4:
            append_no = manchester(i)

        # Appending numbers to y list
        for j in append_no:
            ys.append(j)
    ys.append(0)
    return ys


def plot_error_bar(axis: dict,
                   error: dict,
                   labels: dict = None,
                   show_graph: bool = True) -> None:
    """Plot a graph with an errorbar.
    Parameters
    ----------
    axis: dict[str, Iterable]
        A dictionary of the x and y values
        {
        "x": [1, 2, 3],
        "y": [1, 2, 3]
        }
    error: dict[str, Iterable]
        A dictionary of the x and y error values
        {
        "xerr": [1, 2, 3],
        "yerr": [1, 2, 3]
        }
    labels: dict[str, str]
        A dictionary of the labels of the graph
        Defaults to if None is given:
        {
        "xlabel": "x",
        "ylabel": "y"
        }
    showGraph
        To show the graph or not
    Returns
    -------
    None
    """
    if labels is None:
        labels = {
            "xlabel": "x",
            "ylabel": "y"
        }

    plt.xlabel(labels.get("xlabel"))
    plt.ylabel(labels.get("ylabel"))

    plt.errorbar(axis.get("x"), axis.get("y"), xerr=error.get("xerr"),
                 yerr=error.get("yerr"), capsize=2)

    if show_graph:
        plt.show()
