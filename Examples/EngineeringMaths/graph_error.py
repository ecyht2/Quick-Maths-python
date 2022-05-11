#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import argparse
from Maths import myStats
from InfoandSystem import plot_error_bar

# Arguments
def parse_arguments():
    # Intializing parser
    parser = argparse.ArgumentParser(description='Plots data located in a csv file with an option for error bar')

    # Adding Positional Arguments
    parser.add_argument('x_data', help="The csv file containing the x data", type=str)
    parser.add_argument('y_data', help="The csv file containing the y data", type=str)
    parser.add_argument('mode', help="How the data is being ploted", type=str, choices=["STD", "Range", "Max-Min", "None"], default="STD")
    parser.add_argument('x_label', help="The csv file containing the x data", type=str)
    parser.add_argument('y_label', help="The csv file containing the y data", type=str)

    # Adding Optional Arguments

    # Returning Parsed Arguments
    return parser.parse_args()

# Retriving values from parsed args
args = parse_arguments()
x_data = args.x_data
y_data = args.y_data
mode = args.mode
xLabel = args.x_label
yLabel = args.y_label

# Reading x values
xData = dict()
with open(x_data, "r") as f:
    line = f.readline()
    while not line == '':
        values = line.split(",")
        xData[values[0]] = {"Raw Data": []}
        for value in values[1:]:
            value = value.strip()
            try:
                xData[values[0]]["Raw Data"].append(float(value))
            except ValueError:
                pass
        line = f.readline()
# Reading y values
yData = dict()
with open(y_data, "r") as f:
    line = f.readline()
    while not line == '':
        values = line.split(",")
        yData[values[0]] = {"Raw Data": []}
        for value in values[1:]:
            value = value.strip()
            try:
                yData[values[0]]["Raw Data"].append(float(value))
            except ValueError:
                pass
        line = f.readline()

# Getting Statistics of Raw Data
for key in xData:
    xData[key]["Stats"] = myStats(xData[key]["Raw Data"])
for key in yData:
    yData[key]["Stats"] = myStats(yData[key]["Raw Data"])

# Getting x, xerr, y and yerr values according to mode
x = list()
xerr = list()
y = list()
yerr = list()
if mode == "STD":
    # Looping over xData
    for data in xData.values():
        stats = data["Stats"]
        x.append(stats.mean)
        xerr.append(stats.std)
    # Looping over yData
    for data in yData.values():
        stats = data["Stats"]
        y.append(stats.mean)
        yerr.append(stats.std)
elif mode == "Range":
    # Looping over xData
    for data in xData.values():
        stats = data["Stats"]
        x.append(stats.median)
        xerr.append(stats.range)
    # Looping over yData
    for data in yData.values():
        stats = data["Stats"]
        y.append(stats.median)
        yerr.append(stats.range)
elif mode == "Max-Min":
    # Looping over xData
    for data in xData.values():
        stats = data["Stats"]
        x.append(stats.median)
        xerr.append([x[-1] - stats.min, stats.max - x[-1]])
    # Looping over yData
    for data in yData.values():
        stats = data["Stats"]
        y.append(stats.median)
        yerr.append([y[-1] - stats.min, stats.max - y[-1]])

    # Transposing errors
    # x
    xerr = np.array(xerr)
    xerr = xerr.transpose()
    # y
    yerr = np.array(yerr)
    yerr = yerr.transpose()
elif mode == "None":
    # Looping over xData
    for data in xData.values():
        stats = data["Stats"]
        x.append(stats.mean)
        xerr.append(0)
    # Looping over yData
    for data in yData.values():
        stats = data["Stats"]
        y.append(stats.mean)
        yerr.append(0)

print(x)
print(y)
# Plotting data
plot_error_bar(x, y, xerr, yerr, xLabel, yLabel)
