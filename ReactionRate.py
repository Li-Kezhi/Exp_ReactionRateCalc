#!/usr/bin/env python

"""
Reaction rate post-treatment script
"""

__author__ = "LI Kezhi"
__date__ = "$2017-03-05$"
__version__ = "1.0.1"


import math
import numpy as np
import matplotlib.pyplot as plt

import mpltex


##### Preparation #####
# File names
SOURCE_NAME = './Examples/Industry1.txt'

# Experiment details
NAME = ('NO', 'NH3', 'N2O', 'NO2', 'SO2', 'H2O') # The conversion ratio of 1st gas is calculated!
BACKGROUND = dict(zip(NAME, (500, 500, 0, 0, 0, 0)))
FLOAT_RATE = 100  # ml/min
VOLUMN = 0.1  # ml
SANNING_SPEED = 0.3521689  # min per line
TEMPERATURE = (100, 130, 170, 200,
               230, 260, 300, 350,
               400, 450, 500) # Celsius degree
TIME = (60, 110, 160, 210,
        240, 270, 300, 330,
        360, 390, 420) # min

# File formatting
DATA_ROW = (10, 19, 7, 13, 16, 4) # Correspond to NAME sequence

# Read data
data = np.loadtxt(SOURCE_NAME, skiprows=1, usecols=DATA_ROW)

x = np.zeros_like(data)
for i in xrange(len(x)):
    x[i, 0] = i

##### Plotting #####
@mpltex.presentation_decorator
def plot(data, SOURCE_NAME):
    fig, ax = plt.subplots()
    for index, item in enumerate(NAME):
        ax.plot(x[:, 0], data[:, index], label=item)

    # ax.set_xlim(5, 80)
    ax.set_ylim(-0.1, 550)

    # ax.set_yticks([]) # Hide yticks
    # ax.tick_params(axis='x', top='off', bottom='off')

    ax.legend(loc='best')
    ax.set_ylabel(r'Concentration (ppm)')

    plt.show()

##### Reports #####
# Data manipulation
ave, std = {}, {} # Average concentration and standard error
for item in NAME:
    ave[item] = np.zeros_like(TEMPERATURE, dtype=float)
    std[item] = np.zeros_like(TEMPERATURE, dtype=float)

aveConv = np.zeros_like(TEMPERATURE, dtype=float)
stdConv = np.zeros_like(TEMPERATURE, dtype=float)
ave_k = np.zeros_like(TEMPERATURE, dtype=float)
std_k = np.zeros_like(TEMPERATURE, dtype=float)
ave_lnk = np.zeros_like(TEMPERATURE, dtype=float)
std_lnk = np.zeros_like(TEMPERATURE, dtype=float)
T_inv = np.zeros_like(TEMPERATURE, dtype=float)
linePosition = []
for i in xrange(len(TEMPERATURE)):
    time_i = TIME[i]
    linePosition.append(int(time_i/SANNING_SPEED))

    for index, item in enumerate(NAME):
        ave[item][i] = np.average(data[linePosition[i]-25:linePosition[i],
                                       index])
        std[item][i] = abs(np.std(data[linePosition[i]-25:linePosition[i], 0],
                                  ddof=index))

    # Conversion ratio
    aveConv[i] = (1 - ave[NAME[0]][i] / BACKGROUND[NAME[0]]) * 100  # X = 1 - C / C0
    stdConv[i] = abs(std[NAME[0]][i] / BACKGROUND[NAME[0]] * 100)
    try:
        assert aveConv[i] < 100
    except AssertionError:
        aveConv[i] = 99.99999
    # Reaction rate constant
    try:
        ave_k[i] = -FLOAT_RATE / VOLUMN * math.log(1 - aveConv[i]/100) / 60
          # k = - F/V * ln(1 - X)
    except ValueError:
        ave_k[i] = np.NaN
    std_k[i] = abs(-FLOAT_RATE / VOLUMN / (1 - aveConv[i]/100) * stdConv[i]/100) / 60
      # delta k = - F/V * 1/(1 - X) * delta X
    # lnk - 1/T
    try:
        ave_lnk[i] = math.log(ave_k[i])
        std_lnk[i] = abs(stdConv[i]/100 / ((1-stdConv[i]/100) * math.log(1-stdConv[i]/100)))
    except ValueError:
        ave_lnk[i], std_lnk[i] = np.NaN, np.NaN
    T_inv[i] = 1 / (273.15 + TEMPERATURE[i])

# Generate Report
reportName = SOURCE_NAME.split('.txt')[0] + '_report.txt'
output = file(reportName, 'w')
headLine = 'Temperature(C)  '
for item in NAME:
    headLine += item
    headLine += '(ppm)  StandardError(ppm)  '
headLine += ('ConvertionRatio(%)  StandardError(%)  ' +
             'ReactionRate(s-1)  StandardError(s-1)  ' +
             'lnk  StandardError  ' +
             '1/T\n')
output.write(headLine)
for i in xrange(len(TEMPERATURE)):
    dataLine = '%3d   ' % TEMPERATURE[i]
    for item in NAME:
        dataLine += '%9.4f   %9.4f   ' % (ave[item][i], std[item][i])
    dataLine += '%6.4f   %6.4f   ' % (aveConv[i], stdConv[i])
    dataLine += '%6.4f   %6.4f   ' % (ave_k[i], std_k[i])
    dataLine += '%6.4f   %6.4f   ' % (ave_lnk[i], std_lnk[i])
    dataLine += '%9.7f\n' % T_inv[i]
    output.write(dataLine)
output.close()

plot(data, SOURCE_NAME)
