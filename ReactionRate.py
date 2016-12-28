#!/usr/bin/env python

"""
Reaction rate post-treatment script
"""

__author__ = "LI Kezhi" 
__date__ = "$2016-12-23$"
__version__ = "1.0.0"


import math
import numpy as np
import matplotlib.pyplot as plt

import mpltex


##### Preparation #####
# File names
SOURCE_NAME = './Examples/beta.txt'

# Experiment details
BACKGROUND = {
    'bgCO': 1.045,
    'bgCO2': 0.390695225,
    'bgH2O': 0.02857845
}  # Background concentration
FLOAT_RATE = 25  # ml/min
VOLUMN = 0.1  # ml
SANNING_SPEED = 0.327067686193192  # min per line

# Read data
data = np.loadtxt(SOURCE_NAME, skiprows=1, usecols=(3, 5, 6))

x = np.zeros_like(data)
for i in xrange(len(x)):
    x[i,0] = i

##### Plotting #####
@mpltex.presentation_decorator
def plot(data, SOURCE_NAME):
    fig, ax = plt.subplots()
    ax.plot(x[:,0], data[:,0],label='H2O')
    ax.plot(x[:,0], data[:,1],label='CO')
    ax.plot(x[:,0], data[:,2],label='CO2')

    # ax.set_xlim(5, 80)
    ax.set_ylim(-0.1, 1.2)

    # ax.set_yticks([]) # Hide yticks
    # ax.tick_params(axis='x', top='off', bottom='off')

    ax.legend(loc='best')
    ax.set_ylabel('Concentration (\%)')

    plt.show()

##### Reports #####
# Data manipulation
TEMPERATURE = (70, 80, 90, 100, 110,
             120, 130, 140, 150, 160,
             170, 180, 190, 200, 210,
             220, 240, 260, 280, 300
             )
aveH2O = np.zeros_like(TEMPERATURE, dtype=float)
aveCO = np.zeros_like(TEMPERATURE, dtype=float)
aveCO2 = np.zeros_like(TEMPERATURE, dtype=float)
stdH2O = np.zeros_like(TEMPERATURE, dtype=float)
stdCO = np.zeros_like(TEMPERATURE, dtype=float)
stdCO2 = np.zeros_like(TEMPERATURE, dtype=float)

aveConv = np.zeros_like(TEMPERATURE, dtype=float)
stdConv = np.zeros_like(TEMPERATURE, dtype=float)
ave_k = np.zeros_like(TEMPERATURE, dtype=float)
std_k = np.zeros_like(TEMPERATURE, dtype=float)
ave_lnk = np.zeros_like(TEMPERATURE, dtype=float)
std_lnk = np.zeros_like(TEMPERATURE, dtype=float)
T_inv = np.zeros_like(TEMPERATURE, dtype=float)
linePosition = []
for i in xrange(len(TEMPERATURE)):
    time_i = 60 + i * 30 - 5
    linePosition.append(int(time_i/SANNING_SPEED))
    aveH2O[i] = np.average(data[linePosition[i]-25:linePosition[i],0])
    aveCO[i] = np.average(data[linePosition[i]-25:linePosition[i],1])
    aveCO2[i] = np.average(data[linePosition[i]-25:linePosition[i],2])
    stdH2O[i] = abs(np.std(data[linePosition[i]-25:linePosition[i],0], 
                       ddof=1))
    stdCO[i] = abs(np.std(data[linePosition[i]-25:linePosition[i],1], 
                      ddof=1))
    stdCO2[i] = abs(np.std(data[linePosition[i]-25:linePosition[i],2], 
                       ddof=1))
    # Conversion ratio
    aveConv[i] = (1 - aveCO[i] / BACKGROUND['bgCO']) * 100  # X = 1 - C / C0
    stdConv[i] = abs(stdCO[i] / BACKGROUND['bgCO'] * 100)
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
output.writelines('Temperature(C)  ' +
    'CO(%)  StandardError(%)  ' +
    'CO2(%)  StandardError(%)  ' +
    'H2O(%)  StandardError(%)  ' +
    'ConvertionRatio(%)  StandardError(%)  ' +
    'ReactionRate(s-1)  StandardError(s-1)  ' +
    'lnk  StandardError  ' +
    '1/T\n'
)
for i in xrange(len(TEMPERATURE)):
    output.writelines(
        '%3d   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   \
%6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %9.7f\n'
        % (TEMPERATURE[i], 
        aveCO[i], stdCO[i],
        aveCO2[i], stdCO2[i],
        aveH2O[i], stdH2O[i],
        aveConv[i], stdConv[i], 
        ave_k[i], std_k[i],
        ave_lnk[i], std_lnk[i],
        T_inv[i]
        )
    )
output.close()

plot(data, SOURCE_NAME)