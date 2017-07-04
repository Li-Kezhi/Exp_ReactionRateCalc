#!/usr/bin/env python

"""
Reaction rate post-treatment script for continuous experiments
"""

__author__ = "LI Kezhi"
__date__ = "$2017-07-04$"
__version__ = "1.0.4"


import math
import numpy as np
import matplotlib.pyplot as plt

import mpltex


#####################################
##### User Input #####
# File names
SOURCE_NAME = './Fe.txt'

# Experiment details
NAME = ('NO', 'NH3', 'N2O', 'NO2') # The conversion ratio of 1st gas is calculated!
BACKGROUND = dict(zip(NAME, (503, 552, 1.5, 19)))
FLOAT_RATE = 83.33  # ml/min
MASS = 0.05  # g
BET = 100  # m^2/g
CONCENTRATION = BACKGROUND[NAME[0]] * 1e-6
SCANNING_SPEED = 0.420432044105174  # min per line
def time2temprature(minute):
    '''
    retun the corresponding temperature for a given time
    '''
    if 60 <= minute < 65:
        return 75
    elif 65 <= minute < 290:
        return 75 + (minute - 65)
    elif 290 <= minute < 295:
        return 300
    elif 295 <= minute < 520:
        return 300 - (minute - 295)
    elif 520 <= minute < 525:
        return 75
    else:
        return None

# File formatting
DATA_ROW = (4, 3, 5, 6) # Correspond to NAME sequence

# Report
CALCULATE_REACTION_RATE = False

#####################################


##### Read data #####
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
concentration = {} # concentration
for item in NAME:
    concentration[item] = []
conversionRatio = []
temperature = []
if CALCULATE_REACTION_RATE is True:
    reactionRate = []
    lnr = []
    T_inv = []

for lineNum, line in enumerate(data):
    temperature.append(time2temprature(lineNum * SCANNING_SPEED))
    for i, item in enumerate(NAME):
        concentration[item].append(line[i])

    # Conversion ratio
    conversionRatio.append((1 - concentration[NAME[0]][-1] / BACKGROUND[NAME[0]]) * 100)  # X = 1 - C / C0
    try:
        assert conversionRatio[-1] < 100
    except AssertionError:
        conversionRatio[-1] = 99.99999
    if CALCULATE_REACTION_RATE is True:
        # Reaction rate constant
        try:
            reactionRate.append(
                -FLOAT_RATE * (1E-3/60) * (1/22.4) / (MASS * BET) * math.log(1 - conversionRatio[-1]/100) * CONCENTRATION
            )
            # k = - F/(m * S_BET) * ln(1 - X) * C
        except ValueError:
            reactionRate.append(np.NaN)
        # lnr - 1/T
        try:
            lnr.append(math.log(reactionRate[-1]))
        except ValueError:
            lnr.append(np.NaN)
        if temperature[-1] is None:
            T_inv.append(np.NaN)
        else:
            T_inv.append(1 / (273.15 + temperature[-1]))

# Generate Report
reportName = SOURCE_NAME.split('.txt')[0] + '_report.txt'
output = file(reportName, 'w')
headLine = 'Temperature(C)  '
for item in NAME:
    headLine += item
    headLine += '(ppm)  '
headLine += 'ConvertionRatio(%)  '
if CALCULATE_REACTION_RATE is True:
    headLine += (
        'ReactionRate(mol s-1 m-2)  ' +
        'lnr   ' +
        '1/T'
    )
headLine += '\n'
output.write(headLine)
for i, temp_i in enumerate(temperature):
    if temp_i is None:
        continue
    dataLine = '%3d   ' % temp_i
    for item in NAME:
        dataLine += '%9.4f   ' % concentration[item][i]
    dataLine += '%6.4f   ' % conversionRatio[i]
    if CALCULATE_REACTION_RATE is True:
        dataLine += '%6.4e   ' % reactionRate[i]
        dataLine += '%6.4f   ' % lnr[i]
        dataLine += '%9.7f' % T_inv[i]
    dataLine += '\n'
    output.write(dataLine)
output.close()

plot(data, SOURCE_NAME)
