#!/usr/bin/env python

"""
Reaction rate post-treatment script
"""

__author__ = "LI Kezhi"
__date__ = "$2017-07-04$"
__version__ = "1.0.5"


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
FLOW_RATE = 83.33  # ml/min
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

# Data manipulation
FUNCTIONS_NAME = (
    'conversionratio',
    'reactionrate',
    'lnr',
    'T_inv',
)
FUNCTIONS_DISPLAY = (
    'ConversionRatio(%)',
    'ReactionRate(mol s-1 m-2)',
    'lnr',
    r'1/T',
)
def conversionratio(conc, bg):
    '''
    return conversion ratio (%)
    '''
    result = (1 - conc / bg) * 100
    try:
        assert result < 100
    except AssertionError:
        result = 99.99999
    return result
def reactionrate(conv, conc):
    '''
    return reaction rate (mol s-1 m-2)
    '''
    try:
        result = -FLOW_RATE * (1E-3/60) * (1/22.4) / (MASS * BET) * math.log(1 - conv/100) * conc
        # k = - F/(m * S_BET) * ln(1 - X) * C
    except ValueError:
        result = np.NaN
    return result
def lnr(rate):
    '''
    return log of reaction rate
    '''
    try:
        result = math.log(rate)
    except ValueError:
        result = np.NaN
    return result
def T_inv(temp):
    '''
    return inverse of temperature
    '''
    if temperature[-1] is None:
        result = np.NaN
    else:
        result = 1 / (273.15 + temp)
    return result

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
temperature = []
calculationResults = {}
for item in FUNCTIONS_NAME:
    calculationResults[item] = []

for lineNum, line in enumerate(data):
    temperature.append(
        time2temprature(lineNum * SCANNING_SPEED)
    )
    for i, item in enumerate(NAME):
        concentration[item].append(line[i])
    # Calculations
    calculationResults['conversionratio'].append(
        conversionratio(concentration[NAME[0]][-1], BACKGROUND[NAME[0]])
    )
    calculationResults['reactionrate'].append(
        reactionrate(calculationResults['conversionratio'][-1], concentration[NAME[0]][-1])
    )
    calculationResults['lnr'].append(
        lnr(calculationResults['reactionrate'][-1])
    )
    calculationResults['T_inv'].append(
        T_inv(temperature[-1])
    )

# Generate Report
reportName = SOURCE_NAME.split('.txt')[0] + '_report.txt'
output = file(reportName, 'w')
headLine = 'Temperature(C)  '
for item in NAME:
    headLine += item
    headLine += '(ppm)  '
for item in FUNCTIONS_DISPLAY:
    headLine += item + '   '
headLine += '\n'
output.write(headLine)
for i, temp_i in enumerate(temperature):
    if temp_i is None:
        continue
    dataLine = '%3d   ' % temp_i
    for item in NAME:
        dataLine += '%9.4f   ' % concentration[item][i]
    for item in FUNCTIONS_NAME:
        if item is 'reactionrate':
            dataLine += '%6.4e   ' % calculationResults[item][i]
        else:
            dataLine += '%9.7f   ' % calculationResults[item][i]
    dataLine += '\n'
    output.write(dataLine)
output.close()

plot(data, SOURCE_NAME)
