# -*- coding: mbcs -*-

"""
  ___  ___  ________     ___    ___ ___  ___  ________     ___    ___ ___    ___  _______  ________
 |\  \|\  \|\   __  \   |\  \  /  /|\  \|\  \|\   __  \   |\  \  /  /|\  \  /  /|/  ___  \|\_____  \
 \ \  \\\  \ \  \|\  \  \ \  \/  / | \  \\\  \ \  \|\  \  \ \  \/  / | \  \/  / /__/|_/  /\|____|\ /_
  \ \   __  \ \   __  \  \ \    / / \ \   __  \ \   __  \  \ \    / / \ \    / /|__|//  / /     \|\  \
   \ \  \ \  \ \  \ \  \  /     \/   \ \  \ \  \ \  \ \  \  /     \/   /     \/     /  /_/__   __\_\  \
    \ \__\ \__\ \__\ \__\/  /\   \    \ \__\ \__\ \__\ \__\/  /\   \  /  /\   \    |\________\|\_______\
     \|__|\|__|\|__|\|__/__/ /\ __\    \|__|\|__|\|__|\|__/__/ /\ __\/__/ /\ __\    \|_______|\|_______|
                        |__|/ \|__|                       |__|/ \|__||__|/ \|__|
          01001000 01100001 01111000 01001000 01100001 01111000 01111000 00110010 00110011
"""

#
# ***************************************************************************
# ******    Program to obtain a desire variable from a odb file        ******
# ******    auth: Estefano Muñoz-Moya                                  ******
# ******    LinkTree: https://linktr.ee/estefano23                     ******
# ******    webPage: https://estefano23.github.io/                     ******
# ******    github: estefano23                                         ******
# ******    email: estefano.munoz.moya@gmail.com                       ******
# ***************************************************************************
#

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# introduction
# ------------

# libraries
from __future__ import print_function
from abaqus import *
from abaqusConstants import *
import __main__
import os
import sys
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import odbAccess
import csv
import numpy as np
import datetime
import subprocess
from collections import defaultdict
import errno
import re
###

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# code
#-----------------

pathOdb = sys.argv[-2]
transName = os.path.splitext(sys.argv[-1])[0]
desireExtension = 'NF0'
outPutFile = transName + '_SDV.' + desireExtension
# add the path to the odb file
Step = 'Swelling'
variable = 'SDV2'

# Get the current working directory
current_directory = os.getcwd()

# -------------------------------------------------------------------------------------
# Create a new session
session.mdbData.summary()
o1 = session.openOdb(name=pathOdb)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)

# Get the name of the odb
odbName=session.viewports[session.currentViewportName].odbDisplay.name

# Output file name
odbNameWE = odbName.split('.')[0]
outPutCSVFile = odbNameWE + '.' + desireExtension + '.csv'

# Write the field report
stepsODB = session.odbs[pathOdb].steps
bulkBlock = stepsODB[Step].frames[-1].fieldOutputs[variable].bulkDataBlocks[0]

# Extract the data, element labels, and integration point IDs.
data_values = bulkBlock.data
element_labels = bulkBlock.elementLabels
integration_points = bulkBlock.integrationPoints

# Extract the values, element labels, and integration point IDs as lists.
values = data_values[:, 0]  # gets all values at once
values_list = values.tolist()
element_labels_list = element_labels.tolist()
integration_points_list = integration_points.tolist()

# Build the list of variable data
variableData = [[str(el), str(ip), str(val)] for el, ip, val in zip(element_labels_list, integration_points_list, values_list)]

# Write the data in a new file called outPutFile
with open(outPutFile, 'w') as file:
    for row in variableData:
        row_string = ','.join(row)
        file.write(row_string + '\n')