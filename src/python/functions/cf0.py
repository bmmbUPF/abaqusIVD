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
desireExtension = 'cF0'
outPutFile = transName + '_SDV.' + desireExtension
# add the path to the odb file
Step = 'Swelling'
variable = 'SDV3'

# Get the current working directory
current_directory = os.getcwd()

# -------------------------------------------------------------------------------------
# Create a new session
session.mdbData.summary()
o1 = session.openOdb(name=pathOdb)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)

# check last increment of the step
lastInc = session.odbs[pathOdb].steps[Step].frames[-1]
incrementNumber = lastInc.incrementNumber

# Get the name of the odb
odbName=session.viewports[session.currentViewportName].odbDisplay.name
session.odbData[odbName].setValues(activeFrames=((Step, (incrementNumber, )), ))

# Opening the odb file
odb = session.odbs[pathOdb]
cf = NumberFormat(numDigits=9, precision=17, format=SCIENTIFIC)
session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES, numberFormat=cf)

# Output file name
odbNameWE = odbName.split('.')[0]
outPutCSVFile = odbNameWE + '.' + desireExtension + '.csv'

# avaible frame indices
frameIndices = session.odbs[pathOdb].steps[Step].frames
lastFrame = frameIndices[-1]

# get the step number
stepNumber = session.odbs[pathOdb].steps[Step].number - 1

# Write the field report
session.writeFieldReport(fileName=outPutCSVFile, append=OFF, sortItem='Element Label', odb=odb, step=stepNumber, frame=lastFrame, outputPosition=INTEGRATION_POINT, variable=((variable, INTEGRATION_POINT), ), stepFrame=ALL)

# edit the csv file to deletes all columns except the one with the
# header matching the string stored in the variable variable
# Read the CSV file and process it
with open(outPutCSVFile, 'r') as file:
  variableData = list()
  csvreader = csv.reader(file)
  header = next(csvreader)

  # check the index of the header that contains the string in variable
  indexVariable = [index for index, column_name in enumerate(header) if variable in column_name]
  # check the index of the header that contains the string "Element Label"
  indexElementLabel = [index for index, column_name in enumerate(header) if "Element Label" in column_name]
  # check the index of the header that contains the string "IntPt"
  indexIntPt = [index for index, column_name in enumerate(header) if "IntPt" in column_name]
  
  # extract the column with the header with index indexVariable
  for row in csvreader:
    variableData.append([row[indexElementLabel[0]], row[indexIntPt[0]], row[indexVariable[0]]])

# write the data in a new file called outPutFile
with open(outPutFile, 'w') as file:
  for row in variableData:
    # delete every space in the string
    row = [element.replace(' ', '') for element in row]
    # Join the elements of the list into a single string separated by a delimiter (e.g., comma)
    row_string = ','.join(row)
    file.write(row_string + '\n')
    
# delete the old csv file
os.remove(outPutCSVFile)
