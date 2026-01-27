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

#----------------------------------------------------------------------- --------------------------------------------------------------------------------------------
# code
#-----------------

pathOdb = sys.argv[-2]
transName = os.path.splitext(sys.argv[-1])[0]
desireExtension = 'SP'
variable = 'UVARM1'

# Get the current working directory
current_directory = os.getcwd()

# -------------------------------------------------------------------------------------
# Functions
def extract_time_steps(file_path):
	"""
	Extracts the time steps from a given input file.

	Args:
	file_path (str): Path to the input file.

	Returns:
	dict: A dictionary where keys are step numbers and values are lists containing the index and the time step as a string.
	"""
	time_steps = dict()  # Initialize the dictionary to store time steps
	step = 1  # Initialize the step counter
	with open(file_path, "r") as f:
		lines = f.readlines()  # Read all lines into a list
	for idx, line in enumerate(lines):
		# Check if the line contains '** TIME-STEP'
		if "** TIME-STEP" in line:
			# The next line contains the time step
			next_line_index = idx + 1  # Index of the next line
			next_line = [float(x) for x in lines[next_line_index].split(',') if x.strip()]  # Remove newline character
			time_steps[step] = dict()
			time_steps[step]['time'] = next_line
			# Increment the step counter
			step += 1
        # Check if the line contains the boundary info with a step number
		elif "*BOUNDARY" in line and "SUBMODEL" in line and "STEP" in line:
			# The actual line contains the boundary step as STEP=XX
			globalStep = int(line.split("=")[-1])
			time_steps[step-1]["globalStep"] = globalStep
	return time_steps


def create_step_mapping(stepsODB, global_steps):
    # Create a new dictionary with keys as step numbers (starting from 1)
    stepDict = dict()
    totalTime = 0.0
    # global_steps is a list
    for index, step_name in enumerate(stepsODB.keys()):
        # check if the index + 1 in in the global step
        if index + 1 in global_steps:
          stepDict[index + 1] = dict()
          stepDict[index + 1]['name'] = step_name
          stepDict[index + 1]['n_frames'] = len(stepsODB[stepDict[index + 1]['name']].frames)
          stepDict[index + 1]['frames'] = list()
          for iframe, frame in enumerate(stepsODB[stepDict[index + 1]['name']].frames):
              increment = iframe
              stepTime = float(frame.frameValue)
              stepDict[index + 1]['frames'].append([increment, stepTime])
          totalTime += stepTime
    # remove every key from stepDict that is not in global_steps
    stepDict = {k: v for k, v in stepDict.items() if k in global_steps}
    return stepDict

# -------------------------------------------------------------------------------------
# Read the transport .inp file to get the desired step number that
# we want to extract from mechanical simulation
transNameInp = transName + '.inp'
time_steps = extract_time_steps(transNameInp)
# remove every key without the globalStep key
time_steps = {key: value for key, value in time_steps.items() if 'globalStep' in value}

# -------------------------------------------------------------------------------------
# Create a new session
session.mdbData.summary()
o1 = session.openOdb(name=pathOdb)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)

# Get the name of the odb
odbName=session.viewports[session.currentViewportName].odbDisplay.name

# -------------------------------------------------------------------------------------
# Get the Steps of the ODB
stepsODB = session.odbs[pathOdb].steps

# Build a set of valid global steps from time_steps.
global_steps = list({value['globalStep'] for value in time_steps.values()})

# Here we use the simulated 'steps' dictionary:
stepDict = create_step_mapping(stepsODB, global_steps)

# -------------------------------------------------------------------------------------
# Get the required steps and frames (reqFrames) and build interpMapping
#   - reqFrames: keys are previous simulation steps (2–9)
#   - interpMapping: keys are new simulation steps (from time_steps, e.g. 2–13)
reqFrames = dict()
interpMapping = dict()  # interpMapping[new_step][req_index] = [prev_step, frame] or [prev_step, lower_frame, upper_frame]
tolerance = 1e-6  # tolerance for matching floating-point values

for new_step in sorted(time_steps.keys()):
    # For each new simulation step, get its transport info and determine the previous simulation step.
    ts_info = time_steps[new_step]
    dt = ts_info['time'][0]         # time increment for this new simulation step
    T_final = ts_info['time'][-1]     # final time for this new simulation step
    n = int(round(T_final / dt))
    required_times = [i * dt for i in range(n + 1)]
    
    # Get the previous simulation step (globalStep from transport file)
    prev_step = ts_info['globalStep']  # expected to be in the range 2–9
    
    # Ensure reqFrames for the previous simulation step exists.
    if prev_step not in reqFrames:
        reqFrames[prev_step] = list()
    
    # For the new simulation step, create an entry in interpMapping.
    interpMapping[new_step] = dict()
    
    # Get the list of frames from the previous simulation step (from stepDict).
    # Each frame is stored as: [frame, stepTime, totalTime]
    frames_list = stepDict[prev_step]['frames']
    
    # Compute relative times from the previous simulation frames (relative to the first frame's totalTime)
    base_time = frames_list[0][1]
    rel_times = [frame[1] - base_time for frame in frames_list]
    
    # For each required time from the new simulation, find the corresponding frame(s)
    for req_index, r in enumerate(required_times):
        match_index = None
        for i, t in enumerate(rel_times):
            if abs(t - r) < tolerance:
                match_index = i
                break
        
        if match_index is not None:
            # Exact match found: use that frame.
            frame_item = frames_list[match_index]
            if frame_item not in reqFrames[prev_step]:
                reqFrames[prev_step].append(frame_item)
            interpMapping[new_step][req_index] = [prev_step, frame_item[0]]
        else:
            # No exact match: search for two consecutive frames that bracket r.
            lower = None
            upper = None
            for i in range(len(rel_times) - 1):
                if rel_times[i] < r < rel_times[i + 1]:
                    lower = i
                    upper = i + 1
                    break
            if lower is None or upper is None:
                # Fallback: if r is outside the range, use the first or last frame.
                if r < rel_times[0]:
                    frame_item = frames_list[0]
                    if frame_item not in reqFrames[prev_step]:
                        reqFrames[prev_step].append(frame_item)
                    interpMapping[new_step][req_index] = [prev_step, frame_item[0]]
                else:
                    frame_item = frames_list[-1]
                    if frame_item not in reqFrames[prev_step]:
                        reqFrames[prev_step].append(frame_item)
                    interpMapping[new_step][req_index] = [prev_step, frame_item[0]]
            else:
                # Use both lower and upper frames.
                lower_frame = frames_list[lower]
                upper_frame = frames_list[upper]
                if lower_frame not in reqFrames[prev_step]:
                    reqFrames[prev_step].append(lower_frame)
                if upper_frame not in reqFrames[prev_step]:
                    reqFrames[prev_step].append(upper_frame)
                interpMapping[new_step][req_index] = [prev_step, lower_frame[0], upper_frame[0]]

# Remove duplicate frames within each new simulation step in reqFrames (preserving order)
for new_step in reqFrames:
    unique_frames = list()
    for frame in reqFrames[new_step]:
        if frame not in unique_frames:
            unique_frames.append(frame)
    reqFrames[new_step] = unique_frames

# -------------------------------------------------------------------------------------
# Remove the non-essential tissues of the odb:
session.linkedViewportCommands.setValues(_highlightLinkedViewports=True)
leaf = dgo.LeafFromElementSets(elementSets=("L4-L5-1-1.AF_REBARS", "L4-L5-1-1.BEPBOTTOM", "L4-L5-1-1.BEPTOP", "L4-L5-1-1.REBARS_AI", "L4-L5-1-1.REBARS_AO", "L4-L5-1-1.REBARS_PI", "L4-L5-1-1.REBARS_PO", ))
session.viewports['Viewport: 1'].odbDisplay.displayGroup.remove(leaf=leaf)

# -------------------------------------------------------------------------------------
# create a temporary folder to store the output files
baseName = os.path.basename(transName)  # e.g., 'GENERIC_Transport_DiffReac'
tempFolder = transName + '_' + desireExtension
if not os.path.exists(tempFolder):
    os.makedirs(tempFolder)

# For each new simulation step, process the interpMapping directly.
for new_step in sorted(time_steps.keys()):
    prev_step = time_steps[new_step]['globalStep']
    ts_info = time_steps[new_step]
    dt = ts_info['time'][0]
    T_final = ts_info['time'][-1]
    n = int(round(T_final / dt))
    required_times = [i * dt for i in range(n + 1)]
    
    if new_step not in interpMapping:
        continue
    
    for req_index, mapping in interpMapping[new_step].items():
        # Construct final output filename.
        new_file = os.path.join(tempFolder, baseName + '_SDV.' + str(new_step) + '.' + str(req_index) + '.' + desireExtension)
        
        stepName = stepDict[prev_step]['name']
        
        # --- Case 1: Exact match ---
        if len(mapping) == 2 and req_index != 0:
            exact_frame = mapping[1]
            bulkBlock = stepsODB[stepName].frames[exact_frame].fieldOutputs[variable].bulkDataBlocks[0]
            data_values = bulkBlock.data
            element_labels = bulkBlock.elementLabels
            integration_points = bulkBlock.integrationPoints
            
            values_list = data_values[:, 0].tolist()
            element_labels_list = element_labels.tolist()
            integration_points_list = integration_points.tolist()
            
            variableData = [[str(el), str(ip), str(val)] for el, ip, val in zip(element_labels_list, integration_points_list, values_list)]
            
            with open(new_file, 'w') as fout:
                for row in variableData:
                    fout.write(','.join(row) + '\n')
        
        # --- Case 2: Interpolation needed ---
        elif len(mapping) == 3 and req_index != 0:
            lower_frame = mapping[1]
            upper_frame = mapping[2]
            req_time = required_times[req_index]
            
            # Get frame data from stepDict for interpolation factor.
            frames_list = stepDict[prev_step]['frames']
            base_time = frames_list[0][1]
            lower_data = next((data for data in frames_list if data[0] == lower_frame), None)
            upper_data = next((data for data in frames_list if data[0] == upper_frame), None)
            if lower_data is None or upper_data is None:
                continue
            lower_rel = lower_data[1] - base_time
            upper_rel = upper_data[1] - base_time
            f = (req_time - lower_rel) / (upper_rel - lower_rel) if (upper_rel - lower_rel) != 0 else 0.0
            
            # Extract field data for lower frame.
            bulkBlock_lower = stepsODB[stepName].frames[lower_frame].fieldOutputs[variable].bulkDataBlocks[0]
            lower_values = bulkBlock_lower.data[:, 0].tolist()
            element_labels_list = bulkBlock_lower.elementLabels.tolist()
            integration_points_list = bulkBlock_lower.integrationPoints.tolist()
            
            # Extract field data for upper frame.
            bulkBlock_upper = stepsODB[stepName].frames[upper_frame].fieldOutputs[variable].bulkDataBlocks[0]
            upper_values = bulkBlock_upper.data[:, 0].tolist()
            
            # Interpolate each value.
            variableData = []
            for el, ip, low_val, up_val in zip(element_labels_list, integration_points_list, lower_values, upper_values):
                new_val = float(low_val) + f * (float(up_val) - float(low_val))
                variableData.append([str(el), str(ip), str(new_val)])
            
            with open(new_file, 'w') as fout:
                for row in variableData:
                    fout.write(','.join(row) + '\n')