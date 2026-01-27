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

# ODB path
pathOdb = sys.argv[-2]
#pathOdb = "b308d0disc1_staticIL1b_d_tr_swell.odb"

# desired steps
dSteps  = sys.argv[-1]
#dSteps = 'Lo_N_02_03, Cr_N_02_03, Lo_D_02_03, Cr_D_02_03'

# Define the variable to extract
var_map = {
    'UVARM1':  'NF0',
    'UVARM2':  'NF',
    'UVARM3':  'O2',
    'UVARM4':  'Lact',
    'UVARM5':  'Gluc',
    'UVARM6':  'pH',
    'UVARM7':  'Cell_rho',
    'UVARM8':  'Cell_viab',
    'UVARM9':  'nablaD',
    'UVARM10': 'SP',
    'UVARM11': 'IL1B',
    'UVARM12': 'IL1B Protein',
    'UVARM13': 'TNF',
    'UVARM14': 'TNF Protein',
    'UVARM15': 'Agg IMN',
    'UVARM16': 'Agg IL1B',
    'UVARM17': 'Agg TNF',
    'UVARM18': 'Agg IL1B TNF',
    'UVARM19': 'Col I IMN',
    'UVARM20': 'Col I IL1B',
    'UVARM21': 'Col I TNF',
    'UVARM22': 'Col I IL1B TNF',
    'UVARM23': 'Col II IMN',
    'UVARM24': 'Col II IL1B',
    'UVARM25': 'Col II TNF',
    'UVARM26': 'Col II IL1B TNF',
    'UVARM27': 'MMP3 IMN',
    'UVARM28': 'MMP3 IL1B',
    'UVARM29': 'MMP3 TNF',
    'UVARM30': 'MMP3 IL1B TNF',
    'UVARM31': 'ADAM4 IMN',
    'UVARM32': 'ADAM4 IL1B',
    'UVARM33': 'ADAM4 TNF',
    'UVARM34': 'ADAM4 IL1B TNF'
}
ordered_uvarms = ['UVARM{}'.format(i) for i in range(1, 35)]
visible_set_names = [
    "PART-1-1.DC_AF",
    "PART-1-1.DC_AF-Z1",
    "PART-1-1.DC_AF-Z2",
    "PART-1-1.DC_CEP",
    "PART-1-1.DC_NP",   
    "PART-1-1.DC_NP-Z1",
    "PART-1-1.DC_NP-Z2",
    "PART-1-1.DC_NP-Z3",
    "PART-1-1.DC_NP-Z4",
    "PART-1-1.DC_NP-Z5"
]
# Get the current working directory
current_directory = os.getcwd()

# -------------------------------------------------------------------------------------
# Create a new session
session.mdbData.summary()
o1 = session.openOdb(name=pathOdb)
session.viewports["Viewport: 1"].setValues(displayedObject=o1)
odb = session.odbs[pathOdb]

# Get step/frame info
dir_path  = os.path.dirname(pathOdb)    # → '/home/estefano/Documents/phd/tesis/abqsSim/results/generic/newNablaD_ani'
file_name = os.path.basename(pathOdb)   # → 'GENERIC_Transp.odb'
steps = odb.steps.keys()

# create the folder CSVs if it does not exist in the dir_path
csv_folder = os.path.join(dir_path, 'CSVs')
if not os.path.exists(csv_folder):
    try:
        os.makedirs(csv_folder)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise  # re-raise the exception if it's not about the directory already existing
# build a clean Python list from the comma-separated string -------------------
desired_steps = [s.strip() for s in dSteps.split(',')]          # → ['Lo_N_02_03', 'Cr_N_02_03', 'Lo_D_02_03', 'Cr_D_02_03']

# loop only over the desired steps -------------------------------------------
for stepName in desired_steps:
    if stepName not in steps:
        print('Step {:s} not found in the ODB'.format(stepName))
        continue
    step_number  = list(steps).index(stepName)                  # actual index in the ODB
    stepObj      = odb.steps[stepName]
    last_frame   = len(stepObj.frames) - 1                      # last frame index (0-based)    
    print('Step name: {:<15}  Step number: {:>2}  Last-frame index: {}'.format(stepName, step_number, last_frame))

    # Collect visible nodes
    visible_nodes = set()
    instance = odb.rootAssembly.instances["PART-1-1"]
    for full_set_name in visible_set_names:
        set_name = full_set_name.split(".")[-1]
        try:
            elset = instance.elementSets[set_name]
            for el in elset.elements:
                for node_label in el.connectivity:
                    visible_nodes.add(node_label)
        except KeyError:
            continue

    visible_node_list = sorted(list(visible_nodes))
    print("Total NP node count:", len(visible_node_list))
    print("First 10 node IDs:", visible_node_list[:10])

    # Create report
    csv_filename = os.path.splitext(pathOdb)[0] + '_' + stepName + '.csv'
    csv_filename_avg = os.path.splitext(pathOdb)[0] + '_' + stepName + '_avg.csv'
    
    # the csv should have the path to the CSVs folder
    csv_filename = os.path.join(csv_folder, os.path.basename(csv_filename))
    csv_filename_avg = os.path.join(csv_folder, os.path.basename(csv_filename_avg))
    print("CSV file will be saved as:", csv_filename)
    print("CSV file with averaged data will be saved as:", csv_filename_avg)
    
    nf = NumberFormat(numDigits=9, precision=0, format=ENGINEERING)
    session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES, numberFormat=nf)
    session.writeFieldReport(
        fileName=csv_filename,
        append=OFF,
        sortItem='Node Label',
        odb=odb,
        step=step_number,
        frame=last_frame,
        outputPosition=NODAL,
        variable=tuple([(varname, INTEGRATION_POINT) for varname in ordered_uvarms]),
        stepFrame=SPECIFY
    )

    # ----------------------------
    # Post-process CSV
    # ----------------------------
    with open(csv_filename, 'r') as f:
        reader = csv.reader(f)
        all_rows = list(reader)

    header = all_rows[0]
    node_col_idx = None
    col_map = {}

    # Find column indices
    # ----- revised header scan -----
    node_col_idx = None
    col_map = {}

    for i, col in enumerate(header):
        # catch the node-label column
        if 'Node Label' in col:
            node_col_idx = i
            continue                    # done with this column

        # first token in the header cell (e.g. 'UVARM1', 'UVARM10', ...)
        token = col.strip().split()[0]

        # record the column only if it is an exact match
        if token in ordered_uvarms and token not in col_map:
            col_map[token] = i
    # --------------------------------


    if node_col_idx is None:
        raise RuntimeError("Node Label column not found in CSV header.")

    # Build ordered column list
    ordered_keep_cols = [node_col_idx] + [col_map[var] for var in ordered_uvarms if var in col_map]

    # Filter and reorder
    visible_node_set = set(visible_node_list)
    filtered_rows = [[header[i] for i in ordered_keep_cols]]
    for row in all_rows[1:]:
        try:
            node_label = int(float(row[node_col_idx].strip()))
            if node_label in visible_node_set:
                filtered_rows.append([row[i] for i in ordered_keep_cols])
        except:
            continue

    # Average duplicates
    node_data = defaultdict(list)
    num_vars = len(ordered_keep_cols) - 1

    for row in filtered_rows[1:]:
        try:
            node_label = int(float(row[0]))
            values = [float(row[i]) for i in range(1, len(row))]
            if len(values) == num_vars:
                node_data[node_label].append(values)
        except:
            continue

    averaged_rows = [filtered_rows[0]]
    for node_label in sorted(node_data.keys()):
        all_values = node_data[node_label]
        avg_values = []
        for i in range(num_vars):
            total = 0.0
            for vals in all_values:
                total += vals[i]
            avg_values.append(total / len(all_values))
        averaged_rows.append([node_label] + avg_values)

    # Overwrite the CSV
    with open(csv_filename_avg, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(averaged_rows)

    print("Cleaned CSV saved in place:", csv_filename_avg)
    print("Remaining node count in file:", len(averaged_rows) - 1)