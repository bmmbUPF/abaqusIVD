#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#  ___  ___  ________     ___    ___ ___  ___  ________     ___    ___ ___    ___  _______  ________
# |\  \|\  \|\   __  \   |\  \  /  /|\  \|\  \|\   __  \   |\  \  /  /|\  \  /  /|/  ___  \|\_____  \
# \ \  \\\  \ \  \|\  \  \ \  \/  / | \  \\\  \ \  \|\  \  \ \  \/  / | \  \/  / /__/|_/  /\|____|\ /_
#  \ \   __  \ \   __  \  \ \    / / \ \   __  \ \   __  \  \ \    / / \ \    / /|__|//  / /     \|\  \
#   \ \  \ \  \ \  \ \  \  /     \/   \ \  \ \  \ \  \ \  \  /     \/   /     \/     /  /_/__   __\_\  \
#    \ \__\ \__\ \__\ \__\/  /\   \    \ \__\ \__\ \__\ \__\/  /\   \  /  /\   \    |\________\|\_______\
#     \|__|\|__|\|__|\|__/__/ /\ __\    \|__|\|__|\|__|\|__/__/ /\ __\/__/ /\ __\    \|_______|\|_______|
#                        |__|/ \|__|                       |__|/ \|__||__|/ \|__|
#          01001000 01100001 01111000 01001000 01100001 01111000 01111000 00110010 00110011
#

#
# ******************************************************************************
# ******    Extract Initial solutes values from previous simulation       ******
# ******    auth: Estefano Muñoz-Moya                                     ******
# ******    LinkTree: https://linktr.ee/estefano23                        ******
# ******    webPage: https://estefano23.github.io/                        ******
# ******    github: estefano23                                            ******
# ******    email: estefano.munoz.moya@gmail.com                          ******
# ******************************************************************************
#
#
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Introduction
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
# libraries
from __future__ import print_function
import os
import sys
import argparse
import datetime
import numpy as np  
import pandas as pd
#
###
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# code
# -----------------

# read the file solutes.csv
# headers: "     Node Label", "         UVARM3", "         UVARM4", "         UVARM5"


def main():
    # Read the CSV file
    df = pd.read_csv("solutes.csv")

    # Remove extra whitespace from the column names
    df.columns = df.columns.str.strip()

    # If your file contains exactly these columns:
    # "Node Label", "UVARM3", "UVARM4", "UVARM5"
    # then group by "Node Label" and compute the mean for each solute separately.
    df_averaged = df.groupby("Node Label", as_index=False)[
        ["UVARM3", "UVARM4", "UVARM5"]
    ].mean()

    # Save the averaged values to a new CSV file
    df_averaged.to_csv("solutes_averaged.csv", index=False)
    print("Averaged solute values have been saved to 'solutes_averaged.csv'.")


if __name__ == "__main__":
    main()
