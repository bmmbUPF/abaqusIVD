#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
 ___  ___  ________     ___    ___ ___  ___  ________     ___    ___ ___    ___  _______  ________
|\  \|\  \|\   __  \   |\  \  /  /|\  \|\  \|\   __  \   |\  \  /  /|\  \  /  /|/  ___  \|\_____  \
\ \  \\\  \ \  \|\  \  \ \  \/  / | \  \\\  \ \  \|\  \  \ \  \/  / | \  \/  / /__/|_/  /\|____|\ /_
 \ \   __  \ \   __  \  \ \    / / \ \   __  \ \   __  \  \ \    / / \ \    / /|__|//  / /     \|\  \
  \ \  \ \  \ \  \ \  \  /     \/   \ \  \ \  \ \  \ \  \  /     \/   /     \/     /  /_/__   __\_\  \
   \ \__\ \__\ \__\ \__\/  /\   \    \ \__\ \__\ \__\ \__\/  /\   \  /  /\   \    |\________\|\_______\
    \|__|\|__|\|__|\|__/__/ /\ __\    \|__|\|__|\|__|\|__/__/ /\ __\/__/ /\ __\    \|_______|\|_______|
                       |__|/ \|__|                       |__|/ \|__||__|/ \|__|


"""


def read_ids_from_file(filename):
    with open(filename, 'r') as f:
        content = f.read()
        # Replace newlines with commas and then split the content by comma
        cleaned_content = content.replace('\n', ',').replace(' ', '')
        return [int(id.strip()) for id in cleaned_content.split(',') if id.strip()]



def write_ids_to_file(filename, ids):
    with open(filename, 'w') as f:
        for idx, node_id in enumerate(ids):
            # For formatting, write a newline every 16 IDs
            if idx and idx % 16 == 0:
                f.write('\n')
            f.write(f"{node_id:8},")  # Format the ID with 8 spaces width

# Read node IDs from the source file
node_ids = read_ids_from_file('TransportDisc_nset_IniCondNodes.inp')

# List of node IDs to remove (change this list based on your needs)
ids_to_remove = read_ids_from_file('TransportDisc_nset_ND_SURF_AF.inp')

# Remove undesired node IDs
filtered_ids = [node_id for node_id in node_ids if node_id not in ids_to_remove]

# Write the result to a new file
write_ids_to_file('filtered_file.txt', filtered_ids)
