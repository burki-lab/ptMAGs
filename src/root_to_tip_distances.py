#!/usr/bin/python

## Mahwash Jamy
## Oct 2024

## This script roots a tree at a specified node, and then calculates the root to tip distance for all tips.

## Usage: python root_to_tip_distances.py [tree newick file] [clade name] [output file] 

from ete3 import Tree
import sys

IN=sys.argv[1]
CLADE=sys.argv[2]

TIPS = []

# read tree
t = Tree( IN )

# Get the list of tips that contain the clade name
for leaf in t:
	if CLADE in leaf.name:
		TIPS.append(leaf.name)	

ancestor = t.get_common_ancestor(TIPS)

# Set as outgroup
t.set_outgroup(ancestor)

# Calculate the root to tip distance for all tips
with open(sys.argv[3], 'w') as f:
	for leaf in t.iter_leaves():
        # Compute the distance from the root to the current tip
        	distance = leaf.get_distance(t)
        	f.write(f"{leaf.name}\t{distance}\n")


