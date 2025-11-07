# make paths available:
from enum import Enum
import json
import sys

class Paths(Enum):
    '''This Enum stores the dict with all important paths on UPPMAX.
    '''
    DICT = json.load(open("../../PATHS.json", "r"))
    
# make the slurmscripts accessible
class Scripts(Enum):
    '''This Enum stores the dict with paths to all simple SLURM scripts.
    '''
    DICT = json.load(open("../../SLURMSCRIPTS.json", "r"))

# access some classes from other dirs in this directory.
sys.path.insert(0, "../../src")
