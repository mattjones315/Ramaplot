# -*- coding: utf-8 -*-
"""Entry point for console fastproject script

- Parses arguments
- Loads data from files
- Runs analysis
- Writes outputs to file

Then launches the main pipeline
"""

from __future__ import absolute_import, print_function, division;
import argparse;
import os
import numpy as np;
import Ramaplots
import FileIO

def parseRPArgs(): 
    """ 
    Defines command line argmuments for this Ramaplot program
    """

    parser = argparse.ArgumentParser(prog="Ramaplot", description='Generate Ramachandran plots for a given pdb input file.')

    parser.add_argument("pdb_file", help="Input pdb file")

    args = parser.parse_args()

    return args


def entry(): 
    """
    Entry point for Ramaplot command-line script
    """

    # Read in arguments from command line 
    args = parseRPArgs()

    # Create protein data object to compute dihedral angles on 
    protein = readPDB(args["pdb_file"])
    
    # Compute dihedral angles, saved as a dictionary object
    angles = ComputeAngles(protein)

    # Now plot results
    


