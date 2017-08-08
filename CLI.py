# -*- coding: utf-8 -*-

from __future__ import absolute_import, print_function, division
import argparse
import os
import numpy as np

def parseRPArgs(): 
    """ 
    Defines command line argmuments for this Ramaplot program
    """

    parser = argparse.ArgumentParser(prog="Ramaplot", description='Generate Ramachandran plots for a given pdb input file.')

    parser.add_argument("pdb_file", help="Input pdb file")

    args = parser.parse_args()
    args = vars(args)

    return args



