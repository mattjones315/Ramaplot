"""
Contains functions and objects related to representing proteins in this Ramaplot module.
"""

from __future__ import absolute_import, print_function, division
import argparse
import os
import numpy as np
import CLI


class Atom:

    def __init__(self, serial_number, name, x, y, z):
        """
        In our Atom class, we want to make sure we'll be able to identify the atom
        and locate it spatially. Thus, we'll take into account it's serial number & name 
        as well as its coordinates. These will be grouped into amino acids, which will
        then constitute proteins.
        """

        self.serial_number = serial_number
        self.name = name
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        rep = self.serial_number + " " + self.name + \
                "(" + self.x + ", " + self.y + ", " + self.z + ")
        return rep

class Residue:

    def __init__(self, name, sequence_number):
        self.name = name
        self.seq_num = sequence_number

        self.atoms = []

    def addAtom(self, atom):
        self.atoms.append(atom)

    def __str__(self):
        s = ""
        for a in self.atoms:
            s += str(a) + "\n"
        return s

class Protein:

    def __init__(self, name):
        self.name = name
        self.amino_acids = []
    
    def addResidue(self, aa):
        self.amino_acids.append(aa)

    def __str__(self):
        s = ""
        for r in self.amino_acids:
            s += r.name + " " + r.seq_num + "\n"
            s += str(r) + "\n"
        return s


def readPDB(pdb_file):

    protein = Protein(pdb_file)

    curr_residue = None
    for line in open(pdb_file):
        
        if curr_residue == None:
            curr_residue = Residue(line[17:20], line[22:26].strip())

        # if the record type is an atom, then create new atom
        if line[0:5].strip() == "ATOM":
            print(line[6:11])

            n_atom = Atom(line[6:11], line[12:16], line[30:38], line[38:46], line[46:54]) 

            # If new atom is part of the current residue, add it; else add the 
            # the residue to the protein and start on a new one.
            if line[22:26].strip() == curr_residue.seq_num:
                curr_residue.addAtom(n_atom)

            else:
                protein.addResidue(curr_residue)
                curr_residue = Residue(line[17:20], line[22:26].strip())
                curr_residue.addAtom(n_atom)

    protein.addResidue(curr_residue)

    return protein



       
