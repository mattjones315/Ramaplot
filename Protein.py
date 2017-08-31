"""
Contains functions and objects related to representing proteins in this Ramaplot module.
"""

from __future__ import absolute_import, print_function, division
import argparse
import os
import numpy as np
import CLI
import Vector
import math


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
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def coord(self):
        """
        Return coordinates of atom as a Vector.
        """

        return Vector.Vector([self.x, self.y, self.z])

    def __str__(self):
        rep = self.serial_number + " " + self.name + \
                " (" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ")"
        return rep

class Residue:

    def __init__(self, name, sequence_number):
        self.name = name
        self.seq_num = sequence_number

        self.atoms = []

        # Dihedral angles
        self.omega = 0
        self.phi = 0
        self.psi = 0

    def addAtom(self, atom):
        self.atoms.append(atom)

    def get(self, i):

        assert i >= 0 and i < len(self.atoms)

        return self.atoms[i]

    def compute_dihedral(self, _type,  atoms):
        
        a1 = atoms[0]
        a2 = atoms[1]
        a3 = atoms[2]
        a4 = atoms[3]

        b1 = a2 - a1
        b2 = a2 - a3
        b3 = a4 - a3

        b1 = b1.normalize()
        b2 = b2.normalize()
        b3 = b3.normalize()

        n1 = b1.cross(b2)
        n2 = b2.cross(b3)

        m1 = n1.cross(b2)

        x = n1.dot(n2)
        y = m1.dot(n2)

        dihedral = math.degrees(math.atan2(y, x))

        if _type == "omega":
            self.omega = dihedral
        elif _type == "phi":
            self.phi = dihedral
        elif _type == "psi":
            self.psi = dihedral
        else:
            raise Exception("Dihedral type not recognized")

    @property
    def length(self):
        return len(self.atoms)


    def __str__(self):
        s = ""
        for a in self.atoms:
            s += str(a) + "\n"
        return s

class Protein:

    def __init__(self, name):
        self.name = name
        self.amino_acids = []
        self.omega = {}
        self.phi = {}
        self.psi = {}
    
    def addResidue(self, aa):
        self.amino_acids.append(aa)

    def gather_dihedral_groups(self):
        """
        Gathers groups of atoms that the dihedral angles will be computed from.

        Populates the dihedral group dictionary by appending tuples of 4 coordinates to the
        appropriate group (omega, phi, psi). 
        """

        for a_ii in range(len(self.amino_acids)):
            curr_aa = self.amino_acids[a_ii]

            # If this is the first residue, can't calculate omega or phi
            if curr_aa.seq_num == 0:

                # Also make sure it's not a one residue protein 
                if curr_aa.seq_num < len(self.amino_acids): 
                    next_aa = self.amino_acids[a_ii + 1]
                    self.psi[curr_aa.seq_num] = [curr_aa.get(0).coord(), curr_aa.get(1).coord(),
                                            curr_aa.get(2).coord(), next_aa.get(0).coord()] 

            else:

                # If this is the last residue, can't calculate psi
                if a_ii == (len(self.amino_acids)-1):
                    prev_aa = self.amino_acids[a_ii - 1]
                    self.phi[curr_aa.seq_num] = [prev_aa.get(2).coord(), curr_aa.get(0).coord(), 
                                                curr_aa.get(1).coord(), curr_aa.get(2).coord()]
                    self.omega[curr_aa.seq_num] =[prev_aa.get(1).coord(), prev_aa.get(2).coord(),
                                                curr_aa.get(0).coord(), curr_aa.get(1).coord()]
                            

                else:
                    prev_aa = self.amino_acids[a_ii - 1]
                    next_aa = self.amino_acids[a_ii + 1]

                    self.omega[curr_aa.seq_num] = [prev_aa.get(1).coord(), 
                        prev_aa.get(2).coord(), curr_aa.get(0).coord(), curr_aa.get(1).coord()]
                    self.phi[curr_aa.seq_num] = [prev_aa.get(2).coord(), curr_aa.get(0).coord(), 
                                                curr_aa.get(1).coord(), curr_aa.get(2).coord()]
                    self.psi[curr_aa.seq_num] = [curr_aa.get(0).coord(), curr_aa.get(1).coord(),
                                            curr_aa.get(2).coord(), next_aa.get(0).coord()] 



    def compute_dihedrals(self):
        """
        Computes omega, phi, and psi angles for each residue in the given protein.
        """

        for a_ii in range(len(self.amino_acids)):
            curr_aa = self.amino_acids[a_ii]

            if curr_aa.seq_num in self.omega:
                curr_aa.compute_dihedral("omega", self.omega[curr_aa.seq_num])
            if curr_aa.seq_num in self.phi: 
                curr_aa.compute_dihedral("phi", self.phi[curr_aa.seq_num])
            if curr_aa.seq_num in self.psi:
                curr_aa.compute_dihedral("psi", self.psi[curr_aa.seq_num])


    def get_psi_angles(self):
        
        psi_a = []
        for a_ii in range(len(self.amino_acids)):
            psi_a.append(self.amino_acids[a_ii].psi)
        return psi_a

    def get_phi_angles(self):

        phi_a = []
        for a_ii in range(len(self.amino_acids)):
            phi_a.append(self.amino_acids[a_ii].phi)
        return phi_a

    def get_omega_angles(self):

        omega_a = []
        for a_ii in range(len(self.amino_acids)):
            omega_a.append(self.amino_acids[a_ii].omega)
        return omega_a

    @property
    def length(self):
        return len(self.amino_acids)


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

            n_atom = Atom(line[6:11].strip(), line[12:16].strip(), line[30:38].strip(), line[38:46].strip(), line[46:54].strip()) 

            # If new atom is part of the current residue, add it; else add the 
            # the residue to the protein and start on a new one.
            if line[22:26].strip() == curr_residue.seq_num:
                curr_residue.addAtom(n_atom)

            else:
                protein.addResidue(curr_residue)
                curr_residue = Residue(line[17:20], line[22:26].strip())
                curr_residue.addAtom(n_atom)

    if curr_residue.length > 0:
        protein.addResidue(curr_residue)

    return protein



       
