"""
Main entry point for the Ramaplot package -- will delegate the following tasks:

    - read in pdb file and store as a Protein object
    - calculate the dihedral angles for all amino acids
    - generate Ramachandran plots from the resulting dihedral angles

Future update: read in multiple pdb files and create multiple plots in one run

Ambitious udpate: create a model for predicting unfavorable conmforations of proteins given
only a sequence of amino acids. 

"""

from __future__ import absolute_import, print_function, division
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import CLI
import Protein
import Vector



def main():
    
    # Read in arguments from command line 
    args = CLI.parseRPArgs()

    # Create protein data object to compute dihedral angles on 
    protein = Protein.readPDB(args["pdb_file"])
    
    # Compute dihedral groups, saved as a dictionary object & compute angles
    protein.gather_dihedral_groups()
    protein.compute_dihedrals()


    # Now plot results
    psi_angles = np.array(protein.get_psi_angles())
    phi_angles = np.array(protein.get_phi_angles())
    
    plt.axis((-180, 180, -180, 180))
    plt.plot(phi_angles, psi_angles, 'ro')
    plt.xlabel("Phi")
    plt.ylabel("Psi")
    
    plt.show()





if __name__ == "__main__":
    main();



