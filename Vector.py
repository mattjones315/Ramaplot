"""
Representation of a vector, for vector operations. Naturally, using Numpy would be a more 
traditional and efficient way to deal with vector math, but this is part of the assignment.
"""



from __future__ import absolute_import, print_function, division
import argparse
import os
import numpy as np
import CLI
import Protein

class Vector:

    def __init__(l):
        
        # Represent a list as a Vector -- all operations are on a native Python list
        self.vec = l

