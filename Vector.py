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

    def __init__(self, l):
        
        # Represent a list as a Vector -- all operations are on a native Python list
        self.vec = l

    def add(self, v2)
        
        assert v2.length() == self.length()

        n_l = []
        for i in self.length(): 
            n_l.append(self.get(i) + v2.get(i))
        return Vector(n_l)

    def subtract(self, v2):
        
        assert v2.length() == self.length()

        n_l = []
        for i in self.length():
            n_l.append(self.get(i) - v2.get(i))
        return Vector(n_l)

    def dot(self, v2):
        
        assert v2.length() == self.length()

        prod = 0
        for i in self.length():
            prod += self.get(i) * v2.get(i)
        return prod

    def cross(self, v2):
        #TODO: IMPLEMENT CROSS PROD CAPABILITY
        
    
    def get(self, i):
        return self.vec[i]

    def length(self):
        return len(self.vec)
