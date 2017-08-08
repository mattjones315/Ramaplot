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

    def add(self, v2):
        
        assert v2.length == self.length

        n_l = []
        for i in range(self.length): 
            n_l.append(self.get(i) + v2.get(i))
        return Vector(n_l)

    def subtract(self, v2):
        
        assert v2.length == self.length

        n_l = []
        for i in range(self.length):
            n_l.append(self.get(i) - v2.get(i))
        return Vector(n_l)

    def dot(self, v2):
        
        assert v2.length == self.length

        prod = 0
        for i in range(self.length):
            prod += self.get(i) * v2.get(i)
        return prod

    def cross(self, v2):
        assert self.length == 3 and v2.length == 3

        a1 = self.get(0)
        a2 = self.get(1)
        a3 = self.get(2)

        b1 = v2.get(0)
        b2 = v2.get(1)
        b3 = v2.get(2)

        n_l = [0]*3
        n_l[0] = a2*b3 - a3*b2
        n_l[1] = a3*b1 - a1*b3
        n_l[2] = a1*b2 - a2*b1

        return Vector(n_l)

    def normalize(self):

        norm_factor = self.squared_sum() ** (1/2)
        norm_vec = Vector([norm_factor]*self.length)

        return self / norm_vec

    def place(self, i , val):
        self.vec[i] = val
        
    def get(self, i):
        return self.vec[i]

    def squared_sum(self):  
        total = 0
        for i in range(self.length):
            total += self.get(i)**2

        return total

    @property
    def length(self):
        return len(self.vec)

    @property
    def mean(self):
        return sum(self.vec) / self.length

    @property
    def sd(self):
        muvec = Vector([self.mean]*self.length)
        
        sdvec = self - muvec
        sdvec2 = sdvec * sdvec
        return (sdvec2.mean) ** (1/2)

    def __str__(self):
        return str(self.vec)

    def __repr__(self):
        return str(self.vec)

    def __add__(self, other):
        return self.add(other)

    def __sub__(self, other):
        return self.subtract(other)

    def __mul__(self, other):
        
        assert type(other) == int or type(other) == float or self.length == other.length 

        n_l = []
        for i in range(self.length):
            if type(other) == int:
                n_l.append(self.get(i) * other)
            else:
                n_l.append(self.get(i) * other.get(i))
        return Vector(n_l)

    def __rmul__(self, other):

        n_l = []
        for i in range(self.length):
            n_l.append(self.get(i) * other)

        return Vector(n_l)


    def __div__(self, other):

        assert self.length == other.length

        n_l = []
        for i in range(self.length):
            n_l.append(self.get(i) / other.get(i))
        return Vector(n_l)

    def __truediv__(self, other):

        assert self.length == other.length

        n_l = []
        for i in range(self.length):
            n_l.append(self.get(i) / other.get(i))
        return Vector(n_l)

    def __floordiv__(self, other):
        
        assert self.length == other.length

        n_l = []
        for i in range(self.length):
            n_l.append(self.get(i) // other.get(i))
        return Vector(n_l)
        
