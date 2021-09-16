#!/usr/bin/env python3
from fractions import Fraction

class Vector(list):
    def __init__(self, x):
        try:
            super().__init__(map(Fraction, x))
        except TypeError as te:
            print(x, 'is not iterable or its elements cannot be cast to Fractions')
            exit(2)
    
    def __add__(self, other):
        if not isinstance(other, Vector):
            other = Vector(other)
        assert len(self) == len(other), 'incompatible dimensions'
        return Vector(x+y for x,y in zip(self,other))
    
    def __sub__(self, other):
        if not isinstance(other, Vector):
            other = Vector(other)
        assert len(self) == len(other), 'incompatible dimensions'
        return Vector(x-y for x,y in zip(self,other))
    
    def __mul__(self, scal):
        if not isinstance(scal, Fraction):
            scal = Fraction(scal)
        return Vector(x*scal for x in self)
    
    def __rmul__(self, scal):
        if not isinstance(scal, Fraction):
            scal = Fraction(scal)
        return Vector(x*scal for x in self)
    
    def __truediv__(self, scal):
        if not isinstance(scal, Fraction):
            scal = Fraction(scal)
        return Vector(x/scal for x in self)
    
    def __repr__(self):
        return f'[{", ".join(str(x) for x in self)}]'
    
    def dot(self, other):
        if not isinstance(other, Vector):
            other = Vector(other)
        assert len(self) == len(other), 'incompatible dimensions'
        return sum(map(lambda x: x[0]*x[1], zip(self, other)))

    def n2(self):
        return self.dot(self)
    
    def mu(self, j):
        if not isinstance(j, Vector):
            j = Vector(j)
        assert len(self) == len(j), 'incompatible dimensions'
        return self.dot(j)/self.n2()
    
    def proj(self, other):
        if not isinstance(other, Vector):
            other = Vector(other)
        assert len(self) == len(other), 'incompatible dimensions'
        muij = other.mu(self)
        return muij*other, muij