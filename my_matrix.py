#!/usr/bin/env python3
from my_vector import Vector
from fractions import Fraction

class my_Matrix(list):
    def __init__(self, rows):
        try:
            assert len(set(len(row) for row in rows)) == 1, str(rows) + ' is not rectangular'
            super().__init__(map(Vector, rows))
        except TypeError as te:
            print('cannot cast', rows, 'to list')
            exit(2)
    
    def __repr__(self):
        return '[{}]'.format("\n ".join(Vector.__repr__(row) for row in self))
    
    def Gram_Schmidt(self):
        b = self
        ort_b = []
        mu_mat = []
        for bi in b:
            v = Vector(bi)
            mui = []
            for ort_bj in ort_b:
                projij, muij = bi.proj(ort_bj)
                v = v - projij
                mui.append(muij)
            mu_mat.append(mui + [1] + [0 for _ in range(len(b)-len(ort_b)-1)])
            ort_b.append(v)
        return my_Matrix(ort_b), my_Matrix(mu_mat)
    
    def LLL(self, delta = Fraction(3,4)):
        b = self
        n = len(b)
        ort_b, mu = b.Gram_Schmidt()
        B = [o.n2() for o in ort_b]
        k = 1
        while k < n:
            for j in range(k-1,-1,-1):
                qj = round(mu[k][j])
                b[k] = b[k] - qj*b[j]
                ort_b, mu = b.Gram_Schmidt()
            if B[k] >= (delta - mu[k][k-1]**2)*B[k-1]:
                k += 1
            else:
                b[k], b[k-1] = b[k-1], b[k]
                ort_b, mu = b.Gram_Schmidt()
                B = [o.n2() for o in ort_b]
                k = max(1,k-1)
        return b



