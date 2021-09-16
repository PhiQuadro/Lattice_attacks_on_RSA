#!/usr/bin/env sage
from my_matrix import my_Matrix

class Not_a_Basis(Exception):
    pass

def Coppersmith_univariate(F, beta=1, eps=None, sage_impl=False, verbose=True):
    try:
        F.nvariables()
        print('polynomial is multivariate')
        return None
    except AttributeError:
        pass
    
    N = F.base_ring().cardinality()
    
    if eps == None:
        eps = beta/8
    
    d = F.degree()

    try:
        F /= F.coefficients().pop(-1)
    except ZeroDivisionError:
        p = gcd(N, F.coefficients().pop(-1))
        q = N//int(p)
        print('found factorization by ZeroDivisionError:')
        print(f'p: {p}, q: {q}, p*q == N: {p*q == N}')
        return p,q
    
    F = F.change_ring(ZZ)
    x = F.variables()[0]

    exp = RR(-eps + (beta^2)/d).as_integer_ratio()
    exp = exp[0]/exp[1]

    X = int((N^(exp))/2)
    m = ceil((beta^2)/(d*eps))
    t = floor(d*m*(1/beta - 1))

    exp = RR((beta^2)/d).as_integer_ratio()
    exp = exp[0]/exp[1]
    maxB = int(N^exp)
    if verbose:
        print('\nApplying Coppersmith\'s univariate method with parameters:')
        print(f'X: {X}, X_bits: {int.bit_length(int(X))}')
        print('m: {m}, t: {t}, max theoretical bound: {maxB}')
        print(f'F(x): {F}')
    G = []
    n = d*m + t
    for j in range(d):
        for i in range(m):
            g = (x^j)*(N^(m-i)*(F^(i)))
            row = g.coefficients(sparse=False) + [0]*(n-g.degree()-1)
            G.append([int(r*X^k) if r else int(0) for r,k in zip(row,range(n))])
    for i in range(t):
        h = (x^i)*(F^m)
        row = h.coefficients(sparse=False) + [0]*(n-h.degree()-1)
        G.append([int(r*(X^k)) if r else int(0) for r,k in zip(row,range(n))])
    
    if matrix(ZZ,G).rank() != len(G):
        raise Not_a_Basis()
    
    if verbose:
        print(f'lattice dimension: {len(G)}\n')

    if not sage_impl:
        G = my_Matrix(G).LLL()
    else:
        G = matrix(ZZ,G).LLL()

    f = sum([(int(g)//(X^k))*(x^k) if g else 0 for g,k in zip(G[0],range(n))])

    roots = f.roots()
    for r in roots:
        if r[0] > X:
            print('root > bound!', r[0])
        if beta < 1:
            g = gcd(F(r[0]),N)
            if g > 1 and g != N:
                p = g
                q = N//p
                print('found factorization:')
                print(f'p: {p}, q: {q}, p*q == N: {p*q == N}')
                return p,q

    if beta < 1:
        if verbose:
            print('factorization not found!')
        return roots
    else:
        return roots


    
