from Crypto.Util.number import getPrime
from random import randint

def genkey(nbits,delta=0.284):
    p,q = 5,5
    while gcd(p-1,q-1) > 2 or p == q:
        p = getPrime(nbits//2)
        q = getPrime(nbits//2)
    N = p*q
    phi = (p-1)*(q-1)
    upperbound = floor(N^delta)
    d = 2
    while gcd(d,phi) > 1: 
        d = randint(3,upperbound)
    e = inverse_mod(d,phi)
    return N,e

def build_lattice(F,m,t):
    G = Sequence([], F.parent())
    x,y = F.variables()
    for k in range(m+1):
        for i in range(m-k+1):
            g_ik = (x^i)*(F^k)*(e^(m-k))
            G.append(g_ik)
        for j in range(1,t+1):
            h_jk = (y^j)*(F^k)*(e^(m-k))
            G.append(h_jk)
    return G

def Boneh_Durfee(N,e,delta=0.284):
    A = (N + 1)//2
    R.<x,y> = PolynomialRing(ZZ)
    Rres.<z> = PolynomialRing(ZZ)
    F = x*(A+y) -1
    Y = int(e^(1/2))
    X = int(2*e^delta)
    
    for m in range(9,10):
        t = ceil(m*(1 - 2*delta)/2)
        print(f'm = {m}, t = {t}')

        G = build_lattice(F,m,t)

        G, monomials = G.coefficient_matrix()

        for idx,monomial in enumerate(monomials):
            for jdx in range(G.nrows()):
                G[jdx,idx] *= monomial(X,Y)[0]

        G = matrix(ZZ,G).dense_matrix().LLL()

        G = copy(G[:2])

        for idx,monomial in enumerate(monomials):
            for jdx in [0,1]:
                G[jdx,idx] //= monomial(X,Y)[0]

        pol0 = G[0]*vector(monomials)
        pol1 = G[1]*vector(monomials)

        resx = pol0.resultant(pol1, x)

        if resx.is_zero() or resx.monomials() == [1]:
            continue
        
        resx = resx(z,z)

        roots = resx.roots()

        for y0 in roots:
            f = z*(-2*int(y0[0]) -z) -N
            for p in f.roots():
                if gcd(p[0], N) > 1:
                    p = p[0]
                    q = N//p
                    print('found factorization:')
                    print(f'p: {p}, q: {q}, p*q == N: {p*q == N}')
                    return p,q
    return None

delta = 0.273
N,e = genkey(nbits=1024,delta=delta)
print('done genkey')
print(f'N = {N}\ne = {e}\ndelta = {delta}')
Boneh_Durfee(N,e,delta=delta)