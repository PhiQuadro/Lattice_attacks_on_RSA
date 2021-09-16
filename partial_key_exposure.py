from sage.all import *
from my_coppersmith import Coppersmith_univariate
from Crypto.Util.number import getPrime, getRandomNBitInteger
from random import randint
from partial_p_exposure import partial_p

def genkey(nbits, T=65537, lsb=True, guessbits=0):
    p = getPrime(nbits//2)
    q = getPrime(nbits//2)
    N = p*q
    phi = (p-1)*(q-1)
    if lsb:
        e = 2
        while gcd(e,phi) > 1:
            e = randint(2,T)
        d = pow(e,-1,phi)
        d0bits = ceil(nbits/4)
        d0 = d % 2**d0bits
        guess = (p % 2**(d0bits + guessbits)) >> d0bits
        return N,e,d0,d0bits,guess
    else:
        t = randint(ceil(nbits/4)+nbits//32,nbits//2)
        e = getPrime(t)
        d = pow(e,-1,phi)
        d0 = (d >> (d.bit_length()-t)) << (d.bit_length()-t)
        return N,e,d0,t

def recover_roots(f, power):
    res = ['']
    clone = ['']
    p = 2
    for i in range(1,power+1):
        for idx,r in enumerate(res):
            edited = False
            for bit in ['0','1']:
                result = f(int(bit+r,2))
                if result % p**i == 0:
                    if not edited:
                        clone[idx] = bit+r
                        edited = True
                    else:
                        clone.append(bit+r)
        clone = [c for c in clone if len(c) == i]
        res = clone[:]

    return [int(r,2) for r in res]

def lsb_exposure(nbits, T):
    print(f'testing lsb private key exposure on {nbits} bits modulo\n')
    resp = 'y'
    c = 0
    beta = 1/2
    while resp != 'n':
        eps = 1/(32*2**c)
        Xbits = int(nbits*(beta**2 -eps))
        m = ceil((beta**2)/eps)
        t = floor(m*(1/beta - 1))
        guessbits = nbits//4 - Xbits
        print(f'Need to guess {guessbits} bits of p, or lower eps.')
        print(f'Current dimension of the lattice: {m+t}')
        resp = input('Halve eps? (y/n)').lower()
        print()
        c += 1

    N,e,d0,d0bits,guess = genkey(nbits, T=T, lsb=True, guessbits=guessbits)
    print(f'done genkey:\nN = {N}\ne = {e}\nd0 = {d0}')

    R = PolynomialRing(Zmod(2**d0bits),names=('x',))
    (x,) = R._first_ngens(1)
    s = var('s')
    for k in range(1,T):
        if not k%10:
            print(f'k: {k}/{T}')
        possible_s = solve_mod([e*d0 == 1 + k*(N - s +1)], 2**d0bits)
        for s0 in possible_s:
            f = x**2 -int(s0[0])*x + N
            f = f.change_ring(ZZ)
            roots = recover_roots(f,d0bits)
            for r in roots:
                res = partial_p(r+guess*2**d0bits,N,eps=eps,M=2**(d0bits+guessbits),sage_impl=True,v=False)
                if res != []:
                    return res
    return None

def msb_exposure(nbits):
    print(f'testing msb private key exposure on {nbits} bits modulo\n')
    N,e,d0,t = genkey(nbits, lsb=False)
    print(f'done genkey:\nN = {N}\ne = {e}\nd0 = {d0}\nt = {t}')

    k_tilde = (e*d0 -1)//N

    R = PolynomialRing(GF(e),names=('x',))
    (x,) = R._first_ngens(1)

    for i in range(-40,41):
        if not (i+40)%20:
            print(f'i = {i+40}/80')
        k = k_tilde + i
        if k%e == 0:
            continue
        s = ((1 + k*N +k)*inverse_mod(k,e)) %e
        f = x**2 -int(s)*x +N
        roots = f.roots()
        for p0 in roots:
            res = partial_p(int(p0[0]),N,eps=1/32,M=e,sage_impl=True,v=False)
            if res != []:
                return res
    return None

lsb_exposure(nbits=512, T=100)
print('\n---------------------------\n')
msb_exposure(1024)