from sage.all import *
from my_coppersmith import Coppersmith_univariate
from Crypto.Util.number import getPrime, getStrongPrime
from random import choice, randint
import sys

def pad(m, nbits, mbits):
    return 2**nbits -2**mbits +m

def genkey(nbits):
    p,q = getPrime(nbits//2), getPrime(nbits//2)
    N = p*q
    e = choice([3,5])
    return N,e

if __name__ == '__main__':
    eps = 1/8

    if len(sys.argv) == 1:
        #need to keep parameters low due to slow implementation
        nbits = 80

    else:
        #full power fpLLL
        nbits = 2048

    N,e = genkey(nbits)
    mbits = floor(nbits*(1/e - eps))
    m = randint(1,2**mbits)
    c = pow(pad(m,nbits,mbits),e,N)

    R = PolynomialRing(Zmod(N),names=('x',))
    (x,) = R._first_ngens(1)
    F = (pad(x, nbits, mbits))**e -c

    print(f'N: {N}, m: {m}, c: {c}, mbits: {mbits}')

    if nbits > 80:
        sage_impl = True
    else:
        sage_impl = False

    print(Coppersmith_univariate(F,eps=eps,sage_impl=sage_impl))