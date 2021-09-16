from sage.all import *
from my_coppersmith import Coppersmith_univariate
from Crypto.Util.number import getPrime, getRandomNBitInteger
from random import choice, randint
import sys

def partial_p(p0, N, eps=None, known_bits=None, M=None, sage_impl=False, v=True):
    nbits = int(N).bit_length()
    beta = 0.5
    R = PolynomialRing(Zmod(N),names=('x',))
    (x,) = R._first_ngens(1)
    if M:
        F = x*M + p0
        if v:
            print('\n\nlsb case, F(x) =', F)
    else:
        rootbits = nbits//2 - known_bits
        F = (2**(rootbits))*p0 + x
        if v:
            print('\n\nmsb case, F(x) =', F)
    return Coppersmith_univariate(F,beta,eps=eps,sage_impl=sage_impl,verbose=v)
    
def genkey(nbits):
    N = 1
    while int.bit_length(N) != nbits:
        p,q = getPrime(nbits//2), getPrime(nbits//2)
        N = p*q
    
    #need more than half bits of p, cause eps=beta/8
    k = ceil(nbits/16)

    known_bits = nbits//4 + k

    print(f'N = {N} #{int.bit_length(N)}')
    print(f'p = {p} #{int.bit_length(p)}')
    print(f'q = {q} #{int.bit_length(q)}')

    p_msb = p >> (nbits//4 - k)
    target_msb = p % 2**(nbits//4 - k)

    p_lsb = p % 2**(known_bits)
    target_lsb = p >> (known_bits)

    print('testing partial p exposure')
    print(f'target_msb: {target_msb} #{int.bit_length(target_msb)}')
    print(f'target_lsb: {target_lsb} #{int.bit_length(target_lsb)}')

    return N, p_msb, p_lsb, known_bits

if __name__ == '__main__':
    if len(sys.argv) == 1:
        #low parameters, slow implementation
        N, p_msb, p_lsb, known_bits = genkey(80)
        partial_p(p_msb, N, known_bits=known_bits)
        partial_p(p_lsb, N, M = 2**(known_bits))

    else:
        #full power fpLLL
        N, p_msb, p_lsb, known_bits = genkey(2048)
        partial_p(p_msb, N, known_bits, sage_impl=True)
        partial_p(p_lsb, N, M = 2**(known_bits), sage_impl=True)