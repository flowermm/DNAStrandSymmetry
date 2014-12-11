"""Provides various functons for working with substitution matricies."""
import re
from numpy import *
from fp_check import *

def sparse2matrix(s, dtype = int):
    """Convert a square sparse-matrix string representation to a numpy matrix object.
    String format: matrix dim, # non-zero elements, (index,value)*, "*"
    Indexing assumes row-major ordering."""
    A = re.split("\s+", s.rstrip())
    n = int(A[0])
    M = zeros((n,n))
    for i in range(2, len(A)-1, 2):
        index = int(A[i])
        M[ index // n, index % n ] = dtype(A[i+1])

    return M

def matrix2sparse(M, dtype = int):
    """Convert a matrix to a sparse representation"""
    s = shape(M)
    assert len(s) == 2 and s[0]==s[1]

    L = ["%d %s" % (i,dtype(v)) for i,v in enumerate(M.flat) if not isZero(v)]

    return "{n} {L} {s} *".format(n=s[0], L=len(L), s = " ".join(L))

def reduce_64matrix(M64):
    """Reduce a 64x64 base-triplet substitution matrix to a 4x4 substitution matrix.
    (Note: this will lose the two end-point substitutions.)"""
    M4 = zeros((4,4), dtype = M64.dtype)
    
    for i in range(64):
        b = ((i & 12)  >> 2)
        for j in range(64):
            M4[b, (j & 12) >> 2] += M64[i,j]

    return M4
