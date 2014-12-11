"""Code for reading in different files."""
import re
from RptMatrix import *
from matrix_fun import *

def read_mat(file, target_chr = None):
    """Generator: for returning the contents of a mat file (e.g. hg18_genes.mat) on a given chromosome."""
    wp = open(file)
    for line in wp:
        A = re.split("\s+", line.rstrip())
        if (target_chr and A[1] != target_chr):
            continue
        partNo = int(A[0])
        start = int(A[2])
        finish = int(A[3])
        family = A[5]
        M = sparse2matrix(" ".join(A[8:]))
        
        R = RptMatrix(None, None)
        R.class_name = A[5]
        R.rep_name = A[6]
        R.M = reduce_64matrix(M)
        yield partNo, start, finish, R
