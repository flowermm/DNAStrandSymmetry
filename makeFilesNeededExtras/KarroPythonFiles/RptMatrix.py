"""rpt_matrix.py: A RptMatrix class for storing needed information about individual repeats along with 
their C matrix.  Will store these in psm files (pickled substitution matrix)"""

from RMfileReader import *
import argparse
import pickle
import glob
import re
from numpy import zeros

base2index = {'A':0,'C':1,'G':2,'T':3}                 # Dictionary mapping bases to their index position.
baseComplement = {x:y for x,y in zip("ACGT","TGCA")}   # Dictionary mapping bases to their complements. 
baseSet = {'A','C','G','T'}                            # Set of legal bases
def reverse_complement(s):
    return "".join([baseComplement[c] if c in baseComplement else c for c in s[::-1]])

class RptMatrix:
    def __init__(self, r, M):
        """Class for holding both a repeat object and its calculated substitutinon matrix"""
        if r is None:   # Allows us to create an "empty" object"
            assert M is None   # Empty objects shouldn't be given a matrix
            self.ancestor_coords = (None, None)
            self.start, self.finish = None, None
            self.chr = None
            self.strand = None
            self.rep_name = None
            self.class_name = None
            self.M = zeros((4,4))
        else:
            self.ancestor_coords = r.ancestor_coords()
            self.start, self.finish = r.modern_coords()    # Coordinates
            self.chr = r.modern_chr()
            self.strand = r.modern_strand()
            self.rep_name = r.rep_name()
            self.class_name = r.class_name()        
            self.M = M

    def __getitem__(self, index):
        return self.M[index[0], index[1]]

def align2matrix(s1, s2, maskCpG = False):
    """Compute a matrix from aligned sequences"""
    M = zeros((4,4))

    last_c = False
    for c1,c2 in zip(s1, s2):
        if {c1,c2} <= baseSet:
            if (not maskCpG) or last_c == False or c1 != 'G':
                M[base2index[c1],base2index[c2]] += 1
            last_c == (c1 == 'C')

    return M

def Rpt2Matrix(rptObj, maskCpG = False):
    """repeat object -> RptMat object.
    Count set relative to positive strand."""

    M = zeros((4,4))
    anc_seq = rptObj.ancestor_seq() if rptObj.modern_strand() == '+' else reverse_complement(rptObj.ancestor_seq())
    mod_seq = rptObj.modern_seq() if rptObj.modern_strand() == '+' else reverse_complement(rptObj.modern_seq())


    return RptMatrix(rptObj, align2matrix(anc_seq, mod_seq, maskCpG))

def prmList2psm(prm_list, filter_file = None, filter_fun = lambda r: r.M.sum() >= 40):
    """Create a psm list from a prm list"""
    S = {line.strip() for line in open(filter_file)} if filter_file else {}
    R = [Rpt2Matrix(o) for o in prm_list if ((not S) or (o.class_name() in S))]
    if not filter_fun is None:
        R = [r for r in R if filter_fun(r)]
    return R



def prm2psm(prm_file, psm_file, filter_file = None, filter_fun = lambda r: r.M.sum() >= 40):
    """Create a psm file from a prm file"""
    L = unpickleRepeats(prm_file)
    print("L: ", len(L))
    S = {line.strip() for line in open(filter_file)} if filter_file else {}
    print("S: ", len(S))
    R = [Rpt2Matrix(o) for o in L if ((not S) or (o.class_name() in S))]
    if not filter_fun is None:
        R = [r for r in R if filter_fun(r)]
    print("R: ", len(R))
    pickle.dump(R, open(psm_file, "wb"))

def load_psm(file):
    """Load a psm file and return a list of the repeats"""
    return pickle.load(open(file, "rb"))

def rptsInPartition(rpt_list, start, finish):
    """Return a sublist of all repeats contained within the (start, finish) interval"""
    i = 0
    while i < len(rpt_list) and rpt_list[i].start < start:
        i += 1
    j = i+1
    while j < len(rpt_list) and rpt_list[j].finish <= finish:
        j = j+1
    return rpt_list[i:j]    

def create_psm(prm_location = "/Users/karroje/cache/human/hg18/seq/rmsk", filter_file = "martin_repeats.txt"):
    """Create all psm files locally.  By default: use location on Karro's computer."""
    for file in glob.glob(prm_location + "/*.prm"):
        target = re.sub(".prm$", ".psm", file[file.rfind("/")+1:])
        print("File: ", file, target)
        prm2psm(file, psm_file = target, filter_file = filter_file)



def UCSC2tuple(ucsc_gene):
    """Take gene location description in UCSC format and return it as a stting/int/int tuple"""
    i = ucsc_gene.find(":")
    j = ucsc_gene.find("-", i+1)
    return ucsc_gene[:i], int(ucsc_gene[i+1:j]), int(ucsc_gene[j+1:])

def tuple2UCSC(chr, start, finish):
    """Turn a gene-tuple into UCSC format"""
    return "%s:%d-%d" % (chr,start,finish)

def RptListGenerator(gene_list, family_set = None, base_dir = ".", gene_parser = UCSC2tuple, key_format = tuple2UCSC, limiting_chr = None):
    """Take a list of gene intervals and return generate lists of contained RptMatrix objects for each gene.
    * gene_list: List of genes.
    * base_dir: The directory containing the .psm files.
    * gene_parser: Parse the gene elements to a chr, start, finish tuple.
    * key_format: takes a chromosome, start, and finsish parameter and turns it into a gene key.
    * limiting_chr: A set of chrosomes to be used.  (All used if None or {}.)
    NOTE: If there are overlapping genes, only the first will be used.
    """
    gene_list = sorted([gene_parser(gene) for gene in gene_list])

    if family_set is None:
        family_set = {f.rstrip() for f in open("martin_repeats.txt")}

    current_chr = None
    for chr, start, finish in gene_list:
        if limiting_chr and chr not in limiting_chr:
            continue
        if chr != current_chr:
            print("Loading: %s" % (chr))
            rpt_objects = load_psm("%s/%s.psm" % (base_dir, chr))
            current_chr = chr
            rpt_start = 0
        elif start < rpt_objects[rpt_start]:
            continue

        while rpt_start < len(rpt_objects) and rpt_objects[rpt_start].start < start: rpt_start +=1
        rpt_finish = rpt_start
        while rpt_finish < len(rpt_objects) and rpt_objects[rpt_finish].finish < finish: rpt_finish += 1
        #yield key_format(chr, start, finish), rpt_objects, rpt_start, rpt_finish, family_set
        yield key_format(chr, start, finish), [r for r in rpt_objects[rpt_start:rpt_finish] if (r.class_name in family_set) and (r.M.sum() >= 40)]
        rpt_start = rpt_finish
        
        
            
        
def genes2dic(gene_list, family_set = None, base_dir = ".", gene_parser = UCSC2tuple, key_format = tuple2UCSC, limiting_chr = None):
    """Take a list of gene intervals and return a dictionary mapping genes to contained RptMatrix objects.
    * gene_list: List of genes.
    * base_dir: The directory containing the .psm files.
    * gene_parser: Parse the gene elements to a chr, start, finish tuple.
    * key_format: takes a chromosome, start, and finsish parameter and turns it into a gene key.
    * limiting_chr: A set of chrosomes to be used.  (All used if None or {}.)
    NOTE: If there are overlapping genes, only the first will be used.
    """
    return {key:L for key,L in RptListGenerator(gene_list, family_set, base_dir, gene_parser, key_format, limiting_chr)}


def save_geneDic(gene_list, output_file, base_dir = ".", gene_parser = UCSC2tuple, key_format = tuple2UCSC, limiting_chr = None):
    D = genes2dic(gene_list, family_set = None, base_dir = base_dir, gene_parser = gene_parser, key_format = key_format, limiting_chr = limiting_chr)
    pickle.dump(D, open(output_file, "wb"))

    
def load_geneDic(file):
    return pickle.load(open(file))
