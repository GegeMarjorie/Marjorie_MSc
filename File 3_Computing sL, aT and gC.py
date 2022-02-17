"""
Annex 25. Comparative Sequence Analysis between
the sequence lengths(sL), AT-percentages(aT)  and
GC-percentages (gC) and plotting bar plots Using
python codes.(Aziz, Alhadidi and Mohammed, 2017;
Mallawaarachchi, 2017; Menon, 2018;
Seaman and Buggs, 2020;
Wu et al., 2020, 2020)
"""

from Bio import Entrez
from Bio import SeqIO
from Bio import pairwise2
import os
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
from scipy.cluster.hierarchy import linkage
# import dendrogram
import scipy
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.pairwise2 import format_alignment
from Bio.SeqUtils import seq3
from Bio.SeqUtils import GC
from io import StringIO
import itertools
# import itertools_s
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import _Matrix
import argparse
from statistics import mean
# import dendropy
# import pandas_flavor
# import xlrd
import openpyxl
import seaborn as sns

#sequence lengths
def seq_len(seq):
    length = len(seq)
    return length

# AT% computation
def at(seq):
    count = [N for N in seq if N in ("AT")]
    result = float(round(len(count)/len(seq)*100, 2))

    return result
# GC% computation
def gc(seq):
    count = [N for N in seq if N in 'GC']
    result = float(round(len(count)/len(seq)*100, 2))

    return result




