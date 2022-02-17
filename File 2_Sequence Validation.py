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

Nucleotide = ["A", "C", "G", "T"]


def valid(dhn):
    validseq = dhn.upper()
    for nuc in validseq:
        if nuc not in Nucleotide:
            return false
        return validseq



