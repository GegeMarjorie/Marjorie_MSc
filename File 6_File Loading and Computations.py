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


"""
Annex 29. Reading the MEGA X aligned sequences into python. 
MSA file-> dhn.fas
"""

d= open(r'C:\Users\USER\Desktop\thesis docs uy1\Py_DHN_analysis\DHN Aln\dhn.fas', 'r')
a= d.read()
print(a)

"""
Annex 30Annex 30. Reading distance matrix excel files
"""

df = pd.read_excel(r'C:\Users\USER\Desktop\thesis docs uy1\Py_DHN_analysis\DHN Aln\csvdatam\csvdata.xlsx', index_col = 0, header = 0)
df1 = pd.read_excel(r'C:\Users\USER\Desktop\thesis docs uy1\Py_DHN_analysis\DHN Aln\csvdatam\csvdata2.xlsx', index_col = 0, header = 0)
df3 = pd.read_excel(r'C:\Users\USER\Desktop\thesis docs uy1\Py_DHN_analysis\DHN Aln\csvdatam\csvdata3.xlsx', index_col= 0, header = 0)
print(df3)


"""
Annex 31. Heatmap generation from distance matrix data
"""

#(csvdata3=df3)
hm = df3.plot(cmap='Reds', title = 'DHN Heatmap')
plt.show()

"""
Annex 32. Scatter plot generation from distance matrix
"""

#(df3)
sn = pd.plotting.scatter_matrix(df3, alpha = 0.2, figsize = (12,10), diagonal = 'kde')
plt.show()

"""
Annex 33. Phylo Tree construction with heatmap.
"""

hm= sns.clustermap(df3)
plt.show()


"""
Annex 34. Heatmaps and scatter plot of amino acid 
frequencies of the dehydrin genes
"""

asn = pd.plotting.scatter_matrix(aa, alpha = 0.2, figsize = (12,10))
plt.show()
# AAHM
ahm = sns.heatmap(aa)
plt.show()




