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


# Plots. Bar plots of DHN and sequence lengths, GC % and AT %

w = 0.4

y= ['GmDHN', 'PvDHN', 'SbDHN', 'PsDHN', 'HbDHN', 'OsDHN', 'MsDHN', 'CsDHN', 'ZmDHN', 'R_LbDHN', 'R_CpDHN', 'R_XvDHN', 'TaDHN', 'HvDHN', 'OeDHN', 'Ta2DHN']
sL = [7.30, 10.79, 4.59, 7.85, 6.75, 11.17, 7.11, 10.27, 8.62, 6.72, 23.30, 6.14, 9.87, 6.83, 8.43, 9.05]
gC = [46.30, 40.96, 66.45, 40.67, 46.67, 59.27, 57.10, 46.15, 59.51, 46.58, 40.04, 45.77, 60.28, 60.03, 44.13, 58.01]
aT= [53.70, 59.04, 33.55, 59.87, 53.33, 40.73, 42.90, 53.85, 40.49, 53.42, 59.96, 54.23, 39.72, 39.97, 55.87, 41.99]

plt.bar(sL,y,w,label= "Sequence lengths")
plt.xlabel('Sequence Lengths')
plt.ylabel('Dehydrins')
plt.title('Bar graph showing sequence lengths of dehydrin genes')
plt.legend()
plt.show()

plt.bar(gC,y,w, label="GC percentages")
plt.xlabel('GC percentages')
plt.ylabel('Dehydrins')
plt.title('Bar graph showing GC percentages of dehydrin genes')
plt.legend()
plt.show()

plt.bar(aT,y,w, label="AT percentages")
plt.xlabel('AT percentages')
plt.ylabel('Dehydrins')
plt.title('Bar graph showing AT percentages of dehydrin genes')
plt.legend()
plt.show()

plt.bar(sL,y,w,label= "Sequence lengths")
plt.bar(gC,y,w, label="GC percentages")
plt.bar(aT,y,w, label="AT percentages")
plt.xlabel('Units')
plt.ylabel('Dehydrin genes')
plt.title('Bar graph showing cummulative nucleotide informatons of dehydrin genes')
plt.legend()
plt.show()

bar1=np.arange(len(y))
bar2=[i+w for i in bar1]
bar3=[i+w for i in bar2]

plt.bar(bar1,sL,w,label="Sequence lengths")
plt.bar(bar2,gC,w, label="GC percentages")
plt.bar(bar3,aT, w, label=" AT percentages")
plt.xlabel("Dehydrins")
plt.ylabel("Units")
plt.title("Cummulative bar graph showing dehydrin nucleotides information")
plt.xticks(bar1+w,y)
plt.legend()
plt.show()




