"""
Correlation between sequence Lengths,
AT/GC percentages of the dehydrins.

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


x = np.array(['GmDHN', 'PvDHN', 'SbDHN', 'PsDHN', 'HbDHN', 'OsDHN', 'MsDHN', 'CsDHN', 'ZmDHN', 'R_LbDHN', 'R_CpDHN', 'R_XvDHN', 'TaDHN', 'HvDHN', 'OeDHN', 'Ta2DHN'])
slen = np.array([730, 1079, 459, 785, 675, 1117, 711, 1027, 862, 672, 2330,614, 987, 683, 843, 905])
gcc = np.array([ 46.30, 40.96, 66.45, 40.13, 46.67, 59.27, 57.10, 46.15, 59.51, 46.58, 40.04, 45.77, 60.28, 60.03, 44.13, 58.01])
atc = np.array([ 53.70, 59.04, 33.55, 59.87, 53.33, 40.73, 42.9, 53.85, 40.49, 53.42, 59.96, 54.23, 39.72, 39.97, 55.87, 41.99])
corr1_coef = np.corrcoef(slen,gcc)
corr2_coef = np.corrcoef(slen, atc)
print(corr1_coef)
print(corr2_coef)
# getting regression line
# plot of regression line

def best_fit_slope_and_intercept(xs, ys):
    m = (((mean(xs) * mean(ys)) - mean(xs * ys)) /
         ((mean(xs) * mean(xs)) - mean(xs * xs)))

    b = mean(ys) - m * mean(xs)

    return (m,b)

print(best_fit_slope_and_intercept(slen, gcc))
print(best_fit_slope_and_intercept(slen, atc))

m1 = -0.007089147702876768
b1 = 57.4948395234006
m2 = 0.007648221697438096
b2 = 41.999757585515965

fig, ax = plt.subplots()
ax.plot(slen,gcc, linewidth=0)
ax.plot(slen, m1 * slen + b1, label='Regression_line_gc')
ax.set_xlabel('Sequence Lengths')
ax.set_ylabel('GC percentages')
ax.set_title('GC% Vs Sequence lengths')
ax.legend()
plt.scatter(slen,gcc)
plt.show()

fig, ax = plt.subplots()
ax.plot(slen,atc, linewidth=0)
ax.plot(slen, m2 * slen + b2 , label = 'Regression_line_at')
ax.set_xlabel('Sequence Lengths')
ax.set_ylabel('AT Percentages')
ax.set_title('AC% Vs Sequence lengths')
ax.legend()
plt.scatter(slen, atc)
plt.show()

fig, ax = plt.subplots()
ax.plot(slen,gcc, linewidth=0)
ax.plot(slen, m1 * slen + b1, label='Regression_line_gc')
ax.plot(slen,atc, linewidth=0)
ax.plot(slen, m2 * slen + b2 , label = 'Regression_line_at')
ax.set_xlabel('Sequence Lengths')
ax.set_ylabel('GC/AT Percentages')
ax.set_title(' cummulative graph of GC%/AT% Vs Sequence lengths')
ax.legend()
plt.show()

