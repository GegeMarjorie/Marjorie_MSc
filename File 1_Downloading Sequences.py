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


"""Annex 8. GmDHN download"""

Entrez.email = "marjoirdi@gmail.com"
GmDHN_Ns = "   351724452.fasta"
if not os.path.isfile(GmDHN_Ns):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  351724452, rettype = "fasta", retmode ="text")
    out_handle_1 = open(GmDHN_Ns, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(GmDHN_Ns, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
GmDHN_1 = SeqIO.read(GmDHN_Ns, "fasta")
GmDHN = GmDHN_1.seq


"""Annex 9. PvDHN download"""

Entrez.email = "marjoirdi@gmail.com"
PvDHN_Ns = "  593267520.fasta"
if not os.path.isfile(PvDHN_Ns):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  593267520, rettype = "fasta", retmode ="text")
    out_handle_1 = open(PvDHN_Ns, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(PvDHN_Ns, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
PvDHN_1 = SeqIO.read(PvDHN_Ns, "fasta")
PvDHN = PvDHN_1.seq


"""Annex 10. SbDHN download"""

Entrez.email = "marjoirdi@gmail.com"
SbDHN_Ns = "  1003097071.fasta"
if not os.path.isfile(SbDHN_Ns):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  1003097071, rettype = "fasta", retmode ="text")
    out_handle_1 = open(SbDHN_Ns, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(SbDHN_Ns, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
SbDHN_1 = SeqIO.read(SbDHN_Ns, "fasta")
SbDHN = SbDHN_1.seq

"""Annex 11. PsDHN download"""

Entrez.email = "marjoirdi@gmail.com"
PsDHN_Ns = "  20704.fasta"
if not os.path.isfile(PsDHN_Ns):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  20704, rettype = "fasta", retmode ="text")
    out_handle_1 = open(PsDHN_Ns, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(PsDHN_Ns, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
PsDHN_1 = SeqIO.read(PsDHN_Ns, "fasta")
PsDHN = PsDHN_1.seq


"""Annex 12. HbDHN download"""

Entrez.email = "marjoirdi@gmail.com"
HbDHN_Ns = "  1031989056.fasta"
if not os.path.isfile(HbDHN_Ns):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  1031989056, rettype = "fasta", retmode ="text")
    out_handle_1 = open(HbDHN_Ns, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(HbDHN_Ns, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
HbDHN_1 = SeqIO.read(HbDHN_Ns, "fasta")
HbDHN = HbDHN_1.seq


"""Annex 13. OsDHN  download"""

Entrez.email = "marjoirdi@gmail.com"
OsDHN_Ns = "  55274277.fasta"
if not os.path.isfile(OsDHN_Ns):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  55274277, rettype = "fasta", retmode ="text")
    out_handle_1 = open(OsDHN_Ns, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(OsDHN_Ns, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
OsDHN_1 = SeqIO.read(OsDHN_Ns, "fasta")
OsDHN = OsDHN_1.seq



"""Annex 14. MsDHN download"""

Entrez.email = "marjoirdi@gmail.com"
MsDHN_Ns = "  336396963.fasta"
if not os.path.isfile(MsDHN_Ns):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  336396963, rettype = "fasta", retmode ="text")
    out_handle_1 = open(MsDHN_Ns, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(MsDHN_Ns, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
MsDHN_1 = SeqIO.read(MsDHN_Ns, "fasta")
MsDHN = MsDHN_1.seq


"""Annex 15. CsDHN download"""

Entrez.email = "marjoirdi@gmail.com"
CsDHN_NS = "  363497941.fasta"
if not os.path.isfile(CsDHN_NS):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  363497941, rettype = "fasta", retmode ="text")
    out_handle_1 = open(CsDHN_NS, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(CsDHN_NS, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
CsDHN_1 = SeqIO.read(CsDHN_NS, "fasta")
CsDHN = CsDHN_1.seq

"""Annex 16. ZmDHN download"""

Entrez.email = "marjoirdi@gmail.com"
ZmDHN_NS = "  162461918.fasta"
if not os.path.isfile(ZmDHN_NS):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  162461918, rettype = "fasta", retmode ="text")
    out_handle_1 = open(ZmDHN_NS, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(ZmDHN_NS, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
ZmDHN_1 = SeqIO.read(ZmDHN_NS, "fasta")
ZmDHN = ZmDHN_1.seq


"""Annex 17. R_LbDHN download"""

Entrez.email = "marjoirdi@gmail.com"
R_LbDHN_NS = "  124703009.fasta"
if not os.path.isfile(R_LbDHN_NS):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  124703009, rettype = "fasta", retmode ="text")
    out_handle_1 = open(R_LbDHN_NS, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(R_LbDHN_NS, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
R_LbDHN_1 = SeqIO.read(R_LbDHN_NS, "fasta")
R_LbDHN = R_LbDHN_1.seq


"""Annex 18. R_CpDHN download"""

Entrez.email = "marjoirdi@gmail.com"
R_CpDHN_NS = "  18083.fasta"
if not os.path.isfile(R_CpDHN_NS):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  18083, rettype = "fasta", retmode ="text")
    out_handle_1 = open(R_CpDHN_NS, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(R_CpDHN_NS, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
R_CpDHN_1 = SeqIO.read(R_CpDHN_NS, "fasta")
R_CpDHN = R_CpDHN_1.seq


"""Annex 19. R_XvDHN download"""

Entrez.email = "marjoirdi@gmail.com"
R_XvDHN_NS = " 30349506.fasta"
if not os.path.isfile(R_XvDHN_NS):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  30349506, rettype = "fasta", retmode ="text")
    out_handle_1 = open(R_XvDHN_NS, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(R_XvDHN_NS, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
R_XvDHN_1 = SeqIO.read(R_XvDHN_NS, "fasta")
R_XvDHN = R_XvDHN_1.seq


"""Annex 20. TaDHN download"""

Entrez.email = "marjoirdi@gmail.com"
TaDHN_NS = " 190684062.fasta"
if not os.path.isfile(TaDHN_NS):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  190684062, rettype = "fasta", retmode ="text")
    out_handle_1 = open(TaDHN_NS, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(TaDHN_NS, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
TaDHN_1 = SeqIO.read(TaDHN_NS, "fasta")
TaDHN = TaDHN_1.seq


"""Annex 21. HvDHN download"""

Entrez.email = "marjoirdi@gmail.com"
HvDHN_NS = " 18965.fasta"
if not os.path.isfile(HvDHN_NS):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  18965, rettype = "fasta", retmode ="text")
    out_handle_1 = open(HvDHN_NS, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(HvDHN_NS, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
HvDHN_1 = SeqIO.read(HvDHN_NS, "fasta")
HvDHN = HvDHN_1.seq


"""Annex 22. OeDHN dowmload"""

Entrez.email = "marjoirdi@gmail.com"
OeDHN_NS = " 1278997624.fasta"
if not os.path.isfile(OeDHN_NS):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id =  1278997624, rettype = "fasta", retmode ="text")
    out_handle_1 = open(OeDHN_NS, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(OeDHN_NS, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
OeDHN_1 = SeqIO.read(OeDHN_NS, "fasta")
OeDHN = OeDHN_1.seq


"""Annex 23. Ta2DHN download"""

Entrez.email = "marjoirdi@gmail.com"
Ta2DHN_NS = " 1399053268.fasta"
if not os.path.isfile(Ta2DHN_NS):
    net_handle_1 = Entrez.efetch(db = "nucleotide", id = 1399053268, rettype = "fasta", retmode ="text")
    out_handle_1 = open(Ta2DHN_NS, "w")
    out_handle_1.write(net_handle_1.read())
    out_handle_1.close()
    print()
    print()
    record_1 = SeqIO.read(Ta2DHN_NS, "fasta")
    print() # There is : id, name, description, record, seq, all in record_1
    print()
Ta2DHN_1 = SeqIO.read(Ta2DHN_NS, "fasta")
Ta2DHN = Ta2DHN_1.seq



