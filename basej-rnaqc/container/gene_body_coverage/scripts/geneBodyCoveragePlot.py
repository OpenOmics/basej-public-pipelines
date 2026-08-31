#!/usr/bin/env python3


#import built-in modules
import os,sys
if sys.version_info[0] != 3:
	print("\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " This verion of RSeQC needs python3!\n", file=sys.stderr)
	sys.exit()	

import re
import string
from optparse import OptionParser
import warnings
import collections
import math
from time import strftime
import subprocess
from os.path import basename
import operator

#import third-party modules
from numpy import std,mean
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
import pysam

#import my own modules
from qcmodule import getBamFiles
from qcmodule import mystat
#changes to the paths

#Isai module
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("df_all.tsv",sep = "\t")
df_float = df.iloc[:,1:].apply(lambda row : pd.Series([float(i) for i in row]),axis = 1)

epsilon = 1e-7  # small constant added to avoid division by zero
df_trans = df_float.iloc[:,1:].apply(lambda row : pd.Series([(i -min(row))/(max(row) - min(row) + epsilon) for i in row]),axis = 1)
df_trans['SampleId'] = df['Percentile']
melted = pd.melt(df_trans,id_vars = "SampleId")
	
sns.lineplot(x="variable", y="value",hue = 'SampleId',legend = False,data=melted).set(xlabel='5\'-3\' body percentile', ylabel='Coverage')
plt.savefig('transcript_body_coverage_mqc.png', dpi=300)