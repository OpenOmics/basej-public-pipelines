#!/usr/bin/env python3
'''
Calculate the RNA-seq reads coverage over gene body. 
Note:
1) Only input sorted and indexed BAM file(s). SAM format is not supported.
2) Genes/transcripts with mRNA length < 100 will be skipped (Number specified to "-l" cannot be < 100). 
'''

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

#changing history to this module


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="4.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def valid_name(s):
	'''make sure the string 's' is valid name for R variable'''
	symbols = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.'
	digit = '0123456789'
	rid = '_'.join(i for i in s.split())	#replace space(s) with '_'
	if rid[0] in digit:rid = 'V' + rid
	tmp = ''
	for i in rid:
		if i in symbols:
			tmp = tmp + i
		else:
			tmp = tmp + '_'
	return tmp
	

def printlog (mesg):
	'''print progress into stderr and log file'''
	mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
	LOG=open('log.txt','a')
	print(mesg, file=sys.stderr)
	print(mesg, file=LOG)


def pearson_moment_coefficient(lst):
	'''measure skewness'''
	mid_value = lst[int(len(lst)/2)]
	sigma = std(lst, ddof=1)
	tmp = []
	for i in lst:
		tmp.append(((i - mid_value)/sigma)**3)
	return mean(tmp)

def genebody_percentile(refbed, mRNA_len_cut = 100):
	'''
	return percentile points of gene body
	mRNA length < mRNA_len_cut will be skipped
	'''
	if refbed is None:
		print("You must specify a bed file representing gene model\n", file=sys.stderr)
		exit(0)
	
	g_percentiles = {}
	transcript_count = 0
	for line in open(refbed,'r'):
		try:
			if line.startswith(('#','track','browser')):continue  
			# Parse fields from gene tabls
			fields = line.split()
			chrom     = fields[0]
			tx_start  = int( fields[1] )
			tx_end    = int( fields[2] )
			geneName      = fields[3]
			strand    = fields[5]
			geneID = '_'.join([str(j) for j in (chrom, tx_start, tx_end, geneName, strand)])
				
			exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
			exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
			exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
			exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends))
			transcript_count += 1
		except:
			print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
			continue
		gene_all_base=[]
		mRNA_len =0
		flag=0
		for st,end in zip(exon_starts,exon_ends):
			gene_all_base.extend(list(range(st+1,end+1)))		#1-based coordinates on genome
		if len(gene_all_base) < mRNA_len_cut:
			continue
		g_percentiles[geneID] = (chrom, strand, mystat.percentile_list (gene_all_base))	#get 100 points from each gene's coordinates
	printlog("Total " + str(transcript_count) + ' transcripts loaded')
	return g_percentiles

def genebody_coverage(bam, position_list):
	'''
	position_list is dict returned from genebody_percentile
	position is 1-based genome coordinate
	'''
	samfile = pysam.Samfile(bam, "rb")
	aggreagated_cvg = collections.defaultdict(int)
	
	gene_finished = 0
	res_df = pd.DataFrame()
	
	for transcript_id in position_list:
		
		transcript_id_elements = transcript_id.split("_")
		transcript_id_overall_start = transcript_id_elements[1]
		transcript_id_overall_end = transcript_id_elements[2]
		transcript_id_overall = transcript_id_elements[3]
		#print(transcript_id_overall)
		
		list_values = position_list[transcript_id]
		chrom = list_values[0]
		strand = list_values[1]
		positions = list_values[2]
		coverage = {}
		for i in positions:
			coverage[i] = 0.0
		chrom_start = positions[0]-1
		if chrom_start <0: chrom_start=0
		chrom_end = positions[-1]
		try:
			samfile.pileup(chrom, 1,2)
		except:
			continue
			
		for pileupcolumn in samfile.pileup(chrom, chrom_start, chrom_end, truncate=True):
			ref_pos = pileupcolumn.pos+1
			if ref_pos not in positions:
				continue
			if pileupcolumn.n == 0:
				coverage[ref_pos] = 0
				continue				
			cover_read = 0
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_del: continue
				if pileupread.alignment.is_qcfail:continue 
				if pileupread.alignment.is_secondary:continue 
				if pileupread.alignment.is_unmapped:continue
				if pileupread.alignment.is_duplicate:continue
				cover_read +=1
			coverage[ref_pos] = cover_read
		tmp = [coverage[k] for k in sorted(coverage)]
		if strand == '-':
			tmp = tmp[::-1]
		for i in range(0,len(tmp)):
			aggreagated_cvg[i] += tmp[i]
		
		#print(transcript_id_overall)
		#print(len(tmp))
		total_sum = sum(tmp)
		if total_sum == 0: continue
		#print(tmp)
		mdf = {'transcript_id': transcript_id_overall,'percentile' : np.arange(100)+1, 'coverage': tmp}
		mdf = pd.DataFrame(mdf)
		mdf=mdf.pivot(index='transcript_id', columns='percentile', values='coverage')
		mdf.insert(0, 'transcript_id', transcript_id_overall)

		res_df = pd.concat([res_df, mdf]).reset_index(drop=True)
		

		gene_finished += 1
		
		if gene_finished % 100 == 0:
			print("\t%d transcripts finished\r" % (gene_finished), end=' ', file=sys.stderr)
			
	return 	aggreagated_cvg, res_df
	
	
def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input",action="store",type="string",dest="input_files",help='Input file(s) in BAM format. "-i" takes these input: 1) a single BAM file. 2) "," separated BAM files. 3) directory containing one or more bam files. 4) plain text file containing the path of one or more bam file (Each row is a BAM file path). All BAM files should be sorted and indexed using samtools.')
	parser.add_option("-r","--refgene",action="store",type="string",dest="ref_gene_model",help="Reference gene model in bed format. [required]")
	parser.add_option("-l","--minimum_length",action="store",type="int",default=101, dest="min_mRNA_length",help="Minimum mRNA length (bp). mRNA smaller than \"min_mRNA_length\" will be skipped. default=%default")
	parser.add_option("-f","--format",action="store",type="string",dest="output_format", default='pdf', help="Output file format, 'pdf', 'png' or 'jpeg'. default=%default")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",default='df_sum',help="Prefix of output files(s). [required]")
	(options,args)=parser.parse_args()

    
	if not (options.output_prefix and options.input_files and options.ref_gene_model):
		parser.print_help()
		sys.exit(0)

	if not os.path.exists(options.ref_gene_model):
		print('\n\n' + options.ref_gene_model + " does NOT exists" + '\n', file=sys.stderr)
		#parser.print_help()
		sys.exit(0)
	if options.min_mRNA_length < 100:
		print('The number specified to "-l" cannot be smaller than 100.' + '\n', file=sys.stderr)
		sys.exit(0)
		
	OUT1 = open(options.output_prefix  + ".geneBodyCoverage.tsv"	,'w')
	print("Percentile\t" + '\t'.join([str(i) for i in range(1,101)]), file=OUT1)
		
	printlog("Read BED file (reference gene model) ...")
	gene_percentiles = genebody_percentile(refbed = options.ref_gene_model, mRNA_len_cut = options.min_mRNA_length)
	
	#print(gene_percentiles)
	
	printlog("Get BAM file(s) ...")
	bamfiles = getBamFiles.get_bam_files(options.input_files)
	for f in bamfiles:
		print("\t" + f, file=sys.stderr)
	
	file_container =[]
	
	all_df = pd.DataFrame()
        #Create a master table of all filee
    
	for bamfile in bamfiles:
		printlog("Processing " + basename(bamfile) + ' ...')
		cvg, df_res = genebody_coverage(bamfile, gene_percentiles)
		
		#print(cvg)
		if len(cvg) == 0:
			print("\nCannot get coverage signal from " + basename(bamfile) + ' ! Skip', file=sys.stderr)
			continue
		tmp = valid_name(basename(bamfile).replace('.bam',''))	# scrutinize R identifer
		if file_container.count(tmp) == 0:
			print(tmp + '\t' + '\t'.join([str(cvg[k]) for k in sorted(cvg)]), file=OUT1)
		else:
			print(tmp + '.' + str(file_container.count(tmp)) + '\t' + '\t'.join([str(cvg[k]) for k in sorted(cvg)]), file=OUT1)
		
		#file_container.append(tmp)
		#print(df_res)
		mname = re.sub("_.*", "",basename(bamfile))
		#outfile = "df_genebodypercentile_"+ mname + ".tsv"
		#print(df_res)
		df_res.insert(0, 'File', mname)
		all_df = pd.concat([all_df, df_res]).reset_index(drop=True)
		#df_res.to_csv(outfile,sep = '\t',header = True,index = False,index_label = False)

	OUT1.close()
	all_df.to_csv("df_all.genebodypercentile.tsv",sep = '\t',header = True,index = False,index_label = False)
	
	#Convert into plottable format
	#Read as pandas object
	df = pd.read_csv("df_sum.geneBodyCoverage.tsv",sep = "\t")
	df_float = df.iloc[:,1:].apply(lambda row : pd.Series([float(i) for i in row]),axis = 1)
	
	epsilon = 1e-7  # small constant added to avoid division by zero
	df_trans = df_float.iloc[:,1:].apply(lambda row : pd.Series([(i -min(row))/(max(row) - min(row) + epsilon) for i in row]),axis = 1)
	df_trans['SampleId'] = df['Percentile']
	melted = pd.melt(df_trans,id_vars = "SampleId")
	
	sns.lineplot(x="variable", y="value",hue = 'SampleId',legend = False,data=melted).set(xlabel='5\'-3\' body percentile', ylabel='Coverage')
	plt.savefig('transcript_body_coverage_mqc.png', dpi=300)
	x = df_float.iloc[:,1:].apply(pearson_moment_coefficient,axis = 1)
	df_skew = pd.DataFrame({'SampleId' : list(df['Percentile']),'Skewness' : x })
	df_skew.to_csv("df_skewness.genebodypercentile.tsv",sep = '\t',header = True,index = False,index_label = False)

	

	

	
if __name__ == '__main__':
	main()
	
