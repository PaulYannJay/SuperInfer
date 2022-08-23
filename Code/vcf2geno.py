# -*- coding: utf-8 -*-

import sys
import getopt
import pandas as pd
#import matplotlib.pyplot as plt
import gzip
#import seaborn as sns
#sns.set()
from itertools import combinations
from sklearn.utils import check_random_state
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.metrics.pairwise import distance_metrics
from sklearn.metrics.pairwise import pairwise_distances
#from sklearn.externals.joblib import Parallel, delayed
from joblib import Parallel, delayed
import numpy as np
#import time
#import re

def main(argv):
	global vcfFile
	global GenoFile
	try:
		#opts, args = getopt.getopt(argv,"hg:o:w:k:t:k:",["ifile=","ofile="])
		opts, args = getopt.getopt(argv,"hi:o:")
	except getopt.GetoptError:
		print('Bug1 ! Usage: GenotypeCluster.py -i <file.vcf> -o <OutputFile.geno>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('Usage: GenotypeCluster.py -i <file.vcf> -o <OutputFile.geno>')
			sys.exit()
		elif opt in ("-i"):
			vcfFile = arg
		elif opt in ("-o"):
			GenoFile = arg
	print('Input file is ', vcfFile)
	print('Output file is ', GenoFile)


def formatLine(line):
	new_array=line.split()#Create an array from the genotype at this position (split on white space)
	#print(new_array[0])
	if (new_array[0].startswith('##')):
		next
	elif (new_array[0].startswith('#')):
		Line=new_array[0:2] + new_array[9:len(new_array)+1]
		for element in Line:
			textfile.write(str(element) + " ")
		textfile.write("\n")
	else:
		if "," not in new_array[4]: #Verify that the position is biallelic
			nbIndiv=len(new_array)
			for ind in list(range(10,nbIndiv+1)): #First individual at position 9 in vcf array
				if (new_array[ind-1].startswith('0/0') | new_array[ind-1].startswith('0|0')):
					new_array[ind-1]=0
				elif (new_array[ind-1].startswith('0/1') | new_array[ind-1].startswith('0|1')| new_array[ind-1].startswith('1|0')):
					new_array[ind-1]=1
				elif (new_array[ind-1].startswith('1/1') | new_array[ind-1].startswith('1|1')):
					new_array[ind-1]=2
				elif (new_array[ind-1].startswith('./.') | new_array[ind-1].startswith('.|.')):
					new_array[ind-1]="NaN"
				else:
					print(new_array[ind-1]," bug")
			Line=new_array[0:2] + new_array[9:len(new_array)+1]
			for element in Line:
				textfile.write(str(element) + " ")
			textfile.write("\n")

def vcf2geno(File): #No window overlap #Need to be improved
	'''
	Function to read a vcf file and change the genotype format
	'''
	#with open(File) as infile: #Read line by line (i.e. do not load the file in memory)
	if vcfFile.endswith("gz"):
		with gzip.open(File, "rt") as infile: #Read line by line (i.e. do not load the file in memory)
			for line in infile: #For each genotype position
				formatLine(line)
	else:
		with open(File) as infile: #Read line by line (i.e. do not load the file in memory)
			for line in infile: #For each genotype position
				formatLine(line)

if __name__ == "__main__":
	main(sys.argv[1:])

textfile = open(GenoFile, "w")
vcf2geno(vcfFile)
textfile.close()


