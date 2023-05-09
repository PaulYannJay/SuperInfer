# -*- coding: utf-8 -*-
##### Code by Paul Jay #####

import sys
import os 
import getopt
import pandas as pd
import random
import time
import scipy.stats as stats
#from itertools import combinations
#from itertools import product 
#from sklearn.utils import check_random_state
#from sklearn.preprocessing import StandardScaler
#from sklearn.decomposition import PCA
#from sklearn.cluster import KMeans
#from sklearn.metrics import silhouette_score
#from sklearn.metrics import davies_bouldin_score 
#from sklearn.metrics.pairwise import distance_metrics
#from sklearn.metrics.pairwise import pairwise_distances
#from sklearn.experimental import enable_iterative_imputer
#from sklearn.impute import IterativeImputer
#from sklearn.metrics import pairwise_distances
#from sklearn import mixture
#from sklearn.cluster import AgglomerativeClustering
#from sklearn.manifold import MDS
#from sklearn.manifold import TSNE
import numpy as np
usage='''

//////////////////////////////////////////
Script usage:
python3 AssociationTest.py -g Genotypefile.geno -p PhenotypeFile -o OutputFile [-S SampleList]  [-s startPos -e endPos]

Parameters:
-g [String] The name of the file containinf the genotype
-o [String] The prefix of the output files
-p [String] The name of the file containing the phenotype

Options:
-S [String] The name of the file containing the name of the individuals to use for analyses. By default, all individuals are used. The individual file must contain the name of the individual to use, each in one line.
-s Start position in the scaffold to analyse
-e End position in the scaffold to analyse


Output:

Note: The Genotypefile.geno file should be in the format "Scaf Position GenotypeCodeSample1 GenotypeCodeSample2 ... GenotypeCodeSampleN" with genotype code being 0 for 0/0, 1 for 0/1 and 2 for 1/1 (only biallelic). All position must be genotyped (no "NA"). From a vcf file, to obtain it: 
    1 first, impute the genotype with beagle (eg. java -Xmx20g -jar ~/Software/beagle.05May22.33a.jar  gt=data.vcf.gz out=data.imputed.vcf window=1 overlap=0.5
    2. Transform to .geno file: python3 vcf2geno.py -i data.imputed.vcf.gz -o data.imputed.geno

////////////////////////////////////////////////

'''


def main(argv):
    global GenoFile
    global OutputFile
    global PhenoFile
    global IndFile
    global StartPos 
    global EndPos 
    global optionSubset 
    global optionSubsetPos 
    optionSubset=False
    optionSubsetPos=False
    try:
        opts, args = getopt.getopt(argv,"hg:o:p:s:e:")
    except getopt.GetoptError:
        print('One of the option does not exist !\n', usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-g"):
            GenoFile = arg
        elif opt in ("-p"):
            PhenoFile = arg
        elif opt in ("-S"):
            IndFile = arg
            optionSubset=True
        elif opt in ("-o"):
            OutputFile = arg
        elif opt in ("-s"):
            StartPos = arg
            optionSubsetPos=True
        elif opt in ("-e"):
            EndPos = arg
            optionSubsetPos=True
    print('Input file is ', GenoFile)
    print('Output file is ', OutputFile)

def cramers_V(var1,var2) :
      #dataset=np.vstack((var1, var2))
      #print(dataset)
      ContTab=stats.contingency.crosstab(var1, var2) #Get a contingency table from the genotype array and the phenotype array
      X2 = stats.chi2_contingency(ContTab.count, correction=False)[0] #Calculate chi2
      N = np.sum(ContTab.count) #Number of sample
      minimum_dimension = min(ContTab.count.shape)-1  #minimum number of dim
      result = np.sqrt((X2/N) / minimum_dimension) #Cramer's V
      return result
  
def Perform_CramerV_Association(File, Pheno):
    with open(File) as infile: #Read line by line (i.e. do not load the file in memory)
        next(infile)#Skip header 
        for line in infile: #For each genotype position
            lineArr=line.split()
            Scaff=lineArr[0] #Svaffold
            Pos=lineArr[1] #Position
            GenoArray=np.asarray(lineArr[2:])#Create an array from the genotype at this position (split on white space)
            result=cramers_V(GenoArray, Pheno) #Perform the analyses
            Outfile.write(Scaff+" "+ Pos+" "+ str(result)) #Write the output
            Outfile.write("\n")
            



### Initiation : Create Output Files###
if __name__ == "__main__":
    main(sys.argv[1:])
            
#textfileClustScore = open(OutputFile+".ClustScore", "w")
#textfileCluster = open(OutputFile+".cluster", "w")
#textfileHetero = open(OutputFile+".Hetero", "w")
#textfileHetero.write("Scaffold Start End No.variants No.cluster Cluster Hetero Homo")
#if (optionPI):
#    textfileHetero.write(" Pi\n")
#else:
#    textfileHetero.write("\n")
##textfileBic = open(OutputFile+".Bic", "w") #Used for gaussian mixture
##textfileBic.write("Scaffold Start End No.variants k1 k2 k3 k4 k5 k6\n")
#
#textfileDistance = open(OutputFile+".ClusterDistance", "w")
#textfileDistance.write("Scaffold Start End No.variants No.cluster Cluster1 Cluster2 Distance\n")
#
#if (optionDXY):
#    textfileDxy = open(OutputFile+".Dxy", "w")
#    textfileDxy.write("Scaffold Start End No.variants No.cluster Cluster1 Cluster2 Dxy\n")

def GetPheno(PhenoFile, GenoFile):
    with open(GenoFile) as infile: #Read line by line (i.e. do not load the file in memory)
        first=next(infile)#Skip header 
        header=first.split()
    NbSample=len(header)-2
    print(NbSample)
    Pheno=np.zeros(NbSample)
    Ind=0
    with open(PhenoFile) as infile: #Read line by line (i.e. do not load the file in memory)
        next(infile) #skip header
        for line in infile:
            IndLine=line.split() #tricks to get only the individual name, without the '\n'
            index =  [i for i in range(len(header)) if header[i] ==    IndLine[0]]    #Grep the sample index
            Pheno[index[0]-2]=IndLine[1]
    return Pheno

    #with open(PhenoFile)

def subset_Individual(GenoFile, IndFile):
    textfileGenoSub = open(OutputFile+".GenoSub", "w") #create a temporary geno file to store only the genotype of the focal samples
    with open(GenoFile) as infile: #We open the file line by line, as we are only interested in the first line
        firstline=next(infile) #Header of the geno file
        firstlineArr=firstline.split()
    IndIndex=[] #Indices of the individual to keep
    IndIndex.append(0) #Add the Scaffold name and variant position indices to the column to keep
    IndIndex.append(1)
    with open(IndFile) as infile: #Open the file containing the sample names of the sample to keep for analyses
        for line in infile:
            Ind=line.split()[0] #tricks to get only the individual name, without the '\n'
            index =  [i for i in range(len(firstlineArr)) if firstlineArr[i] ==    Ind]    #Grep the sample index
            IndIndex.append(int(index[0])) #Append the sample indice list
    IndIndex.sort() #Sort the indices (just in case...)
    with open(GenoFile) as infile: #reopen the geno file line by line
        for line in infile:
            linesplit=line.split() #Split the line (each value on the list correspond to the data of a sample
            Array=[linesplit[x]    for x in IndIndex] #Keep only the value that correspond to the indices of the sample to keep
            textfileGenoSub.write(" ".join(Array)) #write in the output file.
            textfileGenoSub.write("\n")
    print("subsetting samples: done !")
    textfileGenoSub.close()

### 

### Start Analyses ###
if (optionSubset):
    start_time = time.time()
    subset_Individual(GenoFile, IndFile)
    print("--- %s seconds:subsetting ---" % (time.time() - start_time))
    GenoFile=OutputFile+".GenoSub"


Outfile = open(OutputFile+".CramersV", "w")
Pheno=GetPheno(PhenoFile, GenoFile)
Perform_CramerV_Association(GenoFile, Pheno)
#write_header(MaxCluster)
#write_headerCluster()
#if (Method in "pca" or optionPCA):
#    OutputFilePCA=OutputFile+".pcaResult"
#    textfilepca = open(OutputFilePCA, "w")
#    write_headerPCA()

#OutputFileMDS=OutputFile+".MdsResult"
#textfilemds = open(OutputFileMDS, "w")
#write_headerMDS()
#
#OutputFileTSNE=OutputFile+".TsneResult"
#textfiletsne = open(OutputFileTSNE, "w")
#write_headerTSNE()

#if (WindType in "variant"):
#    Sliding_window_variant_overlap(GenoFile,WindSize, Slide)
#elif (WindType in "bp"):
#    Sliding_window_bp_overlap(GenoFile,WindSize, Slide)
#else:
#    print("Windows type must be 'variant' or 'bp' ")

if (optionSubset):
    os.remove(OutputFile+".GenoSub")

#textfileClustScore.close()
#textfileCluster.close()
#textfileHetero.close()
#textfileDistance.close()
#if (Method in "pca" or optionPCA):
#    textfilepca.close()
#if (optionDXY):
#    textfileDxy.close()
