# -*- coding: utf-8 -*-
##### Code by Paul Jay #####

import sys
import os 
import getopt
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
python3 ComputeHetero.py -g Genotypefile.geno -o OutputFile -p PositionFile -s SampleFile ]

Parameters:
-g [String] The name of the file containinf the genotype
-p [String] The name of the file containing the region to analyse 
-o [String] The prefix of the output files

Options:
-s [String] The name of the file containing the name of the individuals to use for analyses. By default, all individuals are used. The individual file must contain the name of the individual to use, each in one line.


Output:

Note: The Genotypefile.geno file should be in the format "Scaf Position GenotypeCodeSample1 GenotypeCodeSample2 ... GenotypeCodeSampleN" with genotype code being 0 for 0/0, 1 for 0/1 and 2 for 1/1 (only biallelic). All position must be genotyped (no "NA"). From a vcf file, to obtain it: 
    1 first, impute the genotype with beagle (eg. java -Xmx20g -jar ~/Software/beagle.05May22.33a.jar  gt=data.vcf.gz out=data.imputed.vcf window=1 overlap=0.5
    2. Transform to .geno file: python3 vcf2geno.py -i data.imputed.vcf.gz -o data.imputed.geno

////////////////////////////////////////////////

'''


def main(argv):
    global GenoFile
    global OutputFile
    global IndFile
    global PosFile
    global optionSubset 
    optionSubset=False
    try:
        opts, args = getopt.getopt(argv,"hg:o:p:s:")
    except getopt.GetoptError:
        print('One of the option does not exist !\n', usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-g"):
            GenoFile = arg
        elif opt in ("-s"):
            IndFile = arg
            optionSubset=True
        elif opt in ("-o"):
            OutputFile = arg
        elif opt in ("-p"):
            PosFile = arg
    print('Input file is ', GenoFile)
    print('Output file is ', OutputFile)
           
def computeHetero(GenoFile, Column): #Function to compute the number of heterozygous position for each sample
    Geno=np.loadtxt(GenoFile, skiprows=0) #Open the genotype file (only the predefined colum (i.e. Samples)
    nPos=np.size(Geno, axis=0) #Get the number of position genotyped
    Count1=np.count_nonzero(Geno==1, axis=0) #COunt the number of occurence of '1'
    return nPos, Count1

def subset_Individual(GenoFile, IndFile):
    with open(GenoFile) as infile: #We open the file line by line, as we are only interested in the first line
        firstline=next(infile) #Header of the geno file
        firstlineArr=firstline.split()
        secondline=next(infile) #First line of data (used to get the name of the chromosome analsed
        secondlineArr=secondline.split()
    IndIndex=[] #Indices of the individual to keep
    with open(IndFile) as infile: #Open the file containing the sample names of the sample to keep for analyses
        for line in infile:
            Ind=line.split()[0] #tricks to get only the individual name, without the '\n'
            index =  [i for i in range(len(firstlineArr)) if firstlineArr[i] ==    Ind]    #Grep the sample index
            IndIndex.append(int(index[0])) #Append the sample indice list
    IndIndex.sort() #Sort the indices (just in case...)
    chrom=secondlineArr[0]
    return IndIndex, chrom  #IndIndex contain the index of the position of the samples to be analysed


def getNumberSample_And_Chrom(GenoFile): #Function to extract the number of samples and the focal chromosome
    with open(GenoFile) as infile: #We open the file line by line, as we are only interested in the first line
        firstline=next(infile) #Header of the geno file
        firstlineArr=firstline.split()
        secondline=next(infile)
        secondlineArr=secondline.split()
    IndIndex=range(2, len(firstlineArr)) #Create a sequence containing the index of the sample to analyse
    chrom=secondlineArr[0]
    return IndIndex, chrom 
############ Not used anymore #######################""""
def subset_Pos(GenoFile, StartPos, EndPos): #Function to subset the geno file to only analyse interesting position
    textfileGenoSub = open(OutputFile+"."+StartPos+"-"+EndPos+".geno", "w") #create a temporary geno file to store only the genotypes at the focal position 
    with open(GenoFile) as infile: #open the geno file line by line
        for line in infile:
            linesplit=line.split() #Split the line (each value on the list correspond to the data of a sample
            if (linesplit[1] > StartPos and linesplit[1] < EndPos): #Only keep position that fall within the fical position
                textfileGenoSub.write(line) #write in the temp file.
                textfileGenoSub.write("\n")
    textfileGenoSub.close()
#######################################################

def writeOutputHeader(GenoFile, Column): 
    with open(GenoFile) as infile: #We open the file line by line, as we are only interested in the first line
        firstline=next(infile) #Header of the geno file
        firstlineArr=firstline.split()
        Array=[firstlineArr[x]    for x in Column] #Keep only the value that correspond to the indices of the sample to keep
    Output.write(" ".join([str(firstlineArr[0]),"Start","End","nPos"])+ " ") 
    Output.write(" ".join(Array)+"\n")

def writeOutput(nPos, Count1, StartPos, EndPos): #write the output
    Output.write(" ".join([chrom,str(StartPos),str(EndPos),str(nPos)])) 
    for i in Count1: #Count1 contain the number of hetero position for each sample
        Output.write(" "+str(i))
    Output.write("\n")


def CreateSubSetGenoFiles(PosFile, GenoFile, Column): #FUnction to create subset of geno file
    Position=np.loadtxt(PosFile)
    SortPos=np.sort(Position, axis=0)
    Pos=0
    Geno=np.loadtxt(GenoFile, skiprows=1, usecols=Column, dtype=int)
    Pos=np.loadtxt(GenoFile, skiprows=1, usecols=1)
    for i in Position:
        PosInd=np.where((Pos>i[0]) & (Pos<i[1]))
        GenoSub=Geno[PosInd]
        np.savetxt(OutputFile+"."+str(int(i[0]))+"-"+str(int(i[1]))+".geno", GenoSub, fmt='%d') #        GenoSub=Geno[Geno[0, with open(GenoFile, 'r') as infile: #open the geno file line by line

#'
### 
### Initiation : Create Output Files###
if __name__ == "__main__":
    main(sys.argv[1:])
GenoFileOriginal=GenoFile #Copy the name of the original genofile
Output = open(OutputFile+".SampleHetero", "w")

### Start Analyses ###
if (optionSubset):
    start_time = time.time()
    Column, chrom=subset_Individual(GenoFile, IndFile) #Get the index of the sample to use
    print("--- %s seconds:subsetting samples---" % (time.time() - start_time))
else:
    Column, chrom=getNumberSample_And_Chrom(GenoFile)

writeOutputHeader(GenoFileOriginal, Column)
CreateSubSetGenoFiles(PosFile, GenoFile, Column) #For each position, create first the subsetted array
with open(PosFile) as infile: #We open the position file line by line
   for line in infile: #For each tupple of position
        linesplit=line.split() #Split the line (value on the list correspond to the start and end position
        StartPos=(linesplit[0])
        EndPos=(linesplit[1])
        print(line)
        GenoFile=OutputFile+"."+StartPos+"-"+EndPos+".geno" #The name of the subsetted array containing the genotype at this position
        nPos, Count1=computeHetero(GenoFile, Column)
        writeOutput(nPos, Count1, StartPos, EndPos)
        os.remove(OutputFile+"."+StartPos+"-"+EndPos+".geno")

