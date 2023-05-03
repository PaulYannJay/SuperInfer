# -*- coding: utf-8 -*-
##### Code by Paul Jay #####

import sys
import os 
import os
os.environ["MKL_NUM_THREADS"] = "10" 
os.environ["NUMEXPR_NUM_THREADS"] = "10" 
os.environ["OMP_NUM_THREADS"] = "10" 
import getopt
import random
import time
from itertools import combinations
from itertools import product 
from sklearn.utils import check_random_state
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.metrics import davies_bouldin_score 
from sklearn.metrics.pairwise import distance_metrics
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.metrics import pairwise_distances
from sklearn import mixture
from sklearn.cluster import AgglomerativeClustering
from sklearn.manifold import MDS
from sklearn.manifold import TSNE
import numpy as np
usage='''

//////////////////////////////////////////
Script usage:
python3 GenotypeCluster.py -g Genotypefile.geno -w WindowSize -s WindowSlide -t [variant/bp] -k MaxNumberOfCluster -o OutputFile -p Number_of_pca [-a] 

Parameters:
-g [String] The name of the file containinf the genotype
-w [Integer] The Window size. It defines the size of the windows for analyses
-s [Integer] The extent of the slide of windows
-o [String] The prefix of the output files
-t [String] Defines whether the window size is in variant number or in base pair. For variant number, use "-t variant". For base pair, use "-t bp".
-k [Integer] Defines the maximum number of cluster that will be determined and for which the clustering score will be calculated. For instance, if "-k 6" is used, the script will output, for each window, the clustering score of k=2,k=3,k=4,k=5 and k=6  
-p [Integer] Do the clustering on p principal components. So first the script perform a PCA, and then use p axes to do the clustering.

Options:
-S [String] The name of the file containing the name of the individuals to use for analyses. By default, all individuals are used. The individual file must contain the name of the individual to use, each in one line.

Output:
OutputFile.ClustScore --> The clutering score (Silhouette score or davies bouldin's score) for each k, for each window 
OutputFile.pca --> The position of each sample on each pca axes, for each window
OutputFile.cluster --> The affiliation of each sample (the cluster they belong to), for each number of cluster, for each window
OutputFile.clusterDistance --> The euclidian distance between each cluster, for each number of cluster, for each window
OutputFile.Hetero --> The proportion of heterozygous and homozygous site (only variant) in each cluster, for each number of cluster, for each window. If "-a Pi" is provided, this file also contains the result of Pi calculation

Note: The Genotypefile.geno file should be in the format "Scaf Position GenotypeCodeSample1 GenotypeCodeSample2 ... GenotypeCodeSampleN" with genotype code being 0 for 0/0, 1 for 0/1 and 2 for 1/1 (only biallelic). All position must be genotyped (no "NA"). From a vcf file, to obtain it: 
    1 first, impute the genotype with beagle (eg. java -Xmx20g -jar ~/Software/beagle.05May22.33a.jar  gt=data.vcf.gz out=data.imputed.vcf window=1 overlap=0.5
    2. Transform to .geno file: python3 vcf2geno.py -i data.imputed.vcf.gz -o data.imputed.geno

////////////////////////////////////////////////

'''


def main(argv):
    global GenoFile
    global OutputFile
    global WindType
    global MaxCluster
    global WindSize
    global Slide 
    global NoAxe 
    global optionPCA
    global MaxCompar 
    global optionSingleAxeClust
    global IndFile
    global optionSubset 
    optionSingleAxeClust=False
    optionSubset=False
    options="none"
    NoAxe=10 #Default number of pca axes used for clustering (if "-p" is specified)
    MaxCompar=100 #Default number of comparison used for the computation of Dxy and Pi
    try:
        opts, args = getopt.getopt(argv,"hg:o:w:s:k:t:k:i:p:aS:")
    except getopt.GetoptError:
        print('One of the option does not exist !\n', usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-g"):
            GenoFile = arg
        elif opt in ("-S"):
            IndFile = arg
            optionSubset=True
        elif opt in ("-o"):
            OutputFile = arg
        elif opt in ("-w"):
            WindSize = int(arg)
        elif opt in ("-s"):
            Slide = int(arg)
        elif opt in ("-k"):
            MaxCluster = int(arg)
        elif opt == "-p":
                NoAxe=int(arg)  #Number of pca axes to compute
        elif opt in "-a":
            optionSingleAxeClust=True
        elif opt in ("-t"):
            if (arg in "variant" or arg in "bp"):
                WindType = arg
            else:
                print("Error: Window type must be 'variant' or 'bp'\n", usage)
                sys.exit()
    print('Input file is ', GenoFile)
    print('Output file is ', OutputFile)


def Apply_PCA(Input, n_PC):
    '''
    Function to compute PCA on genomic data (input: array of genotype, number of pca axe to output, output: array of position of sample on x pca axes)
    Note: sample must be in row, variable in column (contrary to .vcf file or similar)
    '''
    pca = PCA(n_components=n_PC) #Define the parameter of the PCA function
    PC=pca.fit_transform(Input)
    return PC
    #PC_Weight = pca.explained_variance_ratio_ #Could be useful in the future...

def Apply_Kmeans(Input,n_cluster):
    '''
    Function to compute k-means (input: array of genotype of n samples or array of position of sample on n pca axes, number of cluster to compute; output: array of cluster attribution for each sample, depending on n_cluster [n_cluster rows, n columns], Center of clusters)
    '''
    kmeans = KMeans(n_clusters=n_cluster, n_init='auto').fit(Input) #Effect of max_iter not properly tested
    wt_kmeansclus = kmeans.predict(Input) #No weight
    return wt_kmeansclus, kmeans.cluster_centers_

#def computeGaussianMixture(Data, MaxCluster):
#    n_components_range = range(1, 7)
#    bic=[]
#    for n_components in n_components_range:
#        gmm = mixture.GaussianMixture(n_components=n_components, covariance_type="full")
#        gmm.fit(Data)
#        bic.append(gmm.bic(Data))
#    #print(bic)
#    return bic

def Define_Cluster_WithPCA(Data, MaxCluster):
    '''
    Calculate the silhouette score for different k on pca Output (first compute pca, then the clustering score) with k between 2 and MaxCluster (input: array of genotype or position on pca axes, Number of cluster ;output: Array of clustering score score for each k, list of cluster attribution, Position of sample on PC, Cluster centers)
    '''
    Clusters = list(range(2,MaxCluster+1)) #Number of cluster to calculate. (e.g. if "-k 5" is provided, the analyses are done for 2,3,4 and 5 clusters
    ScoreList=[]
    ClusterList=[]
    ClusterCenterList=[]
    PCAaxes=min(np.shape(Data)[1], NoAxe) #Number of PCA axis to return
    PC=Apply_PCA(Data,PCAaxes) #Compute PCA
    #Bic=computeGaussianMixture(PC,MaxCluster) #Could be useful in the future
    for cluster in Clusters: #For each number of cluster
        KM, ClusterCenter=Apply_Kmeans(PC,cluster) #Kmeans and cLuster center
        pc="1:"+str(NoAxe)
        KMOut=np.concatenate((cluster,pc,KM), axis=None)
        PCstr=np.repeat([pc], NoAxe)
        #ClusterCenterAxe=np.append(PCstr, ClusterCenter, axis=0)
        ClusterCenterAxe=np.vstack([PCstr, ClusterCenter])
        #KMOut=np.append(cluster,KM) #Append the result to the list # Why using a numpy array ? To be tested
        ClusterList.append(KMOut)#KMOut=np.concatenate((Pref,NoPc[:, None],PC), axis=1)
        ClusterCenterList.append(ClusterCenterAxe) #List of Cluster center
        #ClusterCenterList.append(ClusterCenter) #List of Cluster center
            #return ScoreList,ClusterList,PC,ClusterCenterList,Bic #Same, could be useful
    return ClusterList,PC,ClusterCenterList

def Define_Cluster_WithPCA_SingleAxe(Data, MaxCluster):
    '''
    Calculate the silhouette score for different k on pca Output (first compute pca, then the clustering score) with k between 2 and MaxCluster (input: array of genotype or position on pca axes, Number of cluster ;output: Array of clustering score score for each k, list of cluster attribution, Position of sample on PC, Cluster centers)
    '''
    Clusters = list(range(2,MaxCluster+1)) #Number of cluster to calculate. (e.g. if "-k 5" is provided, the analyses are done for 2,3,4 and 5 clusters
    ScoreList=[]
    ClusterList=[]
    ClusterCenterList=[]
    PCAaxes=min(np.shape(Data)[1], NoAxe) #Number of PCA axis to return
    PC=Apply_PCA(Data,PCAaxes) #Compute PCA
    for pc in range(0,NoAxe):
        for cluster in Clusters: #For each number of cluster
            KM, ClusterCenter=Apply_Kmeans(PC[:,[pc]],cluster) #Kmeans and cLuster center
            KMOut=np.concatenate((cluster,pc+1,KM), axis=None)
            #ClusterCenterAxe=np.concatenate((pc, ClusterCenter), axis=None)
            ClusterCenterAxe=np.append([[pc+1]], ClusterCenter, axis=0)
            ClusterList.append(KMOut)#KMOut=np.concatenate((Pref,NoPc[:, None],PC), axis=1)
            #ClusterCenterList.append(ClusterCenter) #List of Cluster center
            ClusterCenterList.append(ClusterCenterAxe) #List of Cluster center
    return ClusterList,PC,ClusterCenterList

def Compute_Heterozygosity(Geno, WindowPos, Cluster):
    ''' Function to estimate the number of heterozygous and homozygous site in each cluster. (Input: Array of genotype, WindowPosition, List of cluster attribution; output: Directly write the result'''
    Cluster=np.asarray(Cluster)
    #Cluster=Cluster[:,1:] #Grep only the genotype at each loci in the window (i.e. remove the cluster number)
    for vec in Cluster: #For each number of cluster #line contains the cluster attribution of each sample
        line=vec[2:] #Only genotype
        PC=vec[1] #Axis on which analyses where computed
        nbClust=vec[0] #Number of cluster
        myset= set(line) # vector of unique cluster
        for element in myset: #For each cluster
            index_pos_list = [ i for i in range(len(line)) if line[i] ==    element] #Grep de index of the individual belonging to this cluster
            HeteroCountArray=[] #create an empty array that will be used to store the number of heterozygous position in each individual of each cluster
            HomoCountArray=[] #create an empty array that will be used to store the number of homozygous position in each individual of each cluster
            for ind in index_pos_list: #For each individual
                Count1=np.count_nonzero(Geno[ind-1]==1) #Number of Heterozygous position
                Count2=np.count_nonzero(Geno[ind-1]==2) #Number of Homozygous derived allele
                Count0=np.count_nonzero(Geno[ind]==0) #Number of Homozygous ancestral allele
                CountHomo=Count0+Count2 #Number of homozygous site (among the variant site
                HeteroCountArray.append(Count1)#Append the main array with the genotype at this position
                HomoCountArray.append(CountHomo)#Append the main array with the genotype at this position
            NoOfPosition=WindowPos[2]-WindowPos[1] #Number of site considered in this window
            UnvariantSite=NoOfPosition-np.shape(Geno)[0] #Number of invariant site
            meanHetero=np.mean(HeteroCountArray)/NoOfPosition  #Mean number of heterozygous position in this cluster
            meanHomo=(np.mean(HomoCountArray) + UnvariantSite)/NoOfPosition #Mean number of homozygous positions
            Line2write=np.append(WindowPos,[nbClust, PC, element,meanHetero, meanHomo])
            for col in Line2write: #Write the result
                textfileHetero.write(str(col) + " ")
            textfileHetero.write("\n")


def Compute_analyses(Array, CurrScaffold, Start, End):
    ''' This is the main fonction the call the different function to perform different analyses'''
    WindowPos=[CurrScaffold,Start,End, np.shape(Array)[0]] #Position of the window 
    if np.shape(Array)[0]>4: #Verify that their is at least n loci in the window (if not, can't compute PCA with 10 axes)
        #Pos=Array[:,1].astype(float) #Grep the position of the loci in the window 
        Array=Array[:,2:] #Grep only the genotype at each loci in the window (i.e. remove the loci positon and scaffold)
        Array=Array.astype(float).transpose() #Define genotye as float and transpose the matrix
        if(np.isnan(Array).any()): #No "NA" is allowed in the geno file
            print("Fatal error: The window: ", Start, "-", Start+WindSize, " contain NaN values")
            sys.exit()

        if(optionSingleAxeClust):
            Cluster, PC, ClusterCenterList=Define_Cluster_WithPCA_SingleAxe(Array,MaxCluster) #Perform the PCA, the clustering and estimate the score for 2<=k<=MaxCluster
        else:
            Cluster, PC, ClusterCenterList=Define_Cluster_WithPCA(Array,MaxCluster) #Perform the PCA, the clustering and estimate the score for 2<=k<=MaxCluster
        write_pca(PC, WindowPos);
        write_cluster(WindowPos,Cluster)
        write_clusterDistance(ClusterCenterList, WindowPos)
        Compute_Heterozygosity(Array, WindowPos, Cluster)

def Sliding_window_bp_overlap(File, WindSize, Slide): 
    '''
    Function to split the .geno file in window and calculate the silhoutette score or other state in each window (parameter: window size. Input: genotype array, Output: Clustering scores in each window) 
    The function read the geno file line by line, so it can deal with very large file.
    Slide on bp
    '''
    with open(File) as infile: #Read line by line (i.e. do not load the file in memory)
        next(infile)#Skip header 
        secondline=next(infile)# Get the first line of the data to store info
        Array=[] #create an empty array that will be used to store the genotype of individual in each window
        new_array=secondline.split()#Create an array from the genotype at this position (split on white space)
        Array.append(new_array) #push the first loci in the new array
        Start=(int(new_array[1]) // WindSize) * WindSize # Start for sliding window (not 1 or the position of the first variant).  
        CurrScaffold=new_array[0]
        for line in infile: #For each genotype position
            new_array=line.split()#Create an array from the genotype at this position (split on white space)
            if (new_array[0] == CurrScaffold and int(new_array[1]) <= (Start+WindSize)): #If the loci is on the same scaffold and within the current window
                Array.append(new_array)#Append the main array with the genotype at this position
            else: #The variant not fall in the current window. It could be on another scaffold or in a next window on the same scaffold. We write the result of the analyse for this window and move to the next window
                print("Current position: Scaffold=", CurrScaffold, " Position=", Start)
                ArrayNP=np.asarray(Array) #Create a numpy array from the main array 
                End=Start+WindSize
                Compute_analyses(ArrayNP, CurrScaffold, Start, End)
                if (new_array[0] == CurrScaffold): #If the loci was on the same scaffold but not on the current window
                    Array= [x for x in Array if float(x[1])>Start+Slide]
                    NumberEmpty=((int(new_array[1])-(Start+WindSize)) // Slide) #Number of empty window
                    Start=Start+Slide #Change the start position of the window
                    End=Start+WindSize
                    for Wind in list(range(0,NumberEmpty)): #Write empty window
                        ArrayNP=np.asarray(Array) #Create a numpy array from the main array 
                        Compute_analyses(ArrayNP, CurrScaffold, Start, End)
                        Start=Start+Slide #Change the start position of the window
                        End=Start+WindSize
                        Array= [x for x in Array if float(x[1])>Start]
                    Array.append(new_array)#Append the main array with the genotype at this position
                else: #The window fall in a new scaffold
                    Array=[] #Reinitialise the main array (new window)
                    Array.append(new_array) #push the first loci in the new array
                    CurrScaffold = new_array[0] #Define the new scaffold/chromosome. If the loci was on the same scaffold, this does noit change anything.
                    Start=(int(new_array[1]) // WindSize) * WindSize # Start for sliding window (not 1 or the position of the first variant).  
                                
def Sliding_window_variant_overlap(File, WindSize, slide): 
    '''
    Function to split the .geno file in window and calculate the silhoutette score or other state in each window (parameter: window size. Input: genotype array, Output: Clustering scores in each window) 
    The function read the geno file line by line, so it can deal with very large file.
    Slide on variant
    '''
    with open(File) as infile: #Read line by line (i.e. do not load the whole file in memory)
        firstline=next(infile)#Skip header 
        secondline=next(infile)# Get the first line of the data to store info
        Array=[] #create an empty array that will be used to store the genotype of individual in each window
        new_array=secondline.split()#Create an array from the genotype at this position (split on white space)
        Array.append(new_array) #push the first loci in the new array
        Start=int(new_array[1]) # Start for sliding window (not 1, the position of the first variant).  
        CurrScaffold=new_array[0]
        Inc=1 #Increment for sliding window (here based on SNP number). When Inc= Window size, we move to the next window
        for line in infile: #For each loci 
            new_array=line.split()#Create an array from the genotype at this position (split on white space)
            if (new_array[0] == CurrScaffold and Inc < WindSize) : #if the focal loci is on the same scaffold as the previous one and within the same window
                Array.append(new_array)#Append the main array with the genotype at this position
                End=int(new_array[1]) # Record the position of the loci (for output) 
                Inc += 1 #Increment
            else: # The focal loci in on a new scaffold or is on a new window. We write the result of the analyse for this window and move to the next window
                print("Current position: Scaffold=", CurrScaffold, " Position=", Start)
                ArrayNP=np.asarray(Array) #Create a numpy array from the main array 
                #Compute_analyses(ArrayNP, CurrScaffold, Start, WindSize)
                Compute_analyses(ArrayNP, CurrScaffold, Start, End)
                if (new_array[0] != CurrScaffold): #if the focal loci is on a new scaffold 
                    Inc=1 #Restart the loci counter
                    Array=[] #create an empty array that will be used to store the genotype of individual in each window
                    Array.append(new_array) #push the first loci of the scaffold in the new array
                    #Array=np.asarray([new_array]) #create an empty array that will be used to store the genotype of individual in each window
                    CurrScaffold = new_array[0] #Define the new scaffold/chromosome. If the loci was on the same scaffold, this does noit change anything.
                    Start=int(new_array[1]) # Record the position of the loci (for output) 
                else: #if not we just consider a new window in the same scaffold. In this case, we slide the window, so we remove certain loci and change the loci counter
                    Array=Array[Slide:][:] ### Subset the array of genotype to remove the first loci 
                    Array.append(new_array)#Append the main array with the genotype at this locus 
                    Start=int(Array[0][1]) # Start for sliding window (not 1 or the position of the first variant).  
                    Inc=Inc-Slide+1 #modify the loci counter
                End=int(new_array[1]) # Record the position of the loci (for output) #Normally not used

def subset_Individual(GenoFile, IndFile):
    textfileGenoSub = open(OutputFile+".sub.geno", "w") #create a temporary geno file to store only the genotype of the focal samples
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
            
#### Write functions ###
def write_clusterDistance(ClusterCenterList, WindowPos):
    for NCluster in range(len(ClusterCenterList)): #For all number of cluster that might be considered (e.g. 2, 3, 4, and 5 if k was set to 5)
        distance=pairwise_distances(ClusterCenterList[NCluster][1:]) #Calculate the pairwise euclidean distance between the cluster centroid
        ClusterComb=list(combinations(range(len(distance)),2)) #All combination of cluster
        for Comb in ClusterComb: #For each combination of cluster
            for col in WindowPos: #Print the window position
                textfileDistance.write(str(col) + " " )
            textfileDistance.write(str(ClusterCenterList[NCluster][0][0])+ " " + str(len(distance)) + " " + str(Comb[0]) + " " + str(Comb[1]) + " " + str(distance[Comb[0],Comb[1]]) + "\n") #Print the distance between the centroid

def write_cluster(WindowPos, Cluster):
    Cluster=np.asarray(Cluster)
    Pref=np.repeat([WindowPos], Cluster.shape[0], axis=0)
    Cluster=np.concatenate((Pref,Cluster), axis=1)
    for line in Cluster:
        for element in line:
            textfileCluster.write(str(element) + " ")
        textfileCluster.write("\n")

def write_pca(PC, WindowPos):
    PC=PC.transpose()
    NoPc = np.array(np.arange(1,PC.shape[0]+1))
    Pref=np.repeat([WindowPos], PC.shape[0], axis=0)
    PC=np.concatenate((Pref,NoPc[:, None],PC), axis=1)
    for line in PC:
        for element in line:
            textfilepca.write(str(element) + " ")
        textfilepca.write("\n")


def write_headerPCA():
    with open(GenoFile) as infile: #Read line by line (i.e. do not load the file in memory)
        header=next(infile)#Skip header 
        new_array=header.split()#Create an array from the genotype at this position (split on white space)
    textfilepca.write("Scaffold Start End No.variants PC")
    for ind in list(range(2,len(new_array))):
        textfilepca.write(" " + str(new_array[ind]))
    textfilepca.write("\n")

def write_headerCluster():
    with open(GenoFile) as infile: #Read line by line (i.e. do not load the file in memory)
        header=next(infile)#Skip header 
        new_array=header.split()#Create an array from the genotype at this position (split on white space)
    textfileCluster.write("Scaffold Start End No.variants No.cluster Axes")
    for ind in list(range(2,len(new_array))):
        textfileCluster.write(" " + str(new_array[ind]))
    textfileCluster.write("\n")

### Initiation : Create Output Files###
if __name__ == "__main__":
    main(sys.argv[1:])
            
#textfileClustScore = open(OutputFile+".ClustScore", "w")
textfileCluster = open(OutputFile+".cluster", "w")
textfileHetero = open(OutputFile+".Hetero", "w")
textfileHetero.write("Scaffold Start End No.variants No.cluster Axes Cluster Hetero Homo\n")
textfileDistance = open(OutputFile+".ClusterDistance", "w")
textfileDistance.write("Scaffold Start End No.variants Axes No.cluster Cluster1 Cluster2 Distance\n")

### Start Analyses ###
if (optionSubset):
    start_time = time.time()
    subset_Individual(GenoFile, IndFile)
    print("--- %s seconds:subsetting ---" % (time.time() - start_time))
    GenoFile=OutputFile+".sub.geno"

write_headerCluster()
OutputFilePCA=OutputFile+".pcaResult"
textfilepca = open(OutputFilePCA, "w")
write_headerPCA()

if (WindType in "variant"):
    Sliding_window_variant_overlap(GenoFile,WindSize, Slide)
elif (WindType in "bp"):
    Sliding_window_bp_overlap(GenoFile,WindSize, Slide)
else:
    print("Windows type must be 'variant' or 'bp' ")

#if (optionSubset):
#    os.remove(OutputFile+".sub.geno")

#textfileClustScore.close()
textfileCluster.close()
textfileHetero.close()
textfileDistance.close()
textfilepca.close()
