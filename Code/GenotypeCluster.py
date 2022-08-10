# -*- coding: utf-8 -*-
##### Code by Paul Jay #####

import sys
import getopt
import pandas as pd
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
import numpy as np
usage='''

//////////////////////////////////////////
Script usage:
python3 GenotypeCluster.py -g Genotypefile.geno -w WindowSize -s WindowSlide -t [variant/bp] -k MaxNumberOfCluster -o OutputFile [-p Number_of_pca] [-a pca,dxy,pi] [-c Maximum_comparison] [-f Silhouette]

Help: The Genotypefile.geno file is in the format "Scaf Position GenotypeCodeSample1 GenotypeCodeSample2 ... GenotypeCodeSampleN" with genotype code being 0 for 0/0, 1 for 0/1 and 2 for 1/1 (only biallelic). All position must be genotyped (no "NA"). From a vcf file, to obtain it: 
	1 first, impute the genotype with beagle (eg. java -Xmx20g -jar ~/Software/beagle.05May22.33a.jar  gt=data.vcf.gz out=data.imputed.vcf window=1 overlap=0.5
	2. Transform to .geno file: python3 vcf2geno.py -i data.imputed.vcf.gz -o data.imputed.geno
The WindowSize must be an integer. It defines the size of the windows
The WindowSlide must be an integer. It define the extent of the slide of windows
The -t parameter defines wheter the window size is in variant number or in base pair. For variant number, use "-t variant". For base pair, use "-t bp".
The MaxNumberOfCluster must be an integer. It defines the maximum number of cluster that will be determined and for which the clustering score will be calculated. For instance, if "-k 6" is used, the script will output, for each window, the clustering score of k=2,k=3,k=4,k=5 and k=6  
By default, the clustering score used is the Davies Bouldin score. To compute the Silhouette score instead, use "-f Silhouette"
If "-p" is specified, the clustering is done on pca output instead on directly on genotype. So first the script compute pca, and then use Number_of_pca axes to do the clustering. Number_of_pca must be an integer.
If "-a pca" is specified, the clustering is done on genotype, but the pca are also computed and provided as an output
If "-a dxy" is specified, the script compute Dxy between each cluster. Be careful this can be pretty time comsuming. The euclidean distance calculated by the script between all cluster give a very similare result that Dxy and is much faster to compute
If "-a pi" is specified, the script compute Pi of each cluster. Be careful this can be pretty time comsuming 
"-c" determine the number of comparison done for the computation of Dxy and Pi. The default is 100. This allow to avoid to perform all haplotype comparison, which is very time consuming and useless most of the time.

'-o' specify the prefix of the output files

Output:
OutputFile.ClustScore --> The clutering score (Silhouette score or davies bouldin's score) for each k, for each window 
OutputFile.Hetero --> The proportion of heterozygous and homozygous site (only variant) in each cluster, for each number of cluster, for each window. If "-a Pi" is provided, this file also contains the result of Pi calculation
OutputFile.pca --> The position of each sample on each pca axes, for each window
OutputFile.cluster --> The affiliation of each sample (the cluster they belong to), for each number of cluster, for each window
OutputFile.clusterDistance --> The euclidian distance between each cluster, for each number of cluster, for each window
OutputFile.Dxy --> The nucleotide distance between each cluster, for each number of cluster, for each window

////////////////////////////////////////////////

'''


def main(argv):
	global GenoFile
	global OutputFile
	global WindType
	global MaxCluster
	global WindSize
	global Slide 
	global Method 
	global NoAxe 
	#global options 
	global optionPCA
	global ScoreFct 
	global optionDXY
	global MaxCompar 
	global optionPI
	optionPCA=False #Determine whether we perform PCA even if the clustering is based directly on genotype
	optionDXY=False #Determine whether we perform DXY 
	optionPI=False
	Method="direct"
	options="none"
	ScoreFct=davies_bouldin_score
	NoAxe=10 #Default number of pca axes used for clustering (if "-p" is specified)
	MaxCompar=100 #Default number of comparison used for the computation of Dxy and Pi
	try:
		opts, args = getopt.getopt(argv,"hg:o:w:s:k:t:k:m:i:p:a:c:f:")
	except getopt.GetoptError:
		print('One of the option do not exist !\n', usage)
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print(usage)
			sys.exit()
		elif opt in ("-g"):
			GenoFile = arg
		elif opt in ("-o"):
			OutputFile = arg
		elif opt in ("-w"):
			WindSize = int(arg)
		elif opt in ("-s"):
			Slide = int(arg)
		elif opt in ("-k"):
			MaxCluster = int(arg)
		elif opt in ("-c"):
			MaxCompar = int(arg)
		elif opt == "-p":
				Method = "pca" 
				NoAxe=int(arg)  #Number of pca axes to compute
		elif opt in ("-t"):
			if (arg in "variant" or arg in "bp"):
				WindType = arg
			else:
				print("Error: Window type must be 'variant' or 'bp'\n", usage)
				sys.exit()
		elif opt in ("-f"):
			if (arg in "Silhouette"):
				ScoreFct = silhouette_score
			else:
				print("Error: the alternative score function must be 'Silhouette'\n", usage)
				sys.exit()
		elif opt in ("-a"):
			print("Option(s) used:", arg)
			Options=arg.split(',')
			if ("pca" in Options or "dxy" in Options or "pi" in Options):
				if ("pca" in Options):
					optionPCA=True
				if ("dxy" in Options):
					optionDXY=True
				if ("pi" in Options):
					optionPI=True
			else:
				print(arg)
				print("-a must be followed with 'pca','pi' and/or 'dxy'. These are the only option available yet\n", usage)
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
	kmeans = KMeans(n_clusters=n_cluster, random_state=0, max_iter=1000).fit(Input) #Effect of max_iter not properly tested
	wt_kmeansclus = kmeans.predict(Input) #No weight
	return wt_kmeansclus, kmeans.cluster_centers_

def Define_Cluster_and_Score_WithPCA(Data, MaxCluster, ScoreFct):
	'''
	Calculate the silhouette score for different k on pca Output (first compute pca, then the clustering score) with k between 2 and MaxCluster (input: array of genotype or position on pca axes, Number of cluster ;output: Array of clustering score score for each k, list of cluster attribution, Position of sample on PC, Cluster centers)
	'''
	Clusters = list(range(2,MaxCluster+1)) #Number of cluster to calculate. (e.g. if "-k 5" is provided, the analyses are done for 2,3,4 and 5 clusters
	ScoreList=[]
	ClusterList=[]
	ClusterCenterList=[]
	PCAaxes=min(np.shape(Data)[1], NoAxe) #Number of PCA axis to return
	PC=Apply_PCA(Data,PCAaxes) #Compute PCA
	for cluster in Clusters: #For each number of cluster
		KM, ClusterCenter=Apply_Kmeans(PC,cluster) #Kmeans and cLuster center
		KMOut=np.append(cluster,KM) #Append the result to the list # Why using a numpy array ? To be tested
		ClusterList.append(KMOut)#KMOut=np.concatenate((Pref,NoPc[:, None],PC), axis=1)
		ClusterCenterList.append(ClusterCenter) #List of Cluster center
		Score= ScoreFct(PC, KM) #Calculate the score of the clustering (by default, DBS) 
		ScoreList.append(DBS)
	return ScoreList,ClusterList,PC,ClusterCenterList


def Define_Cluster_and_Score_WithoutPCA(Data, MaxCluster, ScoreFct):
	'''
	Calculate the silhouette score for different k, with k between 2 and MaxCluster (input: array of genotype or position on pca axes, Number of cluster ; output: Array of silhouette score, list of cluster attribution, Cluster centers)
	'''
	Clusters = list(range(2,MaxCluster+1))#Number of cluster to calculate. (e.g. if "-k 5" is provided, the analyses are done for 2,3,4 and 5 clusters
	ScoreList=[]
	ClusterList=[]
	ClusterCenterList=[]
	for cluster in Clusters: #For each number of cluster
		KM, ClusterCenter=Apply_Kmeans(Data,cluster) #Calculate Kmean
		KMOut=np.append(cluster,KM) #Store the result
		ClusterList.append(KMOut) #List of cluster centroid
		ClusterCenterList.append(ClusterCenter)
		Score = ScoreFct(Data, KM) #Calculate the score of the clustering (by default, DBS)
		ScoreList.append(Score)
	return ScoreList, ClusterList,ClusterCenterList

def Compute_Heterozygosity(Geno, WindowPos, Cluster):
	''' Function to estimate the number of heterozygous and homozygous site in each cluster. (Input: Array of genotype, WindowPosition, List of cluster attribution; output: Directly write the result'''
	Cluster=np.asarray(Cluster)
	Cluster=Cluster[:,1:] #Grep only the genotype at each loci in the window (i.e. remove the loci positon and scaffold)
	for line in Cluster: #For each number of cluster #line contains the cluster attribution of each sample
		myset= set(line) # vector of unique cluster
		for element in myset: #For each cluster
			index_pos_list = [ i for i in range(len(line)) if line[i] ==	element] #Grep de index of the individual belonging to this cluster
			HeteroCountArray=[] #create an empty array that will be used to store the number of heterozygous position in each individual of each cluster
			HomoCountArray=[] #create an empty array that will be used to store the number of homozygous position in each individual of each cluster
			for ind in index_pos_list: #For each individual
				Count1=np.count_nonzero(Geno[ind]==1) #Number of Heterozygous position
				Count2=np.count_nonzero(Geno[ind]==2) #Number of Homozygous derived allele
				Count0=np.count_nonzero(Geno[ind]==0) #Number of Homozygous ancestral allele
				CountHomo=Count0+Count2 #Number of homozygous site (among the variant site
				HeteroCountArray.append(Count1)#Append the main array with the genotype at this position
				HomoCountArray.append(CountHomo)#Append the main array with the genotype at this position
			NoOfPosition=WindowPos[2]-WindowPos[1] #Number of site considered in this window
			UnvariantSite=NoOfPosition-np.shape(Geno)[0] #Number of invariant site
			#meanHetero=np.mean(HeteroCountArray)/(WindowPos[2]-WindowPos[1]) 
			meanHetero=np.mean(HeteroCountArray)/NoOfPosition  #Mean number of heterozygous position in this cluster
			#meanHomo=np.mean(HomoCountArray)/(WindowPos[2]-WindowPos[1])
			meanHomo=(np.mean(HomoCountArray) + UnvariantSite)/NoOfPosition #Mean number of homozygous positions
			#Line2write=WindowPos.append([len(myset),element,meanHetero])
			Line2write=np.append(WindowPos,[len(myset),element,meanHetero, meanHomo])
			for col in Line2write: #Write the result
				textfileHetero.write(str(col) + " ")
			if (optionPI): #If the option for calculating pi is provided
				Pi=Compute_Pi(Geno[index_pos_list], MaxCompar)/(WindowPos[2]-WindowPos[1]) #Compute Pi
				#Pi=Compute_Pi(GenoSub, MaxCompar)/(WindowPos[2]-WindowPos[1])
				textfileHetero.write(str(Pi) + " ")
			textfileHetero.write("\n")

def Compute_Dxy(Geno, WindowPos, Cluster, MaxCompar):
	''' Function to calculate Dxy. Use several tricks to fasten the calculation and made several approximation because the calculation are for unphased diploid. Could not be 100% accurate, need to be ytested. (Input: Array of genotype, WindowPosition, List of cluster attribution, Maximul number of comparison to perform; output: Directly write the result'''
	Cluster=np.asarray(Cluster)
	Cluster=Cluster[:,1:] #Grep only the genotype at each loci in the window (i.e. remove the loci positon and scaffold)
	for line in Cluster: #For each number of cluster
		myset= set(line) # vector of unique cluster
		CombClust=list(combinations(myset,2)) # All comvbination of two cluster (e.g. 0-1 0-2 1-2)
		#CombClust=combinations(myset,2)
		for item in CombClust: #Pour chaque combination
			index_pos_list0 = [ i for i in range(len(line)) if line[i] ==    item[0]]	#Individual in the first cluster
			index_pos_list1 = [ i for i in range(len(line)) if line[i] ==    item[1]]	#Individual in the second cluster	
			AllComb=list(product(index_pos_list0, index_pos_list1)) #All possble diploid comparison between the samples of the two cluster
			NbCompar=min(MaxCompar, len(AllComb)) #Number of comparison to perform (either the maximum defined in the command line or All)
			RandomComb=random.sample(AllComb, NbCompar) #Random sample of comparison
			DxySum=0 #Dxy sum, to be incremented
			for Comb in RandomComb: #For all comparisons (Comb[1] and Comb[0] represent the index of the sample to compare)
				Diff=abs(Geno[Comb[0]]-Geno[Comb[1]])#Number of difference between the two diploid. However, this is not sufficient to calculate dxy. Indeed, the number of difference between two heterozygote (genotype == "1") is 0, but in fact, their is in average 0.5 difference between the haploid genomes of two diploid heterozygote (Aa vs Aa --> could be A vs A, a vs a, A vs a, or a vs A ). So we need a fast trick to get the true number of difference. here is a way:
				ind1_No2=np.where(Geno[Comb[1]]==2,0,Geno[Comb[1]]) #Replace the homozygote derived by 0. This give a 0,1 array where 0 mean homozygote and 1 heterozygote
				ind0_No2=np.where(Geno[Comb[0]]==2,0,Geno[Comb[0]]) # Same
				maxInd_No2=np.maximum(ind1_No2, ind0_No2) #Only get position where there is one heterozygote individual (0,1 array)
				NbDiff=sum(np.maximum(Diff,maxInd_No2)/2) #Determine the maximum between this array of heterozygote indiv and the difference between the two haploids. In position where their is heterozygote individuals, this is alway 1, elsewhere it could be 2 or 0, which is what we expect. We devide by two two get the mean pairwise difference between haploid genomes, not diploid genomes
				DxySum=DxySum+NbDiff
			Dxy=(DxySum/NbCompar)/(WindowPos[2]-WindowPos[1]) #Mean number of difference between haploid genome devided by number of site
			for col in WindowPos:
				textfileDxy.write(str(col) + " ")
			textfileDxy.write(str(len(myset)) + " ")
			#textfileDxy.write(str(item) + " ")
			textfileDxy.write(str(item[0]) + " " + str(item[1]) + " ")
			textfileDxy.write(str(Dxy) + " ")
			textfileDxy.write("\n")

def Compute_Pi(Geno, MaxCompar):
	''' Function to calculate Pi. Use several tricks to fasten the calculation and made several approximation because the calculation are for unphased diploid. Could not be 100% accurate, need to be ytested. (Input: Array of genotype, Maximul number of comparison to perform; output: Pi'''
	PiSum=0 #To be incremented
	myset=set(range(np.shape(Geno)[0]))
	CombInd=list(combinations(myset,2)) # All combination of two individuals (e.g. 0-1 0-2 1-2)
	NbCompar=min(MaxCompar, len(CombInd))
	if(NbCompar > 0): #More than one sample in the cluster
		for Comb in random.sample(CombInd, NbCompar):
			ind0=np.take(Geno, Comb[0], 0)
			ind1=np.take(Geno, Comb[1], 0)
			Diff=abs(ind0-ind1)#Number of difference between the two diploid. However, this is not sufficient to calculate dxy. Indeed, the number of difference between two heterozygote (genotype == "1") is 0, but in fact, their is in average 0.5 difference between the haploid genomes of two diploid heterozygote (Aa vs Aa --> could be A vs A, a vs a, A vs a, or a vs A ). So we need a fast trick to get the true number of difference. here is a way:
			ind1_No2=np.where(ind1==2,0,ind1) #Replace the homozygote derived by 0. This give a 0,1 array where 0 mean homozygote and 1 heterozygote
			ind0_No2=np.where(ind0==2,0,ind0) # Same
			maxInd_No2=np.maximum(ind1_No2, ind0_No2) #Only get position where there is one heterozygote individual (0,1 array)
			NbDiff=sum(np.maximum(Diff,maxInd_No2)/2) #Determine the maximum between this array of heterozygote indiv and the difference between the two haploids. In position where their is heterozygote individuals, this is alway 1, elsewhere it could be 2 or 0, which is what we expect. We devide by two two get the mean pairwise difference between haploid genomes, not diploid genomes
			PiSum=PiSum+NbDiff
		PiSum=PiSum/NbCompar #Mean Pi in comparison
	else: #If there is only one sample in the cluster, Pi is the number of heterozygous position
		PiSum=np.count_nonzero(Geno==1)#Number of heterozygous position
	return PiSum
		

def Compute_analyses(Array, CurrScaffold, Start, End):
	''' This is the main fonction the call the different function to perform different analyses'''
	WindowPos=[CurrScaffold,Start,End, np.shape(Array)[0]] #Position of the window 
	if np.shape(Array)[0]>4: #Verify that their is at least n loci in the window (if not, can't compute PCA with 10 axes)
		Pos=Array[:,1].astype(float) #Grep the position of the loci in the window 
		Array=Array[:,2:] #Grep only the genotype at each loci in the window (i.e. remove the loci positon and scaffold)
		Array=Array.astype(float).transpose() #Define genotye as float and transpose the matrix
		if(np.isnan(Array).any()): #No "NA" is allowed in the geno file
			print("Fatal error: The window: ", Start, "-", Start+WindSize, " contain NaN values")
			sys.exit()
		if (Method in "pca"): #If the clusteting must be perform on pca output
			ScoreList, Cluster, PC, ClusterCenterList=Define_Cluster_and_Score_WithPCA(Array,MaxCluster, ScoreFct) #Perform the PCA, the clustering and estimate the score for 2<=k<=MaxCluster
			write_pca(PC, WindowPos);
		else: #Perform the clusting directly on genotype
		#	start_time = time.time()
			ScoreList, Cluster, ClusterCenterList=Define_Cluster_and_Score_WithoutPCA(Array,MaxCluster, ScoreFct) #Perform clustering and estimate the silhouette score for 2<=k<=6
		#	print("--- %s seconds:Silhou ---" % (time.time() - start_time))
			if (optionPCA): #if pca in "options', compute pca, but don't use them to compute the cluster
				PC=Apply_PCA(Array,NoAxe)
				write_pca(PC, WindowPos);
		write_cluster(WindowPos,Cluster)
		write_clusterDistance(ClusterCenterList, WindowPos)
		#start_time = time.time()
		Compute_Heterozygosity(Array, WindowPos, Cluster)
		#print("--- %s seconds:Hete ---" % (time.time() - start_time))
		if (optionDXY):
		#	start_time = time.time()
			Compute_Dxy(Array, WindowPos, Cluster,MaxCompar)
		Line=np.concatenate((WindowPos,ScoreList))
	else:
		Line=np.concatenate((WindowPos,["NA"]*(MaxCluster - 1)))
	for element in Line:
		textfileClustScore.write(str(element) + " ")
	textfileClustScore.write("\n")

def Sliding_window_bp_overlap(File, WindSize, Slide): 
	'''
	Function to split the .geno file in window and calculate the silhoutette score or other state in each window (parameter: window size. Input: genotype array, Output: Clustering scores in each window) 
	The function read the geno file line by line, so it can deal with very large file.
	Slide on bp
	'''
	with open(File) as infile: #Read line by line (i.e. do not load the file in memory)
		next(infile)#Skip header 
		secondline=next(infile)# Get the first line of the data to store info
		Array=[secondline.split()] #create an empty array that will be used to store the genotype of individual in each window
		Array=np.asarray(Array) #Create a numpy array from the main array 
		CurrScaffold=Array[0,0]
		Start=(int(Array[0,1]) // WindSize) * WindSize # Start for sliding window (not 1 or the position of the first variant).  
		print(Start)
		for line in infile: #For each genotype position
			new_array=line.split()#Create an array from the genotype at this position (split on white space)
			if (new_array[0] == CurrScaffold and int(new_array[1]) <= (Start+WindSize)): #If the loci is on the same scaffold and within the current window
				Array=np.vstack((Array,[new_array]))
			else: #The variant not fall in the current window. It could be on another scaffold or in a next window on the same scaffold. We write the result of the analyse for this window and move to the next window
				print("Current position: Scaffold=", CurrScaffold, " Position=", Start)
				End=Start+WindSize
				#Compute_analyses(Array, CurrScaffold, Start, WindSize)
				Compute_analyses(Array, CurrScaffold, Start, End)
				if (new_array[0] == CurrScaffold): #If the loci was on the same scaffold but not on the current window
					Array=Array[Array[:,1].astype(float) > Start+Slide,] #Subset the array do remove the position not in the next window
					NumberEmpty=((int(new_array[1])-(Start+WindSize)) // Slide) #Number of empty window
					Start=Start+Slide #Change the start position of the window
					End=Start+WindSize
					for Wind in list(range(0,NumberEmpty)): #Write empty window
						Compute_analyses(Array, CurrScaffold, Start, End)
						Start=Start+Slide #Change the start position of the window
						End=Start+WindSize
						Array=Array[Array[:,1].astype(float) > Start,] #Subset the array do remove the position not in the next window
					Array=np.vstack((Array,[new_array]))#Add the new Postion to the array
				else: #The window fall in a new scaffold
					Array=[] #Reinitialise the main array (new window)
					Array.append(new_array) #push the first loci in the new array
					Array=np.asarray(Array) #Create a numpy array from the main array 
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
			#	Array=np.vstack((Array,[new_array]))
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

#### Write functions ###
def write_clusterDistance(ClusterCenterList, WindowPos):
	for NCluster in range(len(ClusterCenterList)): #For all number of cluster that might be considered (e.g. 2, 3, 4, and 5 if k was set to 5)
		distance=pairwise_distances(ClusterCenterList[NCluster]) #Calculate the pairwise euclidean distance between the cluster centroid
		ClusterComb=list(combinations(range(len(distance)),2)) #All combination of cluster
		for Comb in ClusterComb: #For each combination of cluster
			for col in WindowPos: #Print the window position
				textfileDistance.write(str(col) + " ")
			textfileDistance.write(str(len(distance)) + " " + str(Comb[0]) + " " + str(Comb[1]) + " " + str(distance[Comb[0],Comb[1]]) + "\n") #Print the distance between the centroid

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

def write_header(MaxCluster):
	textfileClustScore.write("Scaffold Start End No.variants")
	Clusters = list(range(2,MaxCluster+1))
	for cluster in Clusters:
		textfileClustScore.write(" k" + str(cluster))
	textfileClustScore.write("\n")

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
	textfileCluster.write("Scaffold Start End No.variants No.cluster")
	for ind in list(range(2,len(new_array))):
		textfileCluster.write(" " + str(new_array[ind]))
	textfileCluster.write("\n")

### Initiation : Create Output File###
if __name__ == "__main__":
	main(sys.argv[1:])
			
textfileClustScore = open(OutputFile+".ClustScore", "w")
textfileCluster = open(OutputFile+".cluster", "w")
textfileHetero = open(OutputFile+".Hetero", "w")
textfileHetero.write("Scaffold Start End No.variants No.cluster Cluster Hetero Homo")
if (optionPI):
	textfileHetero.write(" Pi\n")
else:
	textfileHetero.write("\n")

textfileDistance = open(OutputFile+".ClusterDistance", "w")
textfileDistance.write("Scaffold Start End No.variants No.cluster Cluster1 Cluster2 Distance\n")

if (optionDXY):
	textfileDxy = open(OutputFile+".Dxy", "w")
	textfileDxy.write("Scaffold Start End No.variants No.cluster Cluster1 Cluster2 Dxy\n")

if (Method in "pca" or optionPCA):
	OutputFilePCA=OutputFile+".pcaResult"
	textfilepca = open(OutputFilePCA, "w")
	write_headerPCA()

write_header(MaxCluster)
write_headerCluster()
### 

### Start Analyses ###
if (WindType in "variant"):
	Sliding_window_variant_overlap(GenoFile,WindSize, Slide)
#elif (WindType in "bp" or arg in "bp"):
elif (WindType in "bp"):
	Sliding_window_bp_overlap(GenoFile,WindSize, Slide)
else:
	print("Windows type must be 'variant' or 'bp' ")

