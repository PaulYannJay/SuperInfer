# SuperInfer
Scripts by Paul Jay

SuperInfer is a collection of scripts to dissect the variation on genetic structure along the genome. Its first aim is to detect region of the genome that may include inversion polymorphisms or supergenes, but it can be used to detect any peculiar change in population genomic structure. For instance, it can probably highlight centromeric regions and regions that have been recently introgressed.

Currently, this method is based on population genetic data, e.g. resulting from short-read sequencing, and therefore detects genetic regions with particular genetic structure based on peculiar pattern of linkage disequilibrium among group of samples in population. Consequently, it cannot be used by itself to detect with certainty a chromosomal inversion, as it does not analyse sequence orientation. However, it could provide converging lines of evidence that a genomic region may harbour a chromosomal rearrangement, which may allow to specifically focus on this region for further analyses (e.g. analysing read mapping orientation)  

Currently, the repository include three scripts:
 -  vcf2geno.py This script allows to transform a .vcf or .vcf.gz file into a .geno file, which is the format used by analyses
 -  GenotypeCluster.py This is the main script. It allows to perform the clustering analyses. Use python3 GenotypeCluster -h for a list of options
 -  SlidingCluster.R This script is used for plotting. Depending on the main script output, it may require some modification to get good looking plots.
	
##Usage example

The Genotypefile.geno file should be in the format "Scaffold Variant\_Position GenotypeCodeSample1 GenotypeCodeSample2 ... GenotypeCodeSampleN" with genotype codes being 0 for 0/0, 1 for 0/1 and 2 for 1/1 (only biallelic). All position must be genotyped (no "./." or "NA"). From a vcf file, to get it: 
 - Filter out the position with lot of missing genotypes (e.g.  *bcftools view -e 'F\_MISSING > 0.2' -O z -A Genotypefile\_0.2Miss.vcf.gz GenotypeFile.vcf.gz*)
 - Impute the missing genotypes, for instance with beagle (eg. *java -Xmx20g -jar ~/Software/beagle.05May22.33a.jar  gt=Genotypefile\_0.2Miss.vcf.gz out=Genotypefile\_0.2Miss.imputed.vcf window=1 overlap=0.5*)
 - Transform to .geno file: *python3 vcf2geno.py -i Genotypefile\_0.2Miss.imputed.data.imputed.vcf -o Genotypefile\_0.2Miss.imputed.geno*

Then, use:
 > python3 ../Sex_determining/code/GenotypeCluster.py -g Genotypefile\_0.2Miss.imputed.geno -w 250 -s 25 -t variant -a pca,dxy,pi -c 100 -k 5 -o Genotypefile\_0.2Miss.W250.S25 

to perform all possible analyses: Kmeans cluster, Silhouette score, PCA, Cluster centroid distance, Dxy, Pi, etc. In this example, the script use 250 variant sliding window with 225 overlap between windows (250-25), compute the analyses considering up to 5 clusters, and perform a maximum of 100 comparisons for calculating Pi and Dxy.
E

Note: the main script determine the clusters for a number of cluster up to k. So, if -k 5 is provided, the script splits the samples in two groups, then computes the analyses (e.g. heterozygosity) separately on these two groups of samples and outputs the result, then creates three groups, computes the analyses on these groups and outputs the result, and then do the same for four groups, and finally for five groups. The output files therefore contains a column "No.cluster", which indicate the number of clusters considered by the script for these analyses  (e.g. 2, 3, 4 or 5 if -k 5 is provided), and a column "cluster", which indicate the cluster on which analyses were performed (e.g. if the column "no.cluster" indicate "4", the column "cluster" can be 0, 1, 2 or 3).

Outputs are:
 - OutputFile.Silhou --> The silhouette score for each k, for each window 
 - OutputFile.Hetero --> The proportion of heterozygous and homozygous site (only variant) in each cluster, for each number of cluster, for each window. If "-a Pi" is provided, this file also contains the result of Pi calculation for each cluster.
 - OutputFile.pcaResult --> The position of each sample on each pca axes (by default up to 10), for each window
 - OutputFile.cluster --> The affiliation of each sample (the cluster they belong to), for each number of cluster, for each window
 - OutputFile.clusterDistance --> The euclidian distance between each cluster, for each number of cluster, for each window
 - OutputFile.Dxy --> The nucleotide distance between each cluster, for each number of cluster, for each window

Plots can be obtained with : *Rscript SlidingCluster.R Genotypefile\_0.2Miss.W250.S25*
The plots are: 
 - a. The position of sample on the first pca axis, after a bit of rescaling. This allows to visualise the variation on genetic structure
 - b. The Silhouette score for each value of k. At a supergene or an inversion polymorphism, the silhouette is best for k=2 or k=3
 - c. For k=3, the ratio of the heterozygosity of the cluster with the highest heterozygosity over the heterozygosity of the cluster with the lowest heterozygosity. Since inversion heterozygotes are expected to have a very high heterozygosity, in contrast to the inversion or ancestral homozygotes, peak of this statistic are likely to represent region with putative inversion polymorphism.
 - d. The Dxy of the most distant clusters for k=3. This is expected to be high in non-recombining regions (e.g. inversion polymorphism), and low in recombining regions
 - e. The distance between the centroid of the most distant clusters for k=3. This is similar to Dxy but much quicker to compute. This is expected to be high in non-recombining regions (e.g. inversion polymorphism), and low in recombining regions

Enjoy !
