# SuperInfer
Script by Paul Jay, CNRS

SuperInfer is a collection of scripts to dissect the variation on genetic structure along the genome. Its first aim is to detect region of the genome that may include inversion polymorphisms or supergenes, but it can be used to detect any peculiar change in population genomic structure. For instance, it can probably highlight centremeric region and region that have been recently introgressed.

Currently, this method is based on population genetic data, e.g. resulting from short-read sequencing, and therefore detect particular genetic region based on peculiar pattern of linkage disequilibrium among group of samples in population. Consequently, it cannot be used by itself to detect with certainty a chromosomal inversion, as it does not analyse sequence orientation. However, it could provide converging lines of evidence that a genomic region may harbour a chromosomal rearrangement, which may allow to specifically focus on this region for further analyses (e.g. analysing read mapping orientation)  

Currently, the repository include 3 scripts:
 -  vcf2geno.py This script allows to transform a .vcf or .vcf.gz file into a .geno file, which is the format used by analyses
 -  GenotypeCluster.py This is the main script. It allows to perform the clustering analyses. Use python3 GenotypeCluster -h for a list of options
 -  SlidingCluster.R This script is used for plotting. Depending on the main script output, it may require some modification to get good looking plots.
	
##Usage example:
The Genotypefile.geno file should be in the format "Scaf Position GenotypeCodeSample1 GenotypeCodeSample2 ... GenotypeCodeSampleN" with genotype code being 0 for 0/0, 1 for 0/1 and 2 for 1/1 (only biallelic). All position must be genotyped (no "NA"). From a vcf file, to get it: 
 - Filter out the position with lot of missing genotype (e.g.  *bcftools view -e 'F\_MISSING > 0.2' -O z -A Genotypefile\_0.2Miss.vcf.gz GenotypeFile.vcf.gz*)
 - Impute the missing genotypes, for instance with beagle (eg. *java -Xmx20g -jar ~/Software/beagle.05May22.33a.jar  gt=Genotypefile\_0.2Miss.vcf.gz out=Genotypefile\_0.2Miss.imputed.vcf window=1 overlap=0.5*)
 - Transform to .geno file: *python3 vcf2geno.py -i Genotypefile\_0.2Miss.imputed.data.imputed.vcf -o Genotypefile\_0.2Miss.imputed.geno*

Then, use:
 > python3 ../Sex_determining/code/GenotypeCluster.py -g Genotypefile\_0.2Miss.imputed.geno -w 250 -s 25 -t variant -a pca,dxy,pi -c 100 -k 5 -o Genotypefile\_0.2Miss.W250.S25 
to perform all possible analyses: Kmeans cluster, Silhouette score, PCA, Cluster centroid distance, Dxy, Pi, etc. 
