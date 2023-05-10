library(ggplot2, warn.conflicts = F, quietly = T)
library(cowplot, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(tidyr, warn.conflicts = F, quietly = T)
options(dplyr.summarise.inform = FALSE)
args = commandArgs(trailingOnly=TRUE)
filepca=args[1]
filepca="~/Projects/SuperInfer/Human/1000Gen/PJL.w750.s75.p3.k3"
pca=read.table(paste0(filepca, ".pcaResult"), header=T, stringsAsFactors = F)
Cluster=read.table(paste0(filepca, ".cluster"), header=T, stringsAsFactors = F)
pcaLong=pca[pca$PC==1,] %>% pivot_longer(cols=!c(Scaffold,Start,End,No.variants,PC), names_to="Ind", values_to="value") #Make long data
Cluster3=Cluster[(Cluster$No.cluster==3 & Cluster$Axes==1),] #Keep only cluster data for first axes and with three cluster
Cluster3Long=Cluster3 %>% pivot_longer(cols=!c(Scaffold,Start,End,No.variants,No.cluster, Axes), names_to="Ind", values_to="value")
pcaLong$Cluster=Cluster3Long$value #Add the cluster value to the pca data


pcaLong=pcaLong %>% group_by(Start, End, Cluster) %>% mutate(centroid=mean(value)) %>% mutate(distance=(value-centroid)^2) #Get the cluster centroid position of the cluster of each individuals and determine the distance of this sample to this centroid (distance from own cluster)
pcaLong=pcaLong %>% group_by(Start, End) %>% mutate(top=max(centroid), bottom=min(centroid), mid=median(centroid)) #Get the centroid position of each cluster
pcaLong=pcaLong %>% group_by(Start, End, Ind) %>% mutate(ClosestDist=min(c((value-top)^2, (value-bottom)^2, (value-mid)^2)), MidDist=median(c((value-top)^2, (value-bottom)^2, (value-mid)^2)), maxDist=max(c((value-top)^2, (value-bottom)^2, (value-mid)^2))) # determine the distance of each sample to each cluster (#probably ugly and unifficient). The shortest distanve correspond to their own cluster (so previous step useless...) and the mid distance to the closest other cluster
pcaLongSum=pcaLong %>% group_by(Start, End) %>% summarise(dispersion=sqrt(mean(distance)), DistDisp=sqrt(mean(MidDist))) #for each position, determine the mean distance to the oqn cluster and the mean distance to the closest other cluster
pcaLongSum$ClustScore=pcaLongSum$DistDisp/pcaLongSum$dispersion #Caculate a cluster score, which is the dispersion to their own cluster compared to the distance to the closest cluster.
SignifLineScore=quantile(pcaLongSum$ClustScore, 0.95) #Signif value

pcaLong = pcaLong %>% rowwise %>% mutate(NewColor=ifelse(centroid==bottom, 0, ifelse(centroid==top, 2, 1))) #redefine cluster identity based on position 

Interest=pcaLongSum[pcaLongSum$ClustScore>SignifLineScore,] #Extract the region with score above limit
Interest$Follow=ifelse(lag(Interest$End)>Interest$Start, 0,1) #overlapping windows are indicated with 0
Interest$Follow[1]=0
Interest$ID= cumsum(Interest$Follow) #each time we found non overlappîng window, this ID increase. SO this ID identify different region.
InterestSum=Interest %>% group_by(ID) %>% summarise(MinPos=min(Start), MaxPos=max(End)) # Identify start and end of each regions.

df= data.frame(matrix(nrow = 0, ncol = length(unique(pcaLong$Ind))+3))
colnames(df)=c("Chrom", "Start", "End", unique(pcaLong$Ind))
for (Inv in 1:nrow(InterestSum))
{
    pcaSub=pcaLong[(pcaLong$Start>=InterestSum$MinPos[Inv] & pcaLong$End<=InterestSum$MaxPos[Inv]),] #Subset to get only the intersting pos
    GenoInv=pcaSub %>% group_by(Ind) %>% summarise(Geno=round(mean(NewColor))) %>% pivot_wider(names_from = Ind, values_from = Geno) #Extract the cluster (i.e. phenotype) and average over the region.
    GenoInv=cbind(as.character(pcaSub[1,1]), InterestSum$MinPos[Inv], InterestSum$MaxPos[Inv], GenoInv)
    df[nrow(df)+1,]=GenoInv[1,]  #write.table(GenoInv,paste0(filepca, ".", pcaSub$Scaffold[1], "_", InterestSum$MinPos[Inv], "-",InterestSum$MaxPos[Inv],".V2.InvPhenotype"),  quote=F, row.names = F)
}

write.table(df,paste0(filepca,".V3.InvPhenotype"),  quote=F, row.names = F)
