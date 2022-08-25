library(ggplot2)
library(cowplot)
library(viridis)
library(dplyr)
library(viridis)
library(tidyverse)

ThemeSobr=  theme(
  panel.border = element_blank(),  
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  text = element_text(size=12, face="bold"),
  axis.line = element_line(colour = "black", size=0.5),
  legend.spacing.y= unit(0, 'cm'),
  axis.title = element_text(face="bold"),
  plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=2),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  axis.text = element_text(size = 10),
  legend.background = element_blank(),
  panel.grid.major.x = element_line(size = 0.2, color="grey50", linetype="dashed"),
  panel.grid.minor.x = element_line(size = 0.1, color="grey50", linetype="dashed"),
  strip.text = element_text(face="bold", size=14)
)

#  
Col=c("#dcad34","#bf0c14","#793e74","#11171c","#2a5d9e") #Cycliste sous la pluie, Stuart Franklin
args = commandArgs(trailingOnly=TRUE)
file=args[1]
#PCA
pca=read.table(paste0(file,".pcaResult"), stringsAsFactors = F, header=T)
# pcaLong=pca %>% pivot_longer(cols=starts_with("Ver"), names_to="Ind", values_to="value")
pcaLong=pca %>% pivot_longer(cols=!c(Scaffold,Start,End,No.variants,PC), names_to="Ind", values_to="value")

pcaLong1=subset(pcaLong, pcaLong$PC==1)
#Rotate the PCA# Not usefull most of the time !
# pcaLong1First=pcaLong1[(pcaLong1$Scaffold==pcaLong1$Scaffold[1] &
#                           pcaLong1$Start==pcaLong1$Start[1] &
#                           pcaLong1$End==pcaLong1$End[1]),]
# pcaLong1Lag=pcaLong1 %>% group_by(Scaffold, Start, End) %>% mutate(FirstVal=pcaLong1First$value)
# pcaLong1LagDist=pcaLong1Lag %>% group_by(Scaffold, Start, End) %>% summarise(StandDist=mean(abs(value-FirstVal)), StandDistNeg=mean(abs(value+FirstVal)))
# pcaLong1LagDist= pcaLong1LagDist %>% mutate(WindowName= paste0(Scaffold,Start,End))
# pcaLong1LagDist= pcaLong1LagDist %>% mutate(Reorder=(if_else(StandDist>StandDistNeg, T, F)))
# WindowToInvert=pcaLong1LagDist[pcaLong1LagDist$Reorder==T,]$WindowName
# pcaLong1=pcaLong1 %>% mutate(WindowName= paste0(Scaffold,Start,End))  #Create a column idenfying the window
# pcaLong1[pcaLong1$WindowName %in% WindowToInvert,]$value=-pcaLong1[pcaLong1$WindowName %in% WindowToInvert,]$value
print("Reading and analysing PCA output: Done")

#Sihouete Score
Silhou=read.table(paste0(file,".ClustScore"), stringsAsFactors = F, header=T)
SilhouLong=Silhou %>% group_by(Scaffold, Start, End) %>% pivot_longer(cols = starts_with("k"), names_to="k", values_to="Clust_score")
SilhouLongSub=subset(SilhouLong, (SilhouLong$No.variants>5))
if(Silhou$Score[1]=="Silhouette")
{
	SilhouetteBest=SilhouLongSub %>% group_by(Scaffold,Start,End) %>% slice_max(Clust_score) #Determine the best k in each window
} else if(Silhou$Score[1]=="Davies_bouldin")
{
	SilhouetteBest=SilhouLongSub %>% group_by(Scaffold,Start,End) %>% slice_min(Clust_score) #Determine the best k in each window
}
meanSize=mean(SilhouetteBest$End - SilhouetteBest$Start)
meanSlide=mean(SilhouetteBest$Start - lag(SilhouetteBest$Start), na.rm=T)
FracOverlap=meanSlide/meanSize

print("Reading and analysing ClustScore output: Done")

#Heterozygosity
Het=read.table(paste0(file, ".Hetero"), stringsAsFactors = F, header=T)
Het3=subset(Het, Het$No.cluster==3)
HetMax=Het3 %>% group_by(Scaffold, Start, End) %>% summarise(maxHet=max(Hetero), minHet=min(Hetero),maxHom=max(Homo), minHom=min(Homo))
HetMax$minmaxHet=HetMax$maxHet/HetMax$minHet
Quant95=quantile(HetMax$minmaxHet, 0.95)
HetMaxSup95=HetMax[HetMax$minmaxHet > Quant95,]
print("Reading and analysing Heterozygosity output: Done")

#Dxy
if(file.exists(paste0(file,".Dxy")))
	{
	DxyData=read.table(paste0(file,".Dxy"), stringsAsFactors = F, header=T,fill=T)
	DxyMax=DxyData %>% group_by(Scaffold, Start, End, No.cluster) %>% summarise(MaxDist=max(Dxy))
	DxyMax3Clus=subset(DxyMax, DxyMax$No.cluster==3)
	print("Reading and analysing Dxy output: Done")
}

#Cluster Distance
DistanceData=read.table(paste0(file, ".ClusterDistance"), stringsAsFactors = F, header=T, fill=T)
DistanceMax=DistanceData %>% group_by(Scaffold, Start, End, No.cluster) %>% summarise(MaxDist=max(Distance))
DistanceMax3Clus=subset(DistanceMax, DistanceMax$No.cluster==3)
print("Reading and analysing Distance output: Done")

for (Scaff in unique(Het$Scaffold)) ### Not tested avec Distance ### tester et inclure Dxy
{
  
  base=ggplot(pcaLong1[pcaLong1$Scaffold==Scaff,])
  PCAplot=base+
    geom_line(aes(x=Start+(End-Start)/2, y=value, color=Ind), alpha=0.3)+
    xlab("Position on chr")+ylab("PC1")+
    theme(legend.position = "none")+
    theme(panel.grid.major = element_line(size = 1))+
    ThemeSobr
  
#print("PCAplot: Done")

if (nrow(HetMaxSup95[HetMaxSup95$Scaffold==Scaff,])>0)
{
  base=ggplot(HetMax[HetMax$Scaffold==Scaff,])
  HETplot=base+
    geom_hline(yintercept = Quant95, color="red3", alpha=1,linetype="dashed", size=0.1)+
    geom_line(aes(x=Start, y=minmaxHet), color=Col[1])+
	geom_point(data=HetMaxSup95[HetMaxSup95$Scaffold==Scaff,], aes(x=Start, y=-quantile(HetMax$minmaxHet,0.5)), color="red3", shape=8, size=0.2)+
    xlab("Position on chr")+ylab("MaxHet/MinHet")+
    ThemeSobr
 save_plot(paste0(file,".",Scaff,".png"),HETplot, nrow = 5, base_aspect_ratio = 5)
 }
 else
 {
   base=ggplot(HetMax[HetMax$Scaffold==Scaff,])
  HETplot=base+
    geom_hline(yintercept = Quant95, color="red3", alpha=1,linetype="dashed", size=0.1)+
    geom_line(aes(x=Start, y=minmaxHet), color=Col[1])+
    xlab("Position on chr")+ylab("MaxHet/MinHet")+
    ThemeSobr
 save_plot(paste0(file,".",Scaff,".png"),HETplot, nrow = 5, base_aspect_ratio = 5)
 }

#print("Heteroplot: Done")
 base=ggplot(SilhouLongSub[SilhouLongSub$Scaffold==Scaff,])
  SILHOUplot=base+
    geom_line(aes(x=Start+(End-Start)/2, y=Clust_score, color=k))+
    scale_color_manual(values=Col)+ 
    geom_rect(data=SilhouetteBest[SilhouetteBest$Scaffold==Scaff,], 
              aes(xmin=Start, xmax=End, 
                  ymin=max(SilhouLongSub[SilhouLongSub$Scaffold==Scaff,]$Clust_score)+
                    0.05*max(SilhouLongSub[SilhouLongSub$Scaffold==Scaff,]$Clust_score),
                  ymax=max(SilhouLongSub[SilhouLongSub$Scaffold==Scaff,]$Clust_score)+
                    0.10*max(SilhouLongSub[SilhouLongSub$Scaffold==Scaff,]$Clust_score),
                  fill=k), alpha=FracOverlap)+
    scale_fill_manual(values=Col)+ 
    xlab("Position on chr")+ylab("Clustering score")+
    ThemeSobr
#print("Scoreplot: Done")
  
  base=ggplot(DistanceMax3Clus[DistanceMax3Clus$Scaffold==Scaff,])
  DISTPlot= base+
    geom_line(aes(x=Start+(End-Start)/2, y=MaxDist), color=Col[2])+
    xlab("Position on chr")+ylab("Distance")+
    theme(legend.position = "none")+
    ThemeSobr
  
#print("Distanceplot: Done")

  if(file.exists(paste0(file,".Dxy")))
	{ 
	base=ggplot(dxymax3clus[dxymax3clus$scaffold==scaff,])
  	dxyplot= base+
    geom_line(aes(x=Start+(End-Start)/2, y=MaxDist), color=Col[3])+
    xlab("Position on chr")+ylab("Dxy")+
    theme(legend.position = "none")+
    ThemeSobr
	print("Dxyplot: Done")
	}
  
  if(file.exists(paste0(file,".Dxy")))
	{
 	plots <- align_plots(PCAplot,SILHOUplot,HETplot,DXYPlot,DISTPlot, align = 'vh', axis = 'lrtb')
  	PLOTS=plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],nrow = 5, labels = "auto")
    save_plot(paste0(file,".",Scaff,".png"),PLOTS, nrow = 5, base_aspect_ratio = 5)
	}
  else
	{
 	plots <- align_plots(PCAplot,SILHOUplot,HETplot,DISTPlot, align = 'vh', axis = 'lrtb')
  	PLOTS=plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]],nrow = 4, labels = "auto")
    save_plot(paste0(file,".",Scaff,".png"),PLOTS, nrow = 4, base_aspect_ratio = 5)
	}
}
