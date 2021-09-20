library(tidyverse)
library(ggplot2)

## Normalize the output of Cell2location (Cell type density of each spot) to cell type probablity of each spot.

all_means = read.table("means.txt", header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F) #L-R expression matrix from CellPhoneDB output
expression<-as.data.frame(t(all_means[,12:3855]))
expression$cluster<-rownames(expression)
dis<-read.csv("distance.csv",header=TRUE) #Euclidean distance matrix between spots in ST data.
dis<-dis[,-1]
rownames(dis)<-1:4972
colnames(dis)<-1:4972
prob_novosparc<-read.csv("CellProb.csv",header=TRUE) #Probability matrix from Novosparc
prob_cell2location<-read.csv("W_cell_density_q05.csv",header = TRUE) # Cell type density matrix from Cell2location
colnames(prob_cell2location)=sub('^................','',colnames(prob_cell2location))
prob_novosparc<-prob_novosparc[,-1]
prob_novosparc<-t(prob_novosparc)
colnames(prob_novosparc)<-0:61
rownames(prob_novosparc)<-1:4972
prob_cell2location<-prob_cell2location[,-1]
prob_cell2location<-as.data.frame(t(prob_cell2location))
prob_cell2location$cluster<-as.numeric(rownames(prob_cell2location))
prob_cell2location<-prob_cell2location[order(prob_cell2location$cluster),]
prob_cell2location<-prob_cell2location[,-4973]
prob_cell2location<-t(prob_cell2location)
rownames(prob_cell2location)<-1:4972
prob_cell2location<-sweep(prob_cell2location,2,colSums(prob_cell2location),"/")
colnames(prob_cell2location)<-0:61
rownames(prob_cell2location)<-1:4972
df_novosparc<-data.frame(Distance=seq(1:3844),Group1=seq(1:3844),Group2=seq(1:3844))

##For each cell type, we remain the spots with top 1% probrabilities (50 spots, >=0.99)
for (i in 1:62) {
  prob_novosparc[,i]<-as.numeric(prob_novosparc[,i])
  prob_novosparc[,i][prob_novosparc[,i]<quantile(prob_novosparc[,i],probs = 0.99)]<-NA
}
for (i in 1:62) {
  prob_cell2location[,i]<-as.numeric(prob_cell2location[,i])
  prob_cell2location[,i][prob_cell2location[,i]<quantile(prob_cell2location[,i],probs = 0.99)]<-NA
}

## Caculating the average distance between cell type clusters.
## Novosparc
index=1
for (i in 1:62) {
  a<-c(as.integer(names(na.omit(prob_novosparc[,i]))))
  for (j in 1:62) {
    b<-c(as.integer(names(na.omit(prob_novosparc[,j]))))
    list<-expand.grid(a,b)
    sum=sum(dis[unique(list$Var1),unique(list$Var2)])/2500 # 2500 = 50 cell types * 50 cell types
    df_novosparc$Distance[index]=sum
    df_novosparc$Group1[index]=i
    df_novosparc$Group2[index]=j
    index=index+1
  }
}
## Cell2location
index=1
rm(i)
rm(j)
rm(sum)
df_cell2location<-data.frame(Distance=seq(1:3844),Group1=seq(1:3844),Group2=seq(1:3844))
for (i in 1:62) {
  a<-c(as.integer(names(na.omit(prob_cell2location[,i]))))
  for (j in 1:62) {
    dist = 0
    sum = 0
    b<-c(as.integer(names(na.omit(prob_cell2location[,j]))))
    list<-expand.grid(a,b)
    sum=sum(dis[unique(list$Var1),unique(list$Var2)])/2500
    df_cell2location$Distance[index]=sum
    df_cell2location$Group1[index]=i
    df_cell2location$Group2[index]=j
    index=index+1
  }
}

  
## Caculate the correlation between the L-R expression and the distance

df_cell2location$Group1<-as.character(df_cell2location$Group1-1)
df_cell2location$Group2<-as.character(df_cell2location$Group2-1)
df_novosparc$Group1<-as.character(df_novosparc$Group1-1)
df_novosparc$Group2<-as.character(df_novosparc$Group2-1)
df_cell2location$cluster<-as.character(paste(df_cell2location[,2],df_cell2location[,3],sep="|"))
df_novosparc$cluster<-as.character(paste(df_novosparc[,2],df_novosparc[,3],sep="|"))
corr_cell2location<-merge(df_cell2location,expression,by="cluster")
corr_novosparc<-merge(df_novosparc,expression,by="cluster")
corr_cell2location<-cbind(corr_cell2location[,2],corr_cell2location[,5:1351])
corr_novosparc<-cbind(corr_novosparc[,2],corr_novosparc[,5:1351])
colnames(corr_cell2location)[1]<-"Distance"
colnames(corr_novosparc)[1]<-"Distance"
## For each L-R pair, we used the average expression value in each distance group
corr_cell2location<-corr_cell2location %>% group_by(Distance) %>% summarise(across(everything(), list(mean))) %>% ungroup
corr_novosparc<-corr_novosparc %>% group_by(Distance) %>% summarise(across(everything(), list(mean)))%>%ungroup()
correlation_cell2location<-cor(corr_cell2location[,1],corr_cell2location[,2:1348],method = "spearman")
correlation_novosparc<-cor(corr_novosparc[,1],corr_novosparc[,2:1348],method = "spearman")
mutual<-as.data.frame(cbind(t(correlation_cell2location),t(correlation_novosparc)))
colnames(mutual)<-c("corr_cell2location","corr_novosparc")
# Calculate the quantile of correlation
mutual$quantile_cell2location<-ecdf(mutual$corr_cell2location)(mutual$corr_cell2location)
mutual$quantile_novosparc<-ecdf(mutual$corr_novosparc)(mutual$corr_novosparc)
## Calculate the FDR value for each L-R pair based on the quantile of correlation
## Celll2location
for (i in 1:1347) {
  if (mutual$quantile_cell2location[i]<=0.5){
    mutual$pvalue_cell2location[i]=as.numeric(2*mutual$quantile_cell2location[i])
  } else {
    mutual$pvalue_cell2location[i]=as.numeric(2*(1-mutual$quantile_cell2location[i]))
  }
}
mutual$fdr_cell2location<-p.adjust(mutual$pvalue_cell2location,method = "BH")
## Novosparc
for (i in 1:1347) {
  if (mutual$quantile_novosparc[i]<=0.5){
    mutual$pvalue_novosparc[i]=as.numeric(2*mutual$quantile_novosparc[i])
  } else {
    mutual$pvalue_novosparc[i]=as.numeric(2*(1-mutual$quantile_novosparc[i]))
  }
}
mutual$fdr_novosparc<-p.adjust(mutual$pvalue_novosparc,method = "BH")
mutual$label<-factor(result1$label,levels = c("Mutual Significant","Cell2location Significant","Novosparc Significant","Not Significant"))