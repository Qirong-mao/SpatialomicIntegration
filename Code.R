library(tidyverse)
library(ggplot2)
library(ggrepel)
library(philentropy)
library(ggpubr)
library(ggExtra)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(gridGraphics)
set.seed(62)

## Normalize the output of Cell2location (Cell type density of each spot) to cell type probability of each spot.

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
correlation_all<-as.data.frame(cbind(t(correlation_cell2location),t(correlation_novosparc)))
colnames(correlation_all)<-c("corr_cell2location","corr_novosparc")



## Creating the random distance matrix


## 1. Randomly select 50 spots for each cell type clusters(100 times) and calculate the distance.
random_distance<-as.data.frame(matrix(nrow=3844,ncol=1000))
## You can input the random_distance matrix here to reproduce the result in the thesis, since generating this
## random_distance matrix is a random process, the result might be different.
df<-as.data.frame(matrix(nrow=50,ncol=62))
for (q in 1:1000){
  for (i in 1:62) {
    df[,i]<-sample(1:4971,50)
  }
index=1
for (i in 1:62) {
  a<-c(as.integer(df[,i]))
  for (j in 1:62) {
    sum = 0
    b<-c(as.integer(df[,j]))
    list<-expand.grid(a,b)
    sum=sum(dis[unique(list$Var1),unique(list$Var2)])/2500
    random_distance[index,q]=sum
    index=index+1
  }
}
cat("Permutation Round: ", q,"\n")
}
random_correlation<-as.data.frame(matrix(nrow = 1000,ncol=1347))
for (i in 1:1000){
  random_correlation[i,]<-cor(random_distance1[,i],expression[,1:1347],method = "spearman")
}
# Calculate the quantile of correlation
for (i in 1:1347) {
correlation_all$quantile_cell2location[i]<-ecdf(random_correlation[,i])(correlation_all$corr_cell2location[i])
correlation_all$quantile_novosparc[i]<-ecdf(random_correlation[,i])(correlation_all$corr_novosparc[i])
}
## Calculate the FDR value for each L-R pair based on the quantile of correlation
## Celll2location
for (i in 1:1347) {
  if (correlation_all$quantile_cell2location[i]<=0.5){
    correlation_all$pvalue_cell2location[i]=as.numeric(2*correlation_all$quantile_cell2location[i])
  } else {
    correlation_all$pvalue_cell2location[i]=as.numeric(2*(1-correlation_all$quantile_cell2location[i]))
  }
}
correlation_all$fdr_cell2location<-p.adjust(correlation_all$pvalue_cell2location,method = "BH")
## Novosparc
for (i in 1:1347) {
  if (correlation_all$quantile_novosparc[i]<=0.5){
    correlation_all$pvalue_novosparc[i]=as.numeric(2*correlation_all$quantile_novosparc[i])
  } else {
    correlation_all$pvalue_novosparc[i]=as.numeric(2*(1-correlation_all$quantile_novosparc[i]))
  }
}
correlation_all$fdr_novosparc<-p.adjust(correlation_all$pvalue_novosparc,method = "BH")
## Label
for (i in 1:1347) {
  if (correlation_all$fdr_novosparc[i]<=0.05 && correlation_all$fdr_cell2location[i]<=0.05) {
    correlation_all$label[i]<-"correlation_all Significant"
  }
  else if (correlation_all$fdr_cell2location[i]<=0.05){
    correlation_all$label[i]<-"Cell2location Significant"
  }
  else if (correlation_all$fdr_novosparc[i]<=0.05){
    correlation_all$label[i]<-"Novosparc Significant"
  }
  else {
    correlation_all$label[i]<-"Not Significant"
  }
}
correlation_all$label<-factor(correlation_all$label,levels = c("correlation_all Significant","Cell2location Significant","Novosparc Significant","Not Significant"))
correlation_all<-cbind(correlation_all,all_means[,(2:6)])

##=========================Figures========================================

## Fig.3 J-S Divergence Plot
Prob_sum<-cbind(prob_cell2location,prob_novosparc)
JSD.sum<-JSD(t(Prob_sum))
annotation_col = data.frame(CellType = factor(rep(c("Cell2location","Novosparc"), c(62,62))))
colnames(JSD.sum)[1:62]<-sprintf("novosparc[%s]",seq(1:62))
colnames(JSD.sum)[63:124]=sprintf("cell2location[%s]",seq(1:62))
rownames(JSD.sum)<-colnames(JSD.sum)
rownames(annotation_col)<-rownames(JSD.sum)
pheatmap(JSD.sum,annotation_col = annotation_col,annotation_row = annotation_col,legend = TRUE,show_rownames=FALSE, show_colnames=FALSE,cellwidth = 5,cellheight = 5,legend_breaks = c(0.9, 0.1),legend_labels = c("High J-S Divergence","Low J-S Divergence"),annotation_names_row=FALSE,annotation_names_col=FALSE)

## Fig.4 Correlation plot
scatter=ggplot(correlation_all,aes(x=corr_cell2location,y=corr_novosparc,color=label))+geom_point()+geom_text_repel(aes(label=interacting_pair),max.overlaps = 16)+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA))+theme(legend.title=element_blank(),legend.position = c(0.9,0.1))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+labs(x="Correlation (Cell2location)",y="Correlation (Novosparc)")+theme(axis.title = element_text(face="bold"),legend.text = element_text(face="bold"))
ggMarginal(scatter,type="histogram",color="black",fill="#00AFBB")

## Fig.5 & 6 Gene Ontology plot
genelist_all<-(unique(c(correlation_all$gene_a,correlation_all$gene_b)))
genelist_up<-(unique(c(correlation_up$gene_a,correlation_up$gene_b)))
genelist_down<-(unique(c(correlation_down$gene_a,correlation_down$gene_b)))

up<-enrichGO(na.omit(genelist_up),keyType = "SYMBOL",OrgDb = org.Hs.eg.db,ont = "ALL",pvalueCutoff = 0.05,pAdjustMethod = "BH",universe = na.omit(genelist_all))
dotplot(up, orderBy = "GeneRatio", color = "p.adjust", showCategory = 15, size = NULL,font.size = 10,title="Gene Ontology of Positively distance-correlated L-R pairs", split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")

down<-enrichGO(na.omit(genelist_down),keyType = "SYMBOL",OrgDb = org.Hs.eg.db,ont = "ALL",pvalueCutoff = 0.05,pAdjustMethod = "BH",universe = na.omit(genelist_all))
dotplot(down, orderBy = "GeneRatio", color = "p.adjust", showCategory = 15, size = NULL,font.size = 10,title="Gene Ontology of Negatively distance-correlated L-R pairs", split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")

## Fig.7 Composition Plot

labelmatrix<-read.csv("LabelMatrix.csv",header=TRUE)
pheatmap(t(labelmatrix),cluster_rows = FALSE,cluster_cols = FALSE,fontsize = 12,scale="column",color = colorRampPalette(c("white","red"))(5),cellwidth = 12,cellheight = 12,breaks = seq(0,5,by=1))

## Fig.8 Significant expression between cluster-cluster plot
down_matrix<-read.csv("down_matrix.csv")
LR<-as.data.frame(matrix(as.numeric(down_matrix[3,-(1:13)]),nrow = 62,ncol = 62))
pheatmap(LR,fontsize = 12,color = brewer.pal(8,"Blues"),cellwidth = 10,cellheight = 10,main="WNT5A-FZD3")
grid.text("Receiver", y=0.05, gp=gpar(fontsize=16))
grid.text("Sender", x=0.03, rot=90, gp=gpar(fontsize=16))