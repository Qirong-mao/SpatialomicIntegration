## Normalize the output of Cell2location (Cell type density of each spot) to cell type probablity of each spot.
dis<-read.csv("distance.csv",header=TRUE)
dis<-dis[,-1]
rownames(dis)<-1:4972
colnames(dis)<-1:4972
prob_novosparc<-read.csv("CellProb.csv",header=TRUE)
prob_cell2location<-read.csv("W_cell_density_q05.csv",header = TRUE)
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
rm(i)
rm(j)
rm(dist)
rm(sum)
rm(k)
rm(l)
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
rm(dist)
rm(sum)
rm(k)
rm(l)
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
  