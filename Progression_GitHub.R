library('clusterSim')
library(Rtsne)
library("fpc")
library("dbscan")
library(dplyr)
library(ggpubr)
library(factoextra)
library('proxy')
library(ggpubr)
library("cluster")
library(philentropy)
library(ClusterR)

##############half year progression analysis################
all=read.csv("Newchange_4-8_months.csv")
colnames(all)
table(all$group)
case=all[all$group=="Case",]
control=all[all$group!="Case",]

####1.1-mean imputation#####
colnames(all)[1:10]
e=c()
for (i in 4:799){
  a=sum(is.na(all[,i]))
  if (a>0){
    e=c(e,i)
  }
}
#380
all1=all
for(i in 4:799){
  all1[is.na(all1[,i]), i]=round(mean(all1[,i], na.rm = TRUE),2)
}

e=c()
for (i in 4:799){
  a=sum(is.na(all1[,i]))
  if (a>0){
    e=c(e,i)
  }
}
#0
#all1$sex=as.numeric(all1$sex)

all2=all1
colnames(all2)

all2$summary_km_f
table(all2$group)
all3=all2[-which(all2$summary_km_f<0),]
table(all3$group)
summary(all3$summary_km_f)
table(all3$group)
summary(all3$summary_km_f)
all4=all3
summary(all4$follow_up)

#all4=all4[,-93]
all5=all4


library(ggpubr)
library(factoextra)
# Dimension reduction using PCA
colnames(all5)[1:10]
res.pca <- prcomp(all5[, 5:799],  scale = TRUE)
a=res.pca$rotation[,1]
which.max(a)


# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
# Add Species groups from the original data sett
ind.coord$group <- all5$group
# Data inspection
head(ind.coord)

# Percentage of variance explained by dimensions
eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
getOption("max.print")
eigenvalue

ind.coord[ind.coord$group=="Sub",]$group="Case"
library(lsa)

###three cluster
# output=data.frame("No_of_PC"=0,"No_of_Clusters"=0,"purity"=0,stringsAsFactors=T)
# min=800
# for(n in 2:10){
#   for (i in 1:406){
#     dist_mat <- dist(ind.coord[,1:i], method = "euclidean")
#     set.seed(23)
#     hclust_avg <- hclust(dist_mat, method = "ward.D")
#     cut_avg <- cutree(hclust_avg, k =n)
#     a=as.data.frame(table(ind.coord$group,cut_avg))
#     b=n
#     p=0
#     for(j in 1:b){
#       k=(max(a[a$cut_avg==j,]$Freq))/sum(a[a$cut_avg==j,]$Freq)
#       p=p+k
#     }
#     h=round(p/b,5)
#     output <- rbind(output, data.frame("No_of_PC"=i,"No_of_Clusters"=n,"purity"=h))
#   }
#   #cluster 2:78 euclidean
#   
# }
# output1=output[-1,]
# output2=output1[output1$No_of_Clusters<5,]
# output$No_of_Clusters
# ggplot(data=output2, aes(x=No_of_PC, y=purity, group=No_of_Clusters)) +
#   geom_line(aes(linetype=factor(No_of_Clusters)))
# 
# ###two cluster
# min=800
# for (i in 100:1){
#   dist_mat <- dist(ind.coord[,1:i], method = "euclidean")
#   set.seed(23)
#   hclust_avg <- hclust(dist_mat, method = "ward.D")
#   cut_avg <- cutree(hclust_avg, k =3)
#   a=table(ind.coord$group,cut_avg)
#   print(i)
#   print(a)
#   # if((a[2,2])<min){
#   #   min=a[2,2]
#   #   print(i)
#   #   print(a)
#   # }
# }
# 
# 
# 
# ###four cluster
# min=800
# for (i in 200:100){
#   dist_mat <- dist(ind.coord[,1:i], method = "cosine")
#   set.seed(23)
#   hclust_avg <- hclust(dist_mat, method = "ward.D")
#   cut_avg <- cutree(hclust_avg, k =4)
#   a=table(ind.coord$group,cut_avg)
#   print(i)
#   print(a)
#   if((a[2,2]+a[2,3]+a[2,4])<min){
#     min=a[2,2]+a[2,3]+a[2,4]
#     print(i)
#     print(a)
#   }
# }
# 
# # all9=all5[,c(1,f)]
# # dist_mat <- dist(all9[,2:40], method = "euclidean")
# # set.seed(23)
# # hclust_avg <- hclust(dist_mat, method = "ward.D")
# # cut_avg <- cutree(hclust_avg, k =3)
# # table(all9$group,cut_avg)
# # 

ind.coord$group1="Case"
ind.coord[ind.coord$group=="Control",]$group1="Control"
dist_mat <- dist(ind.coord[,1:74], method = "euclidean")
set.seed(23)
hclust_avg <- hclust(dist_mat, method = "ward.D")
cut_avg <- cutree(hclust_avg, k =3)
table(ind.coord$group,cut_avg)
ind.coord$group1=as.factor(ind.coord$group1)
ind.coord$cluster=cut_avg

table(all5$group,cut_avg)
all5$cluster=cut_avg

library(dendextend)
library(colorspace)
dend <- as.dendrogram(hclust_avg)
# order it the closest we can to the order of the observations:
dend <- rotate(dend, 1:nrow(ind.coord))
# Color the branches based on the clusters:
dend <- color_branches(dend, k=3) 
# Manually match the labels, as much as possible, to the real classification of the flowers:
labels_colors(dend) <-
  rainbow_hcl(2)[sort_levels_values(
    as.numeric(ind.coord[,(ncol(ind.coord))])[order.dendrogram(dend)]
  )]

# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend,hang_height=0.1)
# reduce the size of the labels:
# dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend <- set(dend, "labels_cex", 0.5)
# And plot:
par(mar = c(3,3,3,7))
plot(dend, 
     main = "Clustered
     (the labels give the true diagnosis)", 
     horiz =  FALSE,  nodePar = list(cex = .007))
s=table(ind.coord$group1,cut_avg)
colnames(s) =c("Cluster1","Cluster2","Cluster3")
s
library(ggplot2)
ggplot(ind.coord, aes(x=Dim.1, y =Dim.2, color = factor(cluster))) + geom_point()




all5$parti=substr(all5$hospid,1,nchar(all5$hospid)-4)
all5$parti
library(broom)
all6=all5[all5$group!="Control",]

all6$cluster=as.factor(all6$cluster)
colnames(all6)[1:10]
f=c()
for (i in c(5:799)){
  m=colnames(all6)[i]
  n=noquote(m)
  q=paste(n,"~","cluster")
  formula <- as.formula(q)
  res.aov <- aov(formula, data = all6)
  a=TukeyHSD(res.aov)
  b=sum(a$cluster[,4]<0.05/795)#res.aov$p.valuesum
  #print(m)
  #print(b)
  if(b==3){
    f=c(f,i)
  }
}

colnames(all6)[f]

f1=f
n=1
df=data.frame(matrix(ncol = 4, nrow = length(f1)))
x <- c("Parameters","Cluster1", "Cluster2", "Cluster3")
colnames(df) <- x
for (i in f){
  print(colnames(all6)[i])
  df[n,]$Parameters=colnames(all6)[i]
  mean1=round(mean(all6[all6$cluster==1,i]),2)
  sd1=round(sd(all6[all6$cluster==1,i]),2)
  cluster_1=paste(mean1,"±", sd1)
  print(paste("cluster_1",cluster_1))
  df[n,]$Cluster1=paste(cluster_1)
  mean2=round(mean(all6[all6$cluster==2,i]),2)
  sd2=round(sd(all6[all6$cluster==2,i]),2)
  cluster_2=paste(mean2,"±",sd2)
  print(paste("cluster_2",cluster_2))
  df[n,]$Cluster2=paste(cluster_2)
  mean3=round(mean(all6[all6$cluster==3,i]),2)
  sd3=round(sd(all6[all6$cluster==3,i]),2)
  cluster_3=paste(mean3,"±",sd3)
  print(paste("cluster_3",cluster_3))
  df[n,]$Cluster3=paste(cluster_3)
  n=n+1
}
# write.csv(df,"1.csv",na="",row.names = F)
# table(all6$cluster,all6$group)

###############validate 1 year##############
data=read.csv("Newchange_10-14_months.csv")
colnames(data)
subclinical=read.csv("subclinical.csv",header=F)
colnames(subclinical)=subclinical[1,]
subclinical=subclinical[-1,]
data[data$hospid %in% subclinical$hospid,]$group="Sub"
table(data$group)

####1.1-mean imputation#####
colnames(data)[1:10]
e=c()
for (i in 4:799){
  a=sum(is.na(data[,i]))
  if (a>0){
    e=c(e,i)
  }
}
#380
data1=data
for(i in 4:799){
  data1[is.na(data1[,i]), i]=round(mean(data1[,i], na.rm = TRUE),2)
}

e=c()
for (i in 4:799){
  a=sum(is.na(data1[,i]))
  if (a>0){
    e=c(e,i)
  }
}
#0
#data1$sex=as.numeric(data1$sex)

data2=data1
colnames(data2)

data2$summary_km_f
data3=data2[-which(data2$summary_km_f<0),]
table(data3$group)
summary(data3$summary_km_f)
table(data3$group)
summary(data3$summary_km_f)
data4=data3
summary(data4$follow_up)

#data4=data4[,-93]
data5=data4
library(ggpubr)
library(factoextra)
# Dimension reduction using PCA
colnames(data5)[1:10]
res.pca1 <- prcomp(data5[, 5:799],  scale = TRUE)
# Coordinates of individuals
ind.coord1 <- as.data.frame(get_pca_ind(res.pca1)$coord)
# Add Species groups from the original data sett
ind.coord1$group <- data5$group
# Data inspection
head(ind.coord1)

cluster_center = aggregate(scale(all5[,5:799]),list(cluster=cut_avg),mean)
d=cluster_center[,1:796]
d
# f=dist(ind.coord[554,1:78],ind.coord[ind.coord$cluster==3,1:78],method="euclidean")
# f=dist(ind.coord[554,1:78],ind.coord[ind.coord$cluster==1,1:78],method="euclidean")
# f=dist(ind.coord[554,1:78],ind.coord[ind.coord$cluster==2,1:78],method="euclidean")

# ind.coord$cluster=cut_avg
# ind.coord[554,]$cluster

# ind.coord2=ind.coord1
# ind.coord2$cluster=0

#b_pca =predict(res.pca, data5)
p=data5
p$cluster=0
p$group=data5$group
p[,5:799]=scale(p[,5:799])
for(i in 1:nrow(p)){
  print(i)
  m=dist(d[,-1],p[i,5:799])
  g=which.min(m)
  p[i,]$cluster=g
  print(g)
  
}
#m=dist(d[],scale(p[1,5:799]))

p$group1="Case"
p[p$group=="Control",]$group1="Control"
s=table(p$group1,p$cluster)
colnames(s) =c("Cluster1","Cluster2","Cluster3")
s

data5$cluster=p$cluster
table(data5$cluster,data5$group)
#summary_c_col_d5
#pachy_hapex_r_p0_dot_2
data6=data5[data5$group!="Control",]

f1=f
n=1
df=data.frame(matrix(ncol = 4, nrow = length(f1)))
x <- c("Parameters","Cluster1", "Cluster2", "Cluster3")
colnames(df) <- x
for (i in f){
  print(colnames(data6)[i])
  df[n,]$Parameters=colnames(data6)[i]
  mean1=round(mean(data6[data6$cluster==1,i]),2)
  sd1=round(sd(data6[data6$cluster==1,i]),2)
  cluster_1=paste(mean1,"±", sd1)
  print(paste("cluster_1",cluster_1))
  df[n,]$Cluster1=paste(cluster_1)
  mean2=round(mean(data6[data6$cluster==2,i]),2)
  sd2=round(sd(data6[data6$cluster==2,i]),2)
  cluster_2=paste(mean2,"±",sd2)
  print(paste("cluster_2",cluster_2))
  df[n,]$Cluster2=paste(cluster_2)
  mean3=round(mean(data6[data6$cluster==3,i]),2)
  sd3=round(sd(data6[data6$cluster==3,i]),2)
  cluster_3=paste(mean3,"±",sd3)
  print(paste("cluster_3",cluster_3))
  df[n,]$Cluster3=paste(cluster_3)
  n=n+1
}
# write.csv(df,"1.csv",na="",row.names = F)
# table(data6$cluster,data6$group)

###############validate 1.5 year##############
data=read.csv("Newchange_16-20_months.csv")
colnames(data)
subclinical=read.csv("subclinical.csv",header=F)
colnames(subclinical)=subclinical[1,]
subclinical=subclinical[-1,]
data[data$hospid %in% subclinical$hospid,]$group="Sub"
table(data$group)

####1.1-mean imputation#####
colnames(data)[1:10]
e=c()
for (i in 4:799){
  a=sum(is.na(data[,i]))
  if (a>0){
    e=c(e,i)
  }
}
#380
data1=data
for(i in 4:799){
  data1[is.na(data1[,i]), i]=round(mean(data1[,i], na.rm = TRUE),2)
}

e=c()
for (i in 4:799){
  a=sum(is.na(data1[,i]))
  if (a>0){
    e=c(e,i)
  }
}
#0
#data1$sex=as.numeric(data1$sex)

data2=data1
colnames(data2)

data2$summary_km_f
data3=data2[-which(data2$summary_km_f<0),]
table(data3$group)
summary(data3$summary_km_f)
table(data3$group)
summary(data3$summary_km_f)
data4=data3
summary(data4$follow_up)

#data4=data4[,-93]
data5=data4
library(ggpubr)
library(factoextra)
# Dimension reduction using PCA
colnames(data5)[1:10]
res.pca1 <- prcomp(data5[, 5:799],  scale = TRUE)
# Coordinates of individuals
ind.coord1 <- as.data.frame(get_pca_ind(res.pca1)$coord)
# Add Species groups from the original data sett
ind.coord1$group <- data5$group
# Data inspection
head(ind.coord1)

cluster_center = aggregate(scale(all5[,5:799]),list(cluster=cut_avg),mean)
d=cluster_center[,1:796]
d
# f=dist(ind.coord[554,1:78],ind.coord[ind.coord$cluster==3,1:78],method="euclidean")
# f=dist(ind.coord[554,1:78],ind.coord[ind.coord$cluster==1,1:78],method="euclidean")
# f=dist(ind.coord[554,1:78],ind.coord[ind.coord$cluster==2,1:78],method="euclidean")

# ind.coord$cluster=cut_avg
# ind.coord[554,]$cluster

# ind.coord2=ind.coord1
# ind.coord2$cluster=0

#b_pca =predict(res.pca, data5)
p=data5
p$cluster=0
p$group=data5$group
p[,5:799]=scale(p[,5:799])
for(i in 1:nrow(p)){
  print(i)
  m=dist(d[,-1],p[i,5:799])
  g=which.min(m)
  p[i,]$cluster=g
  print(g)
  
}
#m=dist(d[],scale(p[1,5:799]))

p$group1="Case"
p[p$group=="Control",]$group1="Control"
s=table(p$group1,p$cluster)
colnames(s) =c("Cluster1","Cluster2","Cluster3")
s

data5$cluster=p$cluster
table(data5$cluster,data5$group)
#summary_c_col_d5
#pachy_hapex_r_p0_dot_2
data6=data5[data5$group!="Control",]

f1=f
n=1
df=data.frame(matrix(ncol = 4, nrow = length(f1)))
x <- c("Parameters","Cluster1", "Cluster2", "Cluster3")
colnames(df) <- x
for (i in f){
  print(colnames(data6)[i])
  df[n,]$Parameters=colnames(data6)[i]
  mean1=round(mean(data6[data6$cluster==1,i]),2)
  sd1=round(sd(data6[data6$cluster==1,i]),2)
  cluster_1=paste(mean1,"±", sd1)
  print(paste("cluster_1",cluster_1))
  df[n,]$Cluster1=paste(cluster_1)
  mean2=round(mean(data6[data6$cluster==2,i]),2)
  sd2=round(sd(data6[data6$cluster==2,i]),2)
  cluster_2=paste(mean2,"±",sd2)
  print(paste("cluster_2",cluster_2))
  df[n,]$Cluster2=paste(cluster_2)
  mean3=round(mean(data6[data6$cluster==3,i]),2)
  sd3=round(sd(data6[data6$cluster==3,i]),2)
  cluster_3=paste(mean3,"±",sd3)
  print(paste("cluster_3",cluster_3))
  df[n,]$Cluster3=paste(cluster_3)
  n=n+1
}
# write.csv(df,"1.csv",na="",row.names = F)
# table(data6$cluster,data6$group)
