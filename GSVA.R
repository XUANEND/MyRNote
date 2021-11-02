library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)    #need R 3.6
library(pheatmap)

##########input data 
##基因集数据在pathway.txt文件中，也可以直接赋值sel_gmt
sel_gmt=read.csv("pathway.csv")
head(sel_gmt)
dim(sel_gmt)
sets=as.list(sel_gmt)
sets=lapply(sets, function(x) x[!is.na(x)])
#sets[1]

#########calculation
#基因数据在exprMatrix中，把自己的数据赋值代入即可
dat <- read.csv("expr.csv",header = T, row.names = 1)
exprMatrix=dat
head(exprMatrix)
dim(exprMatrix)
ssgsea_matrix<- gsva(as.matrix(exprMatrix), sets,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_matrix<- gsva(as.matrix(exprMatrix), sets,method='GSVA',kcdf='Gaussian',abs.ranking=TRUE)
gsva_matrix1<- t(scale(t(gsva_matrix)))

#head(gsva_matrix1)
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}
nor_gsva_matrix1 <- normalization(gsva_matrix1) 
dim(nor_gsva_matrix1)

score.gsva.metab=as.data.frame(t(nor_gsva_matrix1))
head(score.gsva.metab)


#k是预设的分组数，为了方便观察而已，可以自行修改或删除
p<-pheatmap(nor_gsva_matrix1,scale="row",show_colnames=F, show_rownames=F, cluster_cols=T, cluster_rows=T,cex=1, clustering_distance_rows="euclidean", cex=1,clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,cutree_col = 3)
p
cluster = cutree(p$tree_col,k=3) #from left to right 1 3 2
table(cluster)

#pdf("GSVA_heatmap.pdf",width=10,height=10)
#p
#dev.off()

#按热图分组提取样本
d=as.data.frame(cluster)
table(d$cluster)
d$sample=row.names(d)
a=as.data.frame(d[d[, "cluster"] == 1,])
a=dat[rownames(a),]
b=as.data.frame(d[d[, "cluster"] == 2,])
b=dat[rownames(b),]
c=as.data.frame(d[d[, "cluster"] == 3,])
c=dat[rownames(c),]