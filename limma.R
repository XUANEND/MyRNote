rm(list=ls())
library(GEOquery)
library(limma)
GSE60291 <- getGEO('GSE60291', destdir=".",getGPL = F)

#下面是表达矩阵
exprSet=exprs(GSE60291[[1]])
library("annotate")
#GSE60291[[1]]
## 下面是分组信息
pdata=pData(GSE60291[[1]])
treatment=factor(unlist(lapply(pdata$title,function(x) strsplit(as.character(x),"-")[[1]][1])))#lapply()函数格式不太懂
#treatment=relevel(treatment,'control')
## 下面做基因注释
platformDB='hgu133plus2.db'
library(platformDB, character.only=TRUE)
probeset <- featureNames(GSE60291[[1]])
#EGID <- as.numeric(lookUp(probeset, platformDB, "ENTREZID"))
SYMBOL <-  lookUp(probeset, platformDB, "SYMBOL")
## 下面对每个基因挑选最大表达量探针
a=cbind(SYMBOL,exprSet)
## remove the duplicated probeset
rmDupID <-function(a=matrix(c(1,1:5,2,2:6,2,3:7),ncol=6)){ #这里a赋值是不必要的，大概是为了程序的可读性加的
  exprSet=a[,-1]
  rowMeans=apply(exprSet,1,function(x) mean(as.numeric(x),na.rm=T))
  a=a[order(rowMeans,decreasing=T),]
  exprSet=a[!duplicated(a[,1]),]
  exprSet=exprSet[!is.na(exprSet[,1]),]
  rownames(exprSet)=exprSet[,1]
  exprSet=exprSet[,-1]
  return(exprSet)
}
exprSet=rmDupID(a)
rn=rownames(exprSet)
exprSet=apply(exprSet,2,as.numeric)
rownames(exprSet)=rn
exprSet[1:4,1:4]
#exprSet=log(exprSet) ## based on e
boxplot(exprSet,las=2)
## 下面用limma包来进行芯片数据差异分析
design=model.matrix(~ treatment)
fit=lmFit(exprSet,design)
fit=eBayes(fit)
#vennDiagram(decideTests(fit))
DEG=topTable(fit,coef=2,n=Inf,adjust='BH')
dim(DEG[abs(DEG[,1])>1.2 & DEG[,5]<0.05,])  ## 806 genes
write.csv(DEG,"ET1-normal.DEG.csv")
