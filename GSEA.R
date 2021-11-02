
# R包的安装
if(!requireNamespace("BiocManager",quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

# R包的导入
library(clusterProfiler) # GSEA分析
library(org.Hs.eg.db)    # 基因名转换
library(enrichplot)  


# 读入数据
df <- read.table("input.txt",sep="\t",header=T,check.names=F)
df = DEG
df$gene = rownames(DEG)
df = df[c("gene","logFC")]
head(df)


# ENTREZ ID注释
df.id <- bitr(df$gene, # 转换的列是df数据框中的genem基因名列
              fromType = "SYMBOL", # 需要转换ID类型
              toType = "ENTREZID", # 转换成的ID类型
              OrgDb = "org.Hs.eg.db" ) # 对应的物种
head(df.id)


names(df) <- c("SYMBOL","logFC") #  重新命名df的列名
easy.df<-merge(df,df.id,by="SYMBOL",all=F) # 以SMBOL列为依据，合并df和df.id
head(easy.df)


sortdf<-easy.df[order(easy.df$logFC, decreasing = T),]
head(sortdf)


gene.expr = sortdf$logFC #把logFC按照从大到小提取出来
head(gene.expr)
names(gene.expr) <- sortdf$ENTREZID#给上面提取的logFC对应上ENTREZID
head(gene.expr)


# 读取预设基因集
hallmarks <- read.gmt("gene_sets.gmt") # 预先准备好的预设基因集gmt文件格式
# annotation
hallmarks.id <- bitr(hallmarks$gene, # 转换的列是df数据框中的genem基因名列
                     fromType = "SYMBOL", # 需要转换ID类型
                     toType = "ENTREZID", # 转换成的ID类型
                     OrgDb = "org.Hs.eg.db" ) # 对应的物种
head(hallmarks.id)


names(hallmarks) <- c("term","SYMBOL") #  重新命名df的列名
easy.hallmarks<-merge(hallmarks,hallmarks.id,by="SYMBOL",all=F) # 以SMBOL列为依据，合并df和df.id
head(easy.hallmarks)


easy.hallmarks <- easy.hallmarks[,-1] # 删除SYMBOL列
head(easy.hallmarks)


# GSEA分析
y <- GSEA(gene.expr,TERM2GENE = easy.hallmarks)
head(y)


for(i in 1:99){
  pdf(paste("img/GSEA_",as.character(i),".pdf",sep = ""))
  p1 <- gseaplot2(y,#数据
            geneSetID = i,#画那一列的信号通路
            title = "GSEA Analysis",#标题
            base_size = 15,#字体大小
            color = "#483D8B",#线条的颜色
            pvalue_table = TRUE,#加不加p值
            rel_heights = c(1.5, 0.5, 1),
            subplots = 1:3,
            ES_geom="line")#是用线，还是用d点
  print(p1)
  dev.off()
}

gseaplot2(y,#数据
          geneSetID = 3,#画那一列的信号通路
          title = "GSEA Analysis",#标题
          base_size = 15,#字体大小
          color = "#483D8B",#线条的颜色
          pvalue_table = TRUE,#加不加p值
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          ES_geom="line")#是用线，还是用d点

gseaplot2(GO,#数据
          geneSetID = 1,#画那一列的信号通路
          title = "GSEA Analysis",#标题
          base_size = 15,#字体大小
          color = "#483D8B",#线条的颜色
          pvalue_table = TRUE,#加不加p值
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          ES_geom="line")#是用线，还是用d点
