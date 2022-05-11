library(ggplot2)
library(Seurat)

## 人数据集Seurat对象创建及过滤

#读取表达矩阵
setwd("/Users/shangmengyue/Documents/RNA-seq/sc")
h_cardiomyogenesis<-read.table('GSE106118_UMI_count_merge.txt.gz',sep='\t',quote = "",fill = T,  
                               comment.char = "!",header = T)
rownames(h_cardiomyogenesis)<-h_cardiomyogenesis[,1]
h_cardiomyogenesis<-h_cardiomyogenesis[,-1]
head(h_cardiomyogenesis[1:4,1:4])

#替换细胞id字符
name<-colnames(h_cardiomyogenesis)
name<-gsub("_","-",name)
colnames(h_cardiomyogenesis)<-name
head(h_cardiomyogenesis[1:4,1:4])

h_cardiomyogenesis <- CreateSeuratObject(counts = h_cardiomyogenesis, project = "humfetal")
h_cardiomyogenesis

#添加注释信息
meta <- read.table('anno1.txt',sep='\t',quote = "",fill = T,
                   comment.char = "!",header = T)
setwd("/Users/shangmengyue/Documents/RNA-seq/sc/")
mappingrate<-read.table('mappingrate.txt',sep='\t',quote = "",fill = T,
                        comment.char = "!",header = T)
mappingrate<-mappingrate[,-c(3:12)]
mappingrate<-as.data.frame(mappingrate[,-1])


meta1 <- meta[colnames(meta)=="timepoint"]
meta2 <- meta[colnames(meta)=="stage"]
meta3 <- meta[colnames(meta)=="tissue"]
rownames(meta1)<-meta$sample.id
rownames(meta2)<-meta$sample.id
rownames(meta3)<-meta$sample.id

colnames(mappingrate)<-"Mapping.Rate"
rownames(mappingrate)<-meta$sample.id
h_cardiomyogenesis <- AddMetaData(object = h_cardiomyogenesis,
                                  metadata = meta1, col.name = "timepoint")
h_cardiomyogenesis <- AddMetaData(object = h_cardiomyogenesis,
                                  metadata = meta2, col.name = "stage")
h_cardiomyogenesis <- AddMetaData(object = h_cardiomyogenesis,
                                  metadata = meta3, col.name = "tissue")
h_cardiomyogenesis <- AddMetaData(object = h_cardiomyogenesis,
                                  metadata = mappingrate, col.name = "Mapping.Rate")

h_cardiomyogenesis[["percent.mt"]] <- PercentageFeatureSet(h_cardiomyogenesis, pattern = "MT")
VlnPlot(h_cardiomyogenesis, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)

#过滤
h_cardiomyogenesis <- subset(h_cardiomyogenesis, subset = nFeature_RNA > 1000 & nCount_RNA > 5000 & percent.mt < 5 )
h_cardiomyogenesis <- subset(h_cardiomyogenesis, subset = Mapping.Rate >0.4)
VlnPlot(h_cardiomyogenesis, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)

head(h_cardiomyogenesis@meta.data,4)


## 人数据集标准化、均一化
h_cardiomyogenesis <- NormalizeData(h_cardiomyogenesis, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(h_cardiomyogenesis)
h_cardiomyogenesis <- ScaleData(h_cardiomyogenesis, features = all.genes)

#保存
# saveRDS(h_cardiomyogenesis,"h_cardiomyogenesis.rds")

## 人数据集降维聚类
h_cardiomyogenesis <- readRDS("h_cardiomyogenesis.rds")
#以上已经normalize/scale完
h_cardiomyogenesis <- FindVariableFeatures(h_cardiomyogenesis, selection.method = "vst", nfeatures = 2500)
h_cardiomyogenesis <- RunPCA(h_cardiomyogenesis, features = VariableFeatures(object = h_cardiomyogenesis))

ElbowPlot(h_cardiomyogenesis)
h_cardiomyogenesis <- FindNeighbors(h_cardiomyogenesis, dims = 1:19,verbose=FALSE) ######这个dim的含义
h_cardiomyogenesis <- FindClusters(h_cardiomyogenesis,verbose=FALSE)  ###resolution降低，分群数变小
h_cardiomyogenesis <- RunUMAP(h_cardiomyogenesis, dims = 1:19)

umap = h_cardiomyogenesis@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% cbind(cluster = h_cardiomyogenesis@meta.data$seurat_clusters)%>% cbind(Developmental_stage =   h_cardiomyogenesis@meta.data$timepoint)

ggscatter(umap,x = "UMAP_1", y = "UMAP_2",
          color = "cluster",
          size = 0.7,
          font.label = c(12, "plain"),
)+ theme_bw()+theme(panel.grid.major = element_line(color="white"),
                    panel.grid.minor = element_line(color="white"),
                    panel.border = element_rect(size = 1.8),
                    legend.title = element_text(size=10),
                    legend.text = element_text(size=25))
ggscatter(umap,x = "UMAP_1", y = "UMAP_2",
          color = "Developmental_stage",
          size = 0.7,
          font.label = c(12, "plain"),
)+ theme_bw()+theme(panel.grid.major = element_line(color="white"),
                    panel.grid.minor = element_line(color="white"),
                    panel.border = element_rect(size = 1.8),
                    legend.title = element_text(size=10),
                    legend.text = element_text(size=25))
## 细胞类型差异基因
logFCfilter=0.5
adjPvalFilter=0.05
hum.markers <- FindAllMarkers(object = h_cardiomyogenesis,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)
sig.markers=hum.markers[(abs(as.numeric(as.vector(hum.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(hum.markers$p_val_adj))<adjPvalFilter),]
head(sig.markers[1:4,1:4])
# write.table (sig.markers,file ="2021hummarkers-umap.txt",sep='\t', quote =FALSE)
library(dplyr)
humtop10 <- hum.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
head(humtop10[1:4,1:4])
DoHeatmap(object = h_cardiomyogenesis, features = humtop10$gene) + NoLegend()

## 基于marker确定细胞类型
markerdf3<- as.data.frame(c("NKX2-5","MYL7","TNNT2"                   #CM
                            ,"EMCN","PECAM1","CDH5","COL1A1","DCN"    #Fib
                            ,"UPK3B","MSLN","WT1"                     #EP
                            ,"IL7R", "GZMA","CD68"                    #
                            ,"HBM", "HBQ1","LYZ"
                            ,"FN1","BGN"
))
colnames(markerdf3)<-"gene"
dot<-DotPlot(h_cardiomyogenesis, features = markerdf3$gene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none")
dot

#添加细胞类型注释
cell_type <- c("Fibroblast(3)","Blood cells","CM(4)","CM(7)","EP(1)","h_Fibroblast(4)","CM(1)","EC/Fibroblast(1)",
               "MONO/MAC","CM(8)","CM(2)","NKT","CM(6)","CM(3)","EC","CM(5)","Fibroblast(1)",
               "Fibroblast(2)","EP(2)","MONO","EC/Fibroblast(2)")
names(cell_type) <- levels(h_cardiomyogenesis)
h_cardiomyogenesis <- RenameIdents(h_cardiomyogenesis, cell_type)
DimPlot(h_cardiomyogenesis, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

cell_type_combine <- c("Fibroblast","Blood cells","CM","CM","EP","Fibroblast","CM","EC/Fibroblast","MONO/MAC","CM"
                       ,"CM","NKT","CM","CM","EC","CM","Fibroblast","Fibroblast",
                       "EP","MONO","EC/Fibroblast")
names(cell_type_combine) <- levels(h_cardiomyogenesis)
h_cardiomyogenesis <- RenameIdents(h_cardiomyogenesis, cell_type_combine)
DimPlot(h_cardiomyogenesis, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# saveRDS(h_cardiomyogenesis,"h_cardiomyogenesis.rds")

## 人细胞类型相关性验证
h_cardiomyogenesis <- readRDS("h_cardiomyogenesis.rds")
meta<-as.data.frame(h_cardiomyogenesis@meta.data)

library(plyr)
fun1<-function(data,xlab,fillc,pos,xname,yname){
  ggplot(data,aes(x=xlab,fill=fillc))+
    geom_bar(position = pos)+
    labs(x=xname,y=yname)+
    coord_flip()+
    theme_minimal()
}
p1<- fun1(meta,meta$stage,fillc = meta$`cell_type_combine`,pos='fill',xname='timepoint',yname='cell-type%')
p1

## 发育阶段的差异分析
Idents(h_cardiomyogenesis) <- "cell_type_combine"
h_cardiomyogenesis<- subset(h_cardiomyogenesis,idents=c("CM")) 
Idents(h_cardiomyogenesis) <- "stage"
library(tidyverse)
logFCfilter=0.5
adjPvalFilter=0.05
hum.markers <- FindAllMarkers(object = h_cardiomyogenesis,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)

all.markers=hum.markers[(abs(as.numeric(as.vector(hum.markers$avg_logFC)))>logFCfilter & as.numeric(as.vector(hum.markers$p_val_adj))<adjPvalFilter),]

humtop10 <- hum.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# write.table (humtop10,file ="hum_stagetop10.txt",sep='\t', quote =FALSE) 
# write.table (hum.markers,file ="hum_stageDEG.txt",sep='\t', quote =FALSE) 


DoHeatmap(object = hpbmccm, features = humtop10$gene) + NoLegend()+scale_fill_gradientn(colors = c("whitesmoke", "lavenderblush3", "violetred4"))


# saveRDS(h_cardiomyogenesis,"h_cardiomyogenesis.rds")

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(Seurat)
library(corrplot)
library(stringr)
library(clusterProfiler)

## 发育阶段的细胞构成
h_cardiomyogenesis <- readRDS("h_cardiomyogenesis.rds")
meta<-as.data.frame(h_cardiomyogenesis@meta.data)

library(plyr)
fun1<-function(data,xlab,fillc,pos,xname,yname){
  ggplot(data,aes(x=xlab,fill=fillc))+
    geom_bar(position = pos)+
    labs(x=xname,y=yname)+
    coord_flip()+
    theme_minimal()
}
p1<- fun1(meta,meta$stage,fillc = meta$`cell_type_combine`,pos='fill',xname='timepoint',yname='cell-type%')
p1

library(monocle)
library(Seurat)

## 人心肌细胞轨迹构建
hum_cm<-readRDS("人cm分群.rds")
expr_matrix <- as(as.matrix(hum_cm@assays$RNA@counts), 'sparseMatrix')

##提取表型信息到p_data(phenotype_data)里面 
p_data <- hum_cm@meta.data 
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(hum_cm),row.names = row.names(hum_cm))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)
#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
#估计size factor和离散度
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

diff <- differentialGeneTest(cds,fullModelFormulaStr="~timepoint",cores=1)
deg <- subset(diff, qval < 0.01)
deg <- deg[order(deg$qval,decreasing=F),]
ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)
#saveRDS(cds,"人cm时间点拟时.rds")
cds <- orderCells(cds, root_state = 4)

Time_diff <- differentialGeneTest(cds[ordergene,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=3, show_rownames=T, return_heatmap=T)
ggsave("new人cm时间点拟时热图.pdf", p, width = 5, height = 10)

# Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
# write.csv(Time_diff, "Time_diff_all.csv", row.names = F)
write.table(Time_diff,file="new人cm时间点拟时genenew.txt",sep="\t",row.names=T,quote=F)
clusters <- cutree(p$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.table(clustering,file="CM_Time_clustering_all_new.txt",sep="\t",row.names=T,quote=F)
print("CM DONE")

#### 人成纤维细胞及内皮细胞/成纤维细胞与心肌细胞轨迹构建过程一致

## 鼠数据集Seurat对象创建及过滤

##小鼠数据集2采样
library(reticulate)
library(rsvd)

geosketch <- import('geosketch')

# Generate some random data.
# X <- replicate(20, rnorm(1000))

# Get top PCs from randomized SVD.
s <- rsvd(X, k=19)
X.pcs <- s$u %*% diag(s$d)

# Sketch 10% of data.
sketch.size <- as.integer(3000)
sketch.indices <- geosketch$gs(X.pcs, sketch.size, one_indexed = TRUE)
print(sketch.indices)
# X_sketch <- X.pcs[A$V1,]

A<- t(as.data.frame(sketch.indices))

X1<-X[A,]
X1<-t(X1)
setwd("/Users/shangmengyue/Documents/RNA-seq/sc/小鼠单独/")
# write.table(X1,file="mus2采样.txt",sep="\t",row.names=T,quote=F)
samplename<-as.data.frame(colnames(X1))
colnames(samplename)<-"cell"
sampleanno<- left_join(samplename,annomerge3,by="cell",all.x = T)
setwd("/Users/shangmengyue/Documents/RNA-seq/sc/小鼠单独/")
# write.table(sampleanno,file="mus2采样anno.txt",sep="\t",row.names=T,quote=F)

## 小鼠数据集整合
setwd("/Users/shangmengyue/Documents/RNA-seq/sc/202012/0512/")
mus1_2<-read.table('mus1m_NAME.txt',sep='\t',quote = "",fill = T,
                   comment.char = "!",header = T)
setwd("/Users/shangmengyue/Documents/RNA-seq/sc/小鼠单独/")
mus2_2<-read.table('mus2采样.txt',sep='\t',quote = "",fill = T,
                   comment.char = "!",header = T)
mus1_2[1:4,1:4]
mus2_2[1:4,1:4]
setwd("/Users/shangmengyue/Documents/RNA-seq/sc/")
mus1_anno<-read.table('mus1-meta202109.txt',sep='\t',quote = "",fill = T,
                      comment.char = "!",header = T)
setwd("/Users/shangmengyue/Documents/RNA-seq/sc/小鼠单独/")
sampleanno<-read.table('mus2采样anno.txt',sep='\t',quote = "",fill = T,
                       comment.char = "!",header = T)
sampleanno<-sampleanno[,-1]
mus <- CreateSeuratObject(counts = cbind(mus1_2, mus2_2), project = "mus", min.cells = 5,meta.data = rbind(mus1_anno,sampleanno)) 
mus[["percent.mt"]] <- PercentageFeatureSet(mus, pattern = "mt")

VlnPlot(mus, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# mus[["percent.mt"]] <- PercentageFeatureSet(mus, pattern = "mt")
mus <- subset(mus, subset = percent.mt < 10 )
VlnPlot(mus, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

mus <- NormalizeData(mus, normalization.method = "LogNormalize", scale.factor = 10000)
mus<- ScaleData(mus,verbose = FALSE)
mus<- FindVariableFeatures(mus, selection.method = "vst",nfeatures = 2000)
mus <- RunPCA( mus, features = VariableFeatures(object =  mus))

setwd("/Users/shangmengyue/Documents/RNA-seq/sc/小鼠单独/")
options(repr.plot.height = 5, repr.plot.width = 12)
pdf(file = "harmony前.pdf",height = 3,width = 8)

p1 <- DimPlot(object = mus, reduction = "pca", pt.size = .1, group.by = "tech")
p2 <- VlnPlot(object = mus, features = "PC_1", group.by = "tech",  pt.size = .1)
plot_grid(p1,p2)
dev.off()

library(harmony)
# options(repr.plot.height = 2.5, repr.plot.width = 6)
mus <- mus %>% 
  RunHarmony("tech", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(mus, 'harmony')
harmony_embeddings[1:5, 1:5]

## 小鼠数据集标准化、分群
# pdf(file = "harmony前后.pdf",height = 6,width = 8)
# plot_grid(p1,p2,p3,p4)
# dev.off()
ElbowPlot(mus)
mus <- FindNeighbors(mus, dims = 1:19 , verbose=FALSE) ######这个dim的含义
mus <- FindClusters(mus,verbose=FALSE)  ###resolution降低，分群数变小
mus <- RunUMAP(mus, dims = 1:19)
mus <- RunTSNE(object = mus, dims = 1:19, do.fast = TRUE)

mus <- FindNeighbors(mus,reduction = "harmony",dims = 1:19 , verbose=FALSE) ######这个dim的含义
mus <- FindClusters(mus,verbose=FALSE)  ###resolution降低，分群数变小
mus <- RunUMAP(mus, reduction = "harmony",dims = 1:19)
mus <- RunTSNE(object = mus, reduction = "harmony",dims = 1:19, do.fast = TRUE)
# pdf(file = "采样harmony后聚类.pdf",height = 3,width = 8)
# p7<-DimPlot(mus, reduction = "umap",pt.size = .1,group.by = "tech")
# p8<-DimPlot(mus, reduction = "tsne",pt.size = .1,group.by = "tech")
# plot_grid(p7,p8)
# dev.off()
logFCfilter=0.5
adjPvalFilter=0.05
markers <- FindAllMarkers(object = mus,
                          only.pos = FALSE,
                          min.pct = 0.25,
                          logfc.threshold = logFCfilter)
sig.markers=markers[(abs(as.numeric(as.vector(markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(markers$p_val_adj))<adjPvalFilter),]
head(sig.markers[1:4,1:4])
# write.table (sig.markers,file ="2021hummarkers-umap.txt",sep='\t', quote =FALSE)
library(dplyr)
mustop10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
head(mustop10[1:4,1:4])
tiff("鼠采样整合marker.tiff",units = "in",width = 10,height = 15,res = 600,compression = "lzw")
DoHeatmap(object = mus, features = mustop10$gene) + NoLegend()
dev.off()
#定义细胞类型
markerdf3<- as.data.frame(c("Tnnt2","Tnnc1","Actn2","Myh6","Tnni3",  #CM
                            "Pecam1","Cd93","Tek", "Cdh5",           #EC
                            "Col3a1","Col1a2","Fn1","Col1a1","Vim",  #Fib
                            "Raldh2","Tcf21","Wt1","Tbx18",          #EP
                            "Gzma","Klrb1c",                         #NKT
                            "Il13","II6","Mcpt8",                    #MAST
                            "Adgre1","Cd14","Ccr2",                  #MONO
                            "Iba1","Fcgr1","C5ar1", "Cd163",         #MAC
                            "Hbb-bt", "Hba-a2","Hbb-y"
                            
))
colnames(markerdf3)<-"gene"
dot<-DotPlot(mus, features = markerdf3$gene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  theme_bw()+
  theme(
    axis.text.y = element_text(size=20),
    axis.text.x.bottom = element_text(angle = 90,size=20,hjust = 1,vjust = 1),
    panel.border = element_rect(size = 1.8),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none",)

pdf(file = "bubble采样.pdf",height = 6,width = 8)
dot
dev.off()
saveRDS(mus,"采样整合小鼠-注释.rds")

##鼠CM发育轨迹构建
mus_cm<-readRDS("鼠cm分群.rds")
expr_matrix <- as(as.matrix(mus_cm@assays$RNA@counts), 'sparseMatrix')

##提取表型信息到p_data(phenotype_data)里面 
p_data <- mus_cm@meta.data 
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(mus_cm),row.names = row.names(mus_cm))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)
#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
#估计size factor和离散度
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# test <- FindVariableFeatures(mus_cm, selection.method = "vst", nfeatures = 2000)
# 
# express_genes <- VariableFeatures(test)
# cds <- setOrderingFilter(cds, express_genes)
diff <- differentialGeneTest(cds,fullModelFormulaStr="~Developmental_stage",cores=1)
deg <- subset(diff, qval < 0.01)
deg <- deg[order(deg$qval,decreasing=F),]
ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 3)

#saveRDS(cds,"鼠cm时间点拟时.rds")
Time_diff <- differentialGeneTest(cds[ordergene,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
# Time_diff <- differentialGeneTest(cds[ordergene,], cores = 1, 
#                                   fullModelFormulaStr = "~sm.ns(State)")

Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=3, show_rownames=T, return_heatmap=T)
ggsave("new鼠cm时间点拟时热图newtest.pdf", p, width = 5, height = 10)

# Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
# write.csv(Time_diff, "Time_diff_all.csv", row.names = F)
write.table(Time_diff,file="new鼠cm时间点拟时genenew.txt",sep="\t",row.names=T,quote=F)
clusters <- cutree(p$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.table(clustering,file="鼠cm_Time_clustering_allnew.txt",sep="\t",row.names=T,quote=F)
print("DONE")
#### 鼠成纤维细胞及内皮细胞/成纤维细胞与心肌细胞轨迹构建过程一致

##人鼠CM时间点映射
hpbmc<-readRDS("人umap_celltype.rds")
Idents(hpbmc) <- "cell_type_combine"
hum_cm<- subset(hpbmc,idents=c("CM")) 
hum_cm <- NormalizeData(hum_cm, verbose = FALSE)
hum_cm <- FindVariableFeatures(hum_cm, selection.method = "vst", 
                               nfeatures = 2000, verbose = FALSE)
hum_cm <- ScaleData(hum_cm, verbose = FALSE)
hum_cm <- RunPCA(hum_cm, features = VariableFeatures(object = hum_cm))

ElbowPlot(hum_cm)
hum_cm <- FindNeighbors(hum_cm, dims = 1:19,verbose=FALSE) ######这个dim的含义
hum_cm <- FindClusters(hum_cm,verbose=FALSE)  ###resolution降低，分群数变小
hum_cm <- RunUMAP(hum_cm, dims = 1:19)
hum_cm <- RunTSNE(object = hum_cm, dims = 1:19, do.fast = TRUE)
DimPlot(hum_cm, reduction = "umap", group.by = "timepoint", pt.size = .5)
DimPlot(hum_cm, reduction = "tsne", group.by = "timepoint", pt.size = .5)




#提取在75%CM(1512*0.80=1209)中表达的gene
# 提取计数
counts <- GetAssayData(object = hum_cm, slot = "counts")

# 根据在每个细胞的计数是否大于0为每个基因输出一个逻辑向量
nonzero <- counts > 0

# 将所有TRUE值相加，如果每个基因的TRUE值超过10个，则返回TRUE。
keep_genes <- Matrix::rowSums(nonzero) >= 1134

# 仅保留那些在10个以上细胞中表达的基因
filtered_counts <- counts[keep_genes, ]
filtered_counts2 <- as.data.frame(filtered_counts)
hum_gene_80percent_cm<-as.data.frame(rownames(filtered_counts2))
write.table(hum_gene_80percent_cm,file="hum_gene_80percent_cm.txt",sep="\t",row.names=T,quote=F)


# 重新赋值给经过过滤的Seurat对象
hum_cm_75p <- CreateSeuratObject(filtered_counts, meta.data = hum_cm@meta.data)




hum_cm_75p <- NormalizeData(hum_cm_75p,normalization.method = "LogNormalize",scale.factor = 100)
all.genes <-rownames(hum_cm_75p)
hum_cm_75p <- ScaleData(hum_cm_75p,features = all.genes)
hum_cm_75p <- FindVariableFeatures(hum_cm_75p, selection.method = "vst", 
                                   nfeatures = 2000, verbose = FALSE)

hum_cm_75p <- RunPCA(hum_cm_75p, features = VariableFeatures(object = hum_cm_75p))
# DimPlot(hum_cm_75p, reduction = "pca", pt.size = 0.5,group.by = "stage") 
# DimPlot(hum_cm_75p, reduction = "pca", pt.size = 0.5,group.by = "timepoint") 


#pc1坐标
cells<-colnames(x=human_cm)
dims <- c(1,2)
reduction <- "pca"
data <- Embeddings(object = human_cm[[reduction]])[cells, dims]
data <- as.data.frame(x = data)
dims <- paste0(Key(object = human_cm[[reduction]]), dims)

plot <- ggplot(data = data) +
  geom_point(
    mapping = aes_string(
      x = dims[1],
      y = dims[2]
    ),
    size = 1.2
  )
plot

library(irlba)
npcs = 50

head(data)
data$cell<-rownames(data)
allmeta<-human_cm@meta.data
allmeta$cell<-rownames(allmeta)
data<- left_join(data,allmeta,by="cell",all.x = T)
rownames(data)<-data$cell

library(ggplot2)
p<-ggplot(data, aes(x = PC_1))+ geom_density(aes(color = Developmental_stage))
# dev.off()

##提取loading value of PC1
cell_embeddings<-as.data.frame(human_cm@reductions$pca@cell.embeddings)
feature_loadings<-as.data.frame(human_cm@reductions$pca@feature.loadings)
feature_loadings_2pc <- feature_loadings[,1:2]


#####鼠1和鼠2均值映射
Idents(mus2_for_pca) <- "Developmental_stage"
e13.5<- subset(mus2_for_pca,idents=c("e13.5")) 
e13.5 <- NormalizeData(e13.5,normalization.method = "LogNormalize",scale.factor = 100)
all.genes <-rownames(e13.5)
e13.5 <- ScaleData(e13.5,features = all.genes)

mus13.5_scale_2000<-as.data.frame(e13.5@assays$RNA@scale.data)
e13_5mean<-as.data.frame(rowMeans(mus13.5_scale_2000))
material1<-as.matrix(e13_5mean)
e13_5mean_data<-crossprod(material1,material2)

e14.5<- subset(mus2_for_pca,idents=c("e14.5")) 
e14.5 <- NormalizeData(e14.5,normalization.method = "LogNormalize",scale.factor = 100)
all.genes <-rownames(e14.5)
e14.5 <- ScaleData(e14.5,features = all.genes)
mus14.5_scale_2000<-as.data.frame(e14.5@assays$RNA@scale.data)
e14_5mean<-as.data.frame(rowMeans(mus14.5_scale_2000))
material1<-as.matrix(e14_5mean)
e14_5mean_data<-crossprod(material1,material2)

mus13.5_scale_2000<-counts_2000[,1:805]
e13_5mean<-as.data.frame(rowMeans(mus13.5_scale_2000))
material1<-as.matrix(e13_5mean)
e13_5mean_data<-crossprod(material1,material2)

mus14.5_scale_2000<-counts_2000[,806:1458]
e14_5mean<-as.data.frame(rowMeans(mus14.5_scale_2000))
material1<-as.matrix(e14_5mean)
e14_5mean_data<-crossprod(material1,material2)
Idents(mus1_for_pca) <- "Developmental_stage"
e8.5<- subset(mus1_for_pca,idents=c("e8.5")) 

mus8.5_scale_2000<-as.data.frame(e8.5@assays$RNA@scale.data)
e8_5mean<-as.data.frame(rowMeans(mus8.5_scale_2000))
material1<-as.matrix(e8_5mean)
e8_5mean_data<-crossprod(material1,material2)

e9.5<- subset(mus1_for_pca,idents=c("e9.5")) 
mus9.5_scale_2000<-as.data.frame(e9.5@assays$RNA@scale.data)
e9_5mean<-as.data.frame(rowMeans(mus9.5_scale_2000))
material1<-as.matrix(e9_5mean)
e9_5mean_data<-crossprod(material1,material2)

e10.5<- subset(mus1_for_pca,idents=c("e10.5")) 
mus10.5_scale_2000<-as.data.frame(e10.5@assays$RNA@scale.data)
e10_5mean<-as.data.frame(rowMeans(mus10.5_scale_2000))
material1<-as.matrix(e10_5mean)
e10_5mean_data<-crossprod(material1,material2)


