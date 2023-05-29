#!/usr/bin/evn R
library(limma)
library(NMF)
library(ggplot2)
library(reshape)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)

lnames <- load("forLimma.RData")
#Expvalues=read.csv("Sarah.csv")
groups=read.csv("groups.csv")
samplenames=c(groups$RID)
miRNAid=read.csv("miRNAlist.csv")
miRNAid=miRNAid[,1]
miRNAid <- gsub("-","_",miRNAid)
#Outliers
# 2267:PBS XXM
# 3    3   PBS XXM      2267   XMP
# 2033:PBS XXF
# 4    4   PBS XXF      2033   XFP
dCtLimmaMatrix <- subset(dCtLimmaMatrix,select=-c(X3,X4))
groups <- groups[-c(3,4),]
samplenames=c(groups$RID)
########
Expvalues <- dCtLimmaMatrix[miRNAid,]
colnames(Expvalues)=samplenames

groups$Treatment <- factor(groups$Treatment, levels = c("PBS XXM", "PBS XXF", "PBS XYM", "PBS XYF", "HDM XXM", "HDM XXF", "HDM XYM","HDM XYF"))
groups$group <- factor(groups$group, levels = c("XMP", "XFP", "YMP", "YFP", "XMH", "XFH", "YMH", "YFH"))

pdf("deltaCt_limma_plots.pdf")
#######################
sampleDists <- dist( t( Expvalues ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(groups$Sample.ID,":",groups$Treatment,sep="")
colnames(sampleDistMatrix) <- paste(groups$Sample.ID,":",groups$Treatment,sep="")
ph <- pheatmap(sampleDistMatrix, color = greenred(256), border_color = NA, main="Sample vs Sample", fontsize_row=6, fontsize_col=6, show_rownames=TRUE, cluster_cols=TRUE, scale="none")
###
pca <- prcomp(t(na.omit(Expvalues)))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
#  write.table(as.data.frame(pca$rotation), file="gene_loadings.txt", quote=FALSE, sep="\t", col.names=NA)
#  write.table(percentVar, file="explained_variance.txt", quote=FALSE, sep="\t", col.names=NA)
data <- cbind(pca$x,data.frame(Treatment=groups$Treatment, name = groups$Sample.ID ))
attr(data, "percentVar") <- percentVar

percentVar <- round(100 * attr(data, "percentVar"))
print(ggplot(data, aes(PC1, PC2, color=Treatment, label=name)) + geom_text_repel(max.overlaps = 14, box.padding = 0.5, min.segment.length = 0, size=2) +
        geom_point(size=3) + ggtitle("PC1 vs PC2") + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) )

mat <- Expvalues
colnames(mat) <- paste(groups$Sample.ID,":",groups$Treatment,sep="")
df <- data.frame(Treatment=groups$Treatment)
rownames(df) <- colnames(mat)
mycols <- brewer.pal(8, "Dark2")[1:length(unique( df$Treatment ))]
ann_colors <- list(Treatment=c(mycols))
names(ann_colors$Treatment) <- unique(df$Treatment)
title <- "deltaCt Values Heatmap"
ph <- pheatmap(mat, color = colorRampPalette(c("navy", "white", "red"))(256), annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=6, fontsize_col=6, show_rownames=TRUE, cluster_cols=FALSE, scale="row")
ph <- pheatmap(mat, color = colorRampPalette(c("navy", "white", "red"))(256), annotation_col = df, annotation_colors = ann_colors, border_color = NA, main=title, fontsize_row=6, fontsize_col=6, show_rownames=TRUE, cluster_cols=TRUE, scale="row")
#######################

aheatmap(Expvalues, color = colorRampPalette(c("blue","white","red"))(20), scale="row", annColors = "Set1", Colv = NULL, Rowv= NULL, annLegend = TRUE, annCol = groups$group, labCol = groups$Treatment, labRow = miRNAid, border_color = "black")
aheatmap(Expvalues, color = colorRampPalette(c("blue","white","red"))(20), scale="row", annColors = "Set1", Colv = NA, Rowv= NA, annLegend = TRUE, annCol = groups$group, labCol = groups$Treatment, labRow = miRNAid, border_color = "black")

### log expression heatmap (green vs red)
### aheatmap(logEV, color = colorRampPalette(c("green","white","red"))(20), scale="row", annColors = "Set1", Colv = NULL, Rowv= NULL, annLegend = TRUE, annCol = groups$group, labCol = groups$Treatment,labRow = miRNAid, border_color = "black")
### aheatmap(logEV, color = colorRampPalette(c("green","white","red"))(20), scale="row", annColors = "Set1", Colv = NA, Rowv= NA, annLegend = TRUE, annCol = groups$group, labCol = groups$Treatment, labRow = miRNAid, border_color = "black")

################ Limma ################
design=model.matrix(~0 + groups$group)
#colnames(design)=c("XFH","XFP", "XMH","XMP", "YFH", "YFP", "YMH", "YMP")
colnames(design)=gsub("groups\\$group","",colnames(design))
rownames(design)=groups$RID

fit=lmFit(Expvalues,design)
contrast.matrix=makeContrasts(YM_HDMvPBS=YMH-YMP, XM_HDMvPBS=XMH-XMP, XF_HDMvPBS=XFH-XFP, YF_HDMvPBS=YFH-YFP, levels=design)
fit=contrasts.fit(fit,contrast.matrix)
fit=eBayes(fit, trend = TRUE)


#RESULTS - top tables of miRNA comparisons
###########################################
#1) HDM vs PBS XYM
YM_HDMvPBS=topTable(fit, coef = "YM_HDMvPBS", sort.by="p", adjust.method = "fdr", genelist=miRNAid, n=Inf)
write.table(YM_HDMvPBS,file="YM_HDMvPBS.deltaCtLimma.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)

volcanoplot(fit, coef = "YM_HDMvPBS", style = "p-value")
#boxplot top 3 miRNAs
boxplot(Expvalues[rownames(YM_HDMvPBS[1,]),] ~ groups$group, main = rownames(YM_HDMvPBS[1,]), ylab = "miRNA expression", xlab=NULL, names=levels(groups$Treatment),las=2, col = ifelse(grepl("F",levels(groups$Treatment)),"pink","lightblue"))
boxplot(Expvalues[rownames(YM_HDMvPBS[2,]),] ~ groups$group, main = rownames(YM_HDMvPBS[2,]), ylab = "miRNA expression", xlab=NULL, names=levels(groups$Treatment),las=2, col = ifelse(grepl("F",levels(groups$Treatment)),"pink","lightblue"))
boxplot(Expvalues[rownames(YM_HDMvPBS[3,]),] ~ groups$group, main = rownames(YM_HDMvPBS[3,]), ylab = "miRNA expression", xlab=NULL, names=levels(groups$Treatment),las=2, col = ifelse(grepl("F",levels(groups$Treatment)),"pink","lightblue"))
###########################################

#2) HDM vs PBS XXM
XM_HDMvPBS=topTable(fit, coef = "XM_HDMvPBS", sort.by="p", adjust.method = "fdr", genelist=miRNAid, n=Inf)
write.table(XM_HDMvPBS,file="XM_HDMvPBS.deltaCtLimma.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)

volcanoplot(fit, coef = "XM_HDMvPBS", style = "p-value")
#boxplot top 3 miRNAs
boxplot(Expvalues[rownames(XM_HDMvPBS[1,]),] ~ groups$group, main = rownames(XM_HDMvPBS[1,]), ylab = "miRNA expression", xlab=NULL, names=levels(groups$Treatment),las=2, col = ifelse(grepl("F",levels(groups$Treatment)),"pink","lightblue"))
boxplot(Expvalues[rownames(XM_HDMvPBS[2,]),] ~ groups$group, main = rownames(XM_HDMvPBS[2,]), ylab = "miRNA expression", xlab=NULL, names=levels(groups$Treatment),las=2, col = ifelse(grepl("F",levels(groups$Treatment)),"pink","lightblue"))
boxplot(Expvalues[rownames(XM_HDMvPBS[3,]),] ~ groups$group, main = rownames(XM_HDMvPBS[3,]), ylab = "miRNA expression", xlab=NULL, names=levels(groups$Treatment),las=2, col = ifelse(grepl("F",levels(groups$Treatment)),"pink","lightblue"))
###########################################

#3) HDM vs PBS XXF
XF_HDMvPBS=topTable(fit, coef = "XF_HDMvPBS", sort.by="p", adjust.method = "BH", genelist=miRNAid, n=Inf)
write.table(XF_HDMvPBS,file="XF_HDMvPBS.deltaCtLimma.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)

volcanoplot(fit, coef = "XF_HDMvPBS", style = "p-value")
#boxplot top 3 miRNAs
boxplot(Expvalues[rownames(XF_HDMvPBS[1,]),] ~ groups$group, main = rownames(XF_HDMvPBS[1,]), ylab = "miRNA expression", xlab=NULL, names=levels(groups$Treatment),las=2, col = ifelse(grepl("F",levels(groups$Treatment)),"pink","lightblue"))
boxplot(Expvalues[rownames(XF_HDMvPBS[2,]),] ~ groups$group, main = rownames(XF_HDMvPBS[2,]), ylab = "miRNA expression", xlab=NULL, names=levels(groups$Treatment),las=2, col = ifelse(grepl("F",levels(groups$Treatment)),"pink","lightblue"))
boxplot(Expvalues[rownames(XF_HDMvPBS[3,]),] ~ groups$group, main = rownames(XF_HDMvPBS[3,]), ylab = "miRNA expression", xlab=NULL, names=levels(groups$Treatment),las=2, col = ifelse(grepl("F",levels(groups$Treatment)),"pink","lightblue"))
###########################################

#4) HDM vs PBS XYF
YF_HDMvPBS=topTable(fit, coef = "YF_HDMvPBS", sort.by="p", adjust.method = "BH", genelist=miRNAid, p.value=1, n=Inf)
write.table(YF_HDMvPBS,file="YF_HDMvPBS.deltaCtLimma.tsv", quote=FALSE, sep="\t", na="NA", col.names=NA)

volcanoplot(fit, coef = "YF_HDMvPBS", style = "p-value")
#boxplot top 3 miRNAs
boxplot(Expvalues[rownames(YF_HDMvPBS[1,]),] ~ groups$group, main = rownames(YF_HDMvPBS[1,]), ylab = "miRNA expression", xlab=NULL, names=levels(groups$Treatment),las=2, col = ifelse(grepl("F",levels(groups$Treatment)),"pink","lightblue"))
boxplot(Expvalues[rownames(YF_HDMvPBS[2,]),] ~ groups$group, main = rownames(YF_HDMvPBS[2,]), ylab = "miRNA expression", xlab=NULL, names=levels(groups$Treatment),las=2, col = ifelse(grepl("F",levels(groups$Treatment)),"pink","lightblue"))
boxplot(Expvalues[rownames(YF_HDMvPBS[3,]),] ~ groups$group, main = rownames(YF_HDMvPBS[3,]), ylab = "miRNA expression", xlab=NULL, names=levels(groups$Treatment),las=2, col = ifelse(grepl("F",levels(groups$Treatment)),"pink","lightblue"))
#################
#let7i
# miRNAid[8]
#Used in manuscript
#boxplot(Expvalues[8,] ~ groups$group, col = c("pink","pink", "lightblue","lightblue"), main = "Let-7i expression in the 4 core genotypes", xlab ="", ylab = "miRNA expression", names=c('HDM XXF', 'PBS XXF', 'HDM XXM', 'PBS XXM','HDM XYF', 'PBS XYF','HDM XYM', 'PBS XYM'))
boxplot(Expvalues[8,] ~ groups$Treatment, col = ifelse(grepl("F",levels(groups$Treatment)),"pink","lightblue"), main = "Let-7i expression in the 4 core genotypes", xlab ="", ylab = "miRNA expression", names=levels(groups$Treatment))
boxplot(Expvalues[8,] ~ groups$Treatment, col = ifelse(grepl("F",levels(groups$Treatment)),"pink","lightblue"), main = "Let-7i expression in the 4 core genotypes", xlab ="", ylab = "miRNA expression", names=levels(groups$Treatment),las=2)
dev.off()
