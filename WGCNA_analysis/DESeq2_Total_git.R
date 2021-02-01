#clean
rm(list = ls())

workingDir="~"
setwd(workingDir)

lnames = load(file = "4+3/Total_gut_signed.RData");
lnames
# TOM, datExpr, MEs, moduleLabels, moduleColors, geneTree

library(DESeq2)

# choose specific condition
# +1,+2,+3/-1,-2,-3 represent Citrus_spp._CLas+/-_midgut; the rest represents C._medica_CLas+/-_Gut
datExpr=as.data.frame(t(datExpr))
names(datExpr) <- c('+1','+2','+3','-1','-2','-3','+1.1','+1.2','+1.3','+1.4','-1.1','-1.2','-1.3','-1.4')
colData <- data.frame(sample=colnames(datExpr), 
                      condition=as.factor(c('positive','positive','positive',
                                            'negative','negative','negative',
                                            'positive','positive','positive','positive',
                                            'negative','negative','negative','negative')))
    

dds <- DESeqDataSetFromMatrix(countData = round(datExpr),
                              colData = colData,
                              design = ~ condition)
# filter count = 0 genes
dds <- dds[rowSums(counts(dds)) > 1,]

# PCA
rld <- rlog(dds)
vsd <- vst(dds)
plotPCA(rld, intgroup="condition")
plotPCA(rld, intgroup="sample")

# Dispersion estimate plot
plotDispEsts(ddsPvN, ylim = c(1e-6, 1e1), main="DESeq2")

# Differential analysis
ddsPvN=DESeq(dds)
res=results(ddsPvN, alpha = 0.01)
head(res)
summary(res)
write.table(res,"DESeq_mgresult.csv", sep = ",", row.names = TRUE)

# List significant DEGs
res_out <- as.data.frame(res[-which(is.na(res$padj)),])
res_out <- res_out[res_out$padj<0.01,]
write.table(res_out,"DESeq_Totalresult_sig<0.01.csv", sep = ",", row.names = TRUE)

# sort result by p-value
res_pval <- res[order(res$padj),]
head(res_pval)
write.table(res_pval,"DESeq_pvalmgresult.csv", sep = ",", row.names = TRUE)

# plot counts for top 6 differential expressed genes
par(mfrow=c(2,3))
plotCounts(ddsPvN, gene="Dcitr09g07140.1.2", intgroup="condition")
plotCounts(ddsPvN, gene="Dcitr10g06510.1.1", intgroup="condition")
plotCounts(ddsPvN, gene="Dcitr01g15840.1.1", intgroup="condition")
plotCounts(ddsPvN, gene="Dcitr12g10450.1.1", intgroup="condition")
plotCounts(ddsPvN, gene="Dcitr00g11910.1.1", intgroup="condition")
plotCounts(ddsPvN, gene="Dcitr06g08180.1.1", intgroup="condition")

# Volcano Plot
# reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue Down, red Up)
with(subset(res, padj<.01 & log2FoldChange<0 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & log2FoldChange>0), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

## heatmap for count matrix
sum(res$padj < 0.01, na.rm=TRUE)
library("pheatmap")
select=order(rowMeans(counts(dds,normalized=FALSE)),decreasing=TRUE)[1:1000]
nt=normTransform(dds) # defaults to log2(x+1)
df=as.data.frame(colData(dds)[,c("names","condition")])

pheatmap(assay(nt)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)

library(geneplotter)
plotMA(res, main="DESeq2", ylim=c(-2,2))


# Venn Diagrams for modules and DEGs

DEGs=rownames(res_out)
# list the module name and its intersect genes number with DEGs
df1=data.frame('module'=c(0),'intersect'=c(0))
for(i in 1:length(names(table(moduleColors)))){
  df1[i,1]=names(table(moduleColors))[i]
  modNames=names(datExpr)[which(moduleColors==df1[i,1])]
  df1[i,2]=length(intersect(DEGs,modNames))
}
df1=df1[with(df1, order(intersect,decreasing = T)),]

# output the genes that in modules which have the highest interception with DEGs
green=names(datExpr)[which(moduleColors=='green')]
blue=names(datExpr)[which(moduleColors=='blue')]
grey=names(datExpr)[which(moduleColors=='grey')]

# Load library
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(4, 'Accent')
# Chart has 3 modules and DEGs
venn.diagram(
  x = list(DEGs, green, blue, grey),
  category.names = c("DEGs" , "green module","blue module","grey module"),
  filename = 'blue_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 700 , 
  width = 750 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)
