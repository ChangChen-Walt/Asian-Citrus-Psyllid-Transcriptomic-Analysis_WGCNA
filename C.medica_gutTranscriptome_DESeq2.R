workingDir="~"
setwd(workingDir)

library(DESeq2)
data=read.table("expression_file_PEN3.0.txt",header=T,stringsAsFactors = F,sep = "\t")
data=within(data, 
            Expression<-data.frame(do.call('rbind', strsplit(as.character(Expression), ',', fixed=TRUE))))
data=cbind(data[,-6],as.matrix(data[,6]))

# choose specific condition
datExpr_mmGut=as.data.frame(
  data[which(data[,2]=='C._medica_CLas-_Gut'),])[,-c(2,3,4,5)]

datExpr_mpGut=as.data.frame(
  data[which(data[,2]=='C._medica_CLas+_Gut'),])[,-c(2,3,4,5)]

datExpr_Gutcnt=cbind(datExpr_mmGut[,-1],as.matrix(datExpr_mpGut[,-1]))
datExpr_Gutcnt[] <- lapply(datExpr_Gutcnt, as.character)
datExpr_Gutcnt[] <- lapply(datExpr_Gutcnt, as.numeric)

rownames(datExpr_Gutcnt)<- datExpr_mmGut[,1]
names(datExpr_Gutcnt) <- c('-1','-2','-3','-4','+1','+2','+3','+4')

colData <- data.frame(names=colnames(datExpr_Gutcnt), 
                      condition=as.factor(c('negative','negative','negative','negative',
                                             'positive','positive','positive','positive')))
    

dds <- DESeqDataSetFromMatrix(countData = round(datExpr_Gutcnt),
                              colData = colData,
                              design = ~ condition)
# filter count = 0 genes
dds <- dds[rowSums(counts(dds)) > 1,]

# PCA
rld <- rlog(dds)
vsd <- vst(dds)
plotPCA(rld, intgroup="condition")

# Dispersion estimate plot
plotDispEsts(ddsPvN, ylim = c(1e-6, 1e1), main="DESeq2 ")

# Differential analysis
ddsPvN=DESeq(dds)
res=results(ddsPvN, alpha = 0.01)
head(res)
summary(res)
res0.01sortedPadj=res[order(res$padj),]
write.table(res,"DESeq_result.csv", sep = ",", row.names = TRUE)

# sort result by p-value
res_pval <- res[order(res$padj),]
head(res_pval)
write.table(res_pval,"DESeq_pvalresult.csv", sep = ",", row.names = TRUE)
# plot counts for top 6 differential expressed genes
par(mfrow=c(2,3))
plotCounts(ddsPvN, gene="Dcitr03g03290.1.2", intgroup="condition")
plotCounts(ddsPvN, gene="Dcitr06g03890.1.1", intgroup="condition")
plotCounts(ddsPvN, gene="Dcitr09g02200.1.1", intgroup="condition")
plotCounts(ddsPvN, gene="Dcitr03g04320.1.1", intgroup="condition")
plotCounts(ddsPvN, gene="Dcitr05g06760.1.1", intgroup="condition")
plotCounts(ddsPvN, gene="Dcitr00g12490.1.1", intgroup="condition")

# List significant DEGs
res_out <- as.data.frame(res[-which(is.na(res$padj)),])
res_out <- res_out[res_out$padj<0.01,]
DEGs=row.names(res_out)

write.table(res_out,"DESeq_Cm_result_sig<0.01.csv", sep = ",", row.names = TRUE)

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
select=order(rowMeans(counts(dds,normalized=FALSE)),decreasing=TRUE)[1:725]
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
write.table(df1,"DESeq_Cm_modulewise.csv", sep = ",", row.names = F)

# bar chart for intersection
mtx=as.matrix(df1)
NonemtyModule_index=which(as.numeric(mtx[,2])>0)

pdf('C.medica module DEGs comparison.pdf')
p=barplot(as.numeric(mtx[,2])[NonemtyModule_index],names.arg=mtx[NonemtyModule_index,1],
          axes = F,las=2,ylim = c(0,500),cex.names=0.8)
box(col = "grey62")
title(main = 'C.medica module DEGs comparison')
axis(side = 2, at = c(0,as.numeric(mtx[,2])), las = 1, col.axis = "grey62", col = "grey62", cex.axis = 0.8)
dev.off()

# output the genes that in modules which have the highest interception with DEGs
yellow=names(datExpr)[which(moduleColors=='yellow')]
blue=names(datExpr)[which(moduleColors=='blue')]
turquoise=names(datExpr)[which(moduleColors=='turquoise')]
pink=names(datExpr)[which(moduleColors=='pink')]
skyblue3=names(datExpr)[which(moduleColors=='skyblue3')]


for(color in names(table(moduleColors))){
  colorG=names(datExpr)[which(moduleColors==color)]
  write(colorG,file=color,sep = "\t")
}


VennMtx=df1[which(as.numeric(df1[,2])>0),]
ls=c()
for(i in 1:nrow(VennMtx) ){
  ls[[VennMtx[i,1]]]=names(datExpr)[which(moduleColors==VennMtx[i,1])]
}
ls[["DEGs"]]=DEGs

library(nVennR)
png("Venn_C.medica.png",
    width = 11, height = 8, units = "in", pointsize = 12,
    bg = "white", res = 500)
Venn <- plotVenn(ls[-4],nCycles=2000)
showSVG(Venn,setColors=c("blue","turquoise","yellow","pink","skyblue","black"),outFile ="Venn_C.medica.svg")
dev.off()

library(UpSetR)
png("Venn_C.medica.png",
    width = 11, height = 8, units = "in", pointsize = 12,
    bg = "white", res = 500)
upset(fromList(ls[-4]), nsets = 6, 
      sets.bar.color =c("turquoise","blue", "gold","black","pink","skyblue3"),
      main.bar.color =c("turquoise","blue","gold","black","black","pink","skyblue3","black","black","black","black"),
      text.scale = c(1.3, 1.3, 1, 1, 2, 2),
      order.by = "freq")
dev.off()

# Load library
library(VennDiagram)
# Chart has 3 modules and DEGs
venn.diagram(
  x = ls[-c(4,6)],
  category.names = names(ls[-c(4,6)]),
  filename = 'C.medica_Venn.png',
  main="Venn diagram for module DEGs",
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 700 , 
  width = 900 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c(names(ls[-c(4,6)])[-5],"black"),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)

save(res_out,ls,
     file = "C._medica_DEGs.RData")

## barchart for module size
df2=data.frame('module'=c(0),'size'=c(0))
for(i in 1:length(names(table(moduleColors)))){
  df2[i,1]=names(table(moduleColors))[i]
  df2[i,2]=table(moduleColors)[i]
}
df2=df2[with(df2, order(size,decreasing = T)),]
#write.table(df2,"WGCNA_Cm_moduleSize.csv", sep = ",", row.names = F)

mtx2=as.matrix(df2)
SmallModule_index=which(as.numeric(mtx2[,2])<100)
pdf('Citron module size comparison.pdf')
p=barplot(as.numeric(mtx2[,2]),names.arg=mtx2[,1],
          axes = F,las=2,ylim = c(0,4910),cex.names=0.5)
p1=abline(h=100,col = "red")
box(col = "grey62")
title(main = 'Citron module size comparison',sub= paste(length(SmallModule_index),'Modules size < 100'))
axis(side = 2, at = c(0,as.numeric(mtx2[,2])), las = 1, col.axis = "grey62", col = "grey62", cex.axis = 0.8)
dev.off()

## Build GO profiles of small (<100) modules with DE genes
SmallModuleBin=df2[SmallModule_index,1]
NonemtyModuleBin=df1[NonemtyModule_index,1]
# no DE in small modules
length(intersect(NonemtyModuleBin,SmallModuleBin))

