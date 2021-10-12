workingDir="~"
setwd(workingDir)

library(WGCNA)
library(flashClust)
# increase the process speed
allowWGCNAThreads()

data=read.table("expression_file_PEN3.0.txt",header=F,stringsAsFactors = F,sep = "\t")
data=within(data, 
            V6<-data.frame(do.call('rbind', strsplit(as.character(V6), ',', fixed=TRUE))))

# choose sample condition
datExpr_mmGut=as.data.frame(t(
  data[which(data[,2]=='C._medica_CLas-_Gut'),]))[-c(2,3,4,5),]
names(datExpr_mmGut) <- lapply(datExpr_mmGut[1, ], as.character)
datExpr_mmGut <- datExpr_mmGut[-1,] 

datExpr_mpGut=as.data.frame(t(
  data[which(data[,2]=='C._medica_CLas+_Gut'),]))[-c(2,3,4,5),]
names(datExpr_mpGut) <- lapply(datExpr_mpGut[1, ], as.character)
datExpr_mpGut <- datExpr_mpGut[-1,] 

datExpr_mg=rbind(datExpr_mpGut,datExpr_mmGut)
datExpr_mg[] <- lapply(datExpr_mg, as.character)
datExpr_mg[] <- lapply(datExpr_mg, as.numeric)
rownames(datExpr_mg)<-c("C.medica_postive_1","C.medica_postive_2","C.medica_postive_3","C.medica_postive_4",
                        "C.medica_negative_1","C.medica_negative_2","C.medica_negative_3","C.medica_negative_4")
## data pre-processing
# check for missing data
gsg = goodSamplesGenes(datExpr_mg, verbose = 3)
gsg$allOK
if (!gsg$allOK) {
  
  # Optionally, print the gene and sample names that were removed: 
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr_mg)[!gsg$goodGenes], collapse = ", "))); 
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr_mg)[!gsg$goodSamples], collapse = ", "))); 
  # Remove the offending genes and samples from the data:
  datExpr_mg = datExpr_mg[gsg$goodSamples, gsg$goodGenes]
  
}
# form an adjacency matrix for samples
A=adjacency(t(datExpr_mg),type="distance")
# this measures the sample network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized conn
Z.k=scale(k)
# define the threshold to find outliers
# here means 5 standardized divation from other samples
thresholdZ.k=-5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
# calculate the cluster tree using flashClust
sampleTree = flashClust(as.dist(1-A), method = "average")

png("Sample clustering.png",
     width = 12, height = 9, units = "in", pointsize = 12,
     bg = "white", res = 100)
par(cex = 0.6)
plotDendroAndColors(sampleTree,colors=outlierColor,groupLabels = "Outliers",
                    autoColorHeight=F,colorHeight=0.05,
                    main = "Sample clustering to detect outliers")
dev.off()
# manual sample cut for C.medica gut negative

# filter out the data with great variance. 
# But since we are interested in samples variance which may influence the detox gene expression, it shouldn't be processed

#variances <- apply(datExpr, 2, var)
#min <- min(variances, na.rm=TRUE)
#max <- max(variances, na.rm=TRUE)
#cutoff <- min + ((max-min) / 2)
#filter <- which(apply(datExpr, 2, var) > cutoff)
#datExpr <- datExpr[,-filter]
datExpr=datExpr_mg
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

## pick soft threshold
powers = c(c(1:10),seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers,networkType = "signed" ,verbose = 5) 
# Plot the results:
png("C.medica_sft.png",
    width = 9, height = 5, units = "in", pointsize = 12,
    bg = "white", res = 500)

par(mfrow = c(1,2))
cex1 = 0.5
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red",cex=0.7)
text(1.5, 0.88, "0.9",col="red",cex=0.6)

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()
## one-step network construction
# eigengene correaltion higher than 0.75 will be merged
meringthreshold=0.25
net = blockwiseModules(datExpr, power = sft$powerEstimate,maxBlockSize = 20000,
                       corType = "bicor",
                       networkType = "signed",
                       mergeCutHeight = meringthreshold,
                       minModuleSize = 30, reassignThreshold = 0, 
                       numericLabels = TRUE, pamRespectsDendro = FALSE, 
                       saveTOMs = TRUE, saveTOMFileBase = "DcitriTOM", verbose = 3)


# important parameters
moduleLabels = net$colors 
moduleColors = labels2colors(net$colors) 
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Reorder MEs such that similar ones are next to each other.
MEs = orderMEs(MEs0)
geneTree = net$dendrograms[[1]]
power = sft$powerEstimate
#table(net$colors)

# Plot the dendrogram and the module colors underneath 
png("C.medica_modPlot.png",
    width = 9, height = 5, units = "in", pointsize = 12,
    bg = "white", res = 500)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], 
                    "Module colors", dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# visualizing ME and Plot the dendrogram & heatmap
png('C.medica_ME_dendrogram.png',
    width=11,height=8,units = "in", pointsize = 12,
    bg = "white", res = 500)
par(cex = 1.0) 
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()

png('C.medica_ME_heatmap.png',
    width=11,height=8,units = "in", pointsize = 12,
    bg = "white", res = 500)
par(cex = 1.0) 
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), 
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

# merge highly correlated genes
merge=mergeCloseModules(datExpr,moduleColors,cutHeight = 0.35)
MEsManual=merge$newMEs
moduleColorsManual=merge$colors
# show the effect of merging
datColors=data.frame(moduleColors,moduleColorsManual)

png('C.medica_modulePlot.png',
    width=20,height=12,units = "in", pointsize = 12,
    bg = "white", res = 500)
plotDendroAndColors(geneTree,colors =datColors, 
                    groupLabels = c("moduleColorsAutomatic","moduleColorsMerged"),
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

png('C.medica_MEPlot_automatic.png',
    width=11,height=8,units = "in", pointsize = 12,
    bg = "white", res = 500)
plotEigengeneNetworks(MEs, "", 
                      marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
dev.off()

# names (colors) of the modules
modNames = substring(names(MEs), 3)

# eigengene-based module connectivity
nSamples=nrow(datExpr)
datKME = signedKME(datExpr,MEs)
KMEPvalue = as.data.frame(corPvalueStudent(as.matrix(datKME), nSamples))
names(KMEPvalue) = paste("p.KME", modNames, sep="")

# plot module genes and KME(module menbership)
png('C.medica_KME_automatic_TBPG.png',
    width=11,height=8,units = "in", pointsize = 12,
    bg = "white", res = 500)
selectModules=c("turquoise","blue","pink","grey")
# c(4,6) 4 row 6 col
par(mfrow=c(2,2))
for(module in selectModules){
  column=match(module, modNames)
  moduleGenes= moduleColors==module
  plot(abs(datKME[moduleGenes, column]),
       col=module,
       main = paste(module,"module"),
       ylim=c(0,1),
       xlab = "genes size",
       ylab = "KME",
       cex.main = 1.2, cex.lab = 1.2, cex.axis = 0.8)
}
dev.off()

# barplot for MEs across samples
png('C.medica_ExpAcrossSample_automatic_TBP.png',
    width=30,height=20,units = "in", pointsize = 12,
    bg = "white", res = 500)
selectModules=c("MEturquoise","MEblue","MEpink")
# c(4,6) 4 row 6 col
par(mfrow=c(2,2),mar=c(15,29,4,2))
for(module in selectModules){
  barplot(MEs[,module],
          horiz=T,
          names.arg=row.names(MEs),
          axisnames=T,
          col=substring(module, 3),
          main = paste(substring(module, 3),"module",sep=" "),
          las=2,
          xlim=c(-0.5,0.6),
          xlab = "Module 1st PC expression level",
          cex.main = 3, cex.lab = 2, cex.axis = 1,cex=3)
}
dev.off()

# save parameters
save(datExpr,MEs,MEsManual,moduleLabels,moduleColorsManual, geneTree,power,
     file = "C._medica_gut_signed.RData")

# visualizing TOM which is used to get a sense of module detection 
TOM=TOMsimilarityFromExpr(datExpr,networkType = "signed",power=power)
dissTOM = 1- TOM
plotTOM = dissTOM^7;
diag(plotTOM) = NA;

png("TOM.png",
     width = 600, height = 600, units = "px", pointsize = 12,
     quality = 75,
     bg = "white", res = NA)
TOMplot(plotTOM, geneTree, moduleColorsManual, main = "Network heatmap plot, all genes")
dev.off()

# Multi-dimentional scaling
png("MDS.png",
    width = 11, height = 8, units = "in", pointsize = 12,
    bg = "white", res = 500)
cmd1=cmdscale(as.dist(dissTOM),2)
plot(cmd1, col=as.character(colorh1), main="MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")
dev.off()

## hub genes analysis and export to cytoscape
hubs=chooseTopHubInEachModule(datExpr,moduleColors,power=power)
hub_df=data.frame(hubs)
hub_df$freq=data.frame(table(moduleColors)[-which(names(table(moduleColors))=="grey")])[,2]
write.csv(hub_df,file="C.medica_hubgene.csv",quote = F)

modules = c("blue")
probes = names(datExpr) 
inModule = is.finite(match(moduleColorsManual, modules)); 
modProbes = probes[inModule]
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Top 30 hub genes
nTop = 30
IMConn = softConnectivity(datExpr[, modProbes])
top = (   (-IMConn) <= nTop)
cyt = exportNetworkToCytoscape(modTOM[top,top],
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""), 
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""), 
                               weighted = TRUE, 
                               threshold = 0.1, 
                               nodeNames = modProbes[top],
                               nodeAttr = moduleColorsManual[inModule][top])
# extract module transcriptiDs
totalColor=names(table(moduleColors))
for(modules in totalColor){
  inModule = is.finite(match(moduleColors, modules)); 
  modProbes = probes[inModule]
  write.table(modProbes,paste(modules,".txt",sep=""),col.names = F, row.names = F,quote=F)
}