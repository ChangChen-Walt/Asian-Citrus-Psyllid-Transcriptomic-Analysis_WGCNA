workingDir="~"
setwd(workingDir)
# If have the .RData file can directly reload the coefficients
#lnames = load(file = "Total_gut_signed.RData")
#lnames
library(WGCNA)
allowWGCNAThreads()
# data entry
data=read.table("expression_file_PEN3.0.txt",header=F,stringsAsFactors = F,sep = "\t")

# choose specific condition
# choose Citrus_spp._CLas-_midgut
datExpr_mGut=as.data.frame(data[which(data[,2]=='Citrus_spp._CLas-_midgut'),])
V6=datExpr_mGut[,6]
datExpr_mGut=within(datExpr_mGut, 
                     V6<-data.frame(do.call('rbind', strsplit(as.character(V6), ',', fixed=TRUE))))
datExpr_mGut=as.data.frame(t(datExpr_mGut)[-c(2,3,4,5),])
names(datExpr_mGut) <- lapply(datExpr_mGut[1, ], as.character)
datExpr_mGut <- datExpr_mGut[-1,] 
# choose Citrus_spp._CLas+_midgut
datExpr_pGut=as.data.frame(data[which(data[,2]=='Citrus_spp._CLas+_midgut'),])
V6=datExpr_pGut[,6]
datExpr_pGut=within(datExpr_pGut, 
                     V6<-data.frame(do.call('rbind', strsplit(as.character(V6), ',', fixed=TRUE))))
datExpr_pGut=as.data.frame(t(datExpr_pGut)[-c(2,3,4,5),])
names(datExpr_pGut) <- lapply(datExpr_pGut[1, ], as.character)
datExpr_pGut <- datExpr_pGut[-1,] 
# choose C._medica_CLas-_Gut
datExpr_mmGut=as.data.frame(data[which(data[,2]=='C._medica_CLas-_Gut'),])
V6=datExpr_mmGut[,6]
datExpr_mmGut=within(datExpr_mmGut, 
                     V6<-data.frame(do.call('rbind', strsplit(as.character(V6), ',', fixed=TRUE))))
datExpr_mmGut=as.data.frame(t(datExpr_mmGut)[-c(2,3,4,5),])
names(datExpr_mmGut) <- lapply(datExpr_mmGut[1, ], as.character)
datExpr_mmGut <- datExpr_mmGut[-1,] 
# choose C._medica_CLas+_Gut
datExpr_mpGut=as.data.frame(data[which(data[,2]=='C._medica_CLas+_Gut'),])
V6=datExpr_mpGut[,6]
datExpr_mpGut=within(datExpr_mpGut, 
                     V6<-data.frame(do.call('rbind', strsplit(as.character(V6), ',', fixed=TRUE))))
datExpr_mpGut=as.data.frame(t(datExpr_mpGut)[-c(2,3,4,5),])
names(datExpr_mpGut) <- lapply(datExpr_mpGut[1, ], as.character)
datExpr_mpGut <- datExpr_mpGut[-1,] 
# combine the 4 groups of replicates in order
datExpr_g=rbind(datExpr_pGut,datExpr_mGut,datExpr_mpGut,datExpr_mmGut)
datExpr_g[] <- lapply(datExpr_g, as.character)
datExpr_g[] <- lapply(datExpr_g, as.numeric)
row.names(datExpr_g) <-c('Mg+r1','Mg+r2','Mg+r3','Mg-r1','Mg-r2','Mg-r3',
                        'Cmg+r1','Cmg+r2','Cmg+r3','Cmg+r4','Cmg-r1','Cmg-r2','Cmg-r3','Cmg-r4')
# check for missing data
gsg = goodSamplesGenes(datExpr_g, verbose = 3)
gsg$allOK
if (!gsg$allOK) {
        
        # Optionally, print the gene and sample names that were removed: 
        if (sum(!gsg$goodGenes)>0)
                printFlush(paste("Removing genes:", paste(names(datExpr_g)[!gsg$goodGenes], collapse = ", "))); 
        if (sum(!gsg$goodSamples)>0)
                printFlush(paste("Removing samples:", paste(rownames(datExpr_g)[!gsg$goodSamples], collapse = ", "))); 
        # Remove the offending genes and samples from the data:
        datExpr_g = datExpr_g[gsg$goodSamples, gsg$goodGenes]
        
}
sampleTree = hclust(dist(datExpr_g), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 50000, col = "red");
# Determine cluster under the line
clust=cutreeStatic(sampleTree,cutHeight =10000)
#clust=cutree(sampleTree,h=10000)
table(clust)
# clust 1 contains the samples we want to keep
keepSamples = (clust==0)
# change variable
datExpr = datExpr_g[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# pick soft threshold
sizeGrWindow(9, 5) 
par(mfrow = c(1,2))
cex1 = 0.5
powers = c(c(1:10),seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 
# Plot the results:
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# datExpr must be numeric (signed/unsigned)
net = blockwiseModules(datExpr, power = sft$powerEstimate,
                       TOMType = "signed", minModuleSize = 30, reassignThreshold = 0, 
                       mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, 
                       saveTOMs = TRUE, saveTOMFileBase = "DcitriTOM", verbose = 3)
# open a graphics window 
sizeGrWindow(12, 9) 
# Convert labels to colors for plotting 
mergedColors = labels2colors(net$colors) 
# Plot the dendrogram and the module colors underneath 
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], 
                    "Module colors", dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)

# save important coefficient for further analysis
moduleLabels = net$colors 
moduleColors = labels2colors(net$colors) 
MEs = net$MEs
geneTree = net$dendrograms[[1]]
power = sft$powerEstimate
TOM = TOMsimilarityFromExpr(datExpr, power = power)
save(TOM,datExpr, MEs, moduleLabels, moduleColors, geneTree, power,file = "Total_gut_signed.RData")

# visualizing intermodular connectivity
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET = orderMEs(MEs)
sizeGrWindow(5,7.5); 
par(cex = 0.5) 
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM 
dissTOM=1-TOM
nSelect = 10000 
# For reproducibility, we set the random seed 
set.seed(10); 
select = sample(nGenes, size = nSelect); 
selectTOM = dissTOM[select, select]; 
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster. 
selectTree = hclust(as.dist(selectTOM), method = "average") 
selectColors = moduleColors[select]; # Open a graphical window sizeGrWindow(9,9) 
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot 
plotDiss = selectTOM^7; 
diag(plotDiss) = NA; 
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

# export to cytoscape
# TOM = TOMsimilarityFromExpr(datExpr, power = power)
modules = c("blue")
probes = names(datExpr) 
inModule = is.finite(match(moduleColors, modules)); 
modProbes = probes[inModule]
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
#Top 30 hub genes
nTop = 30
IMConn = softConnectivity(datExpr[, modProbes])
top = (rank(-IMConn) <= nTop)
cyt = exportNetworkToCytoscape(modTOM[top,top],
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""), 
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""), 
                               weighted = TRUE, 
                               threshold = 0.1, 
                               nodeNames = modProbes[top],
                               nodeAttr = moduleColors[inModule][top])
# extract module transcriptiDs
totalColor=names(table(moduleColors))
for(modules in totalColor){
        inModule = is.finite(match(moduleColors, modules)); 
        modProbes = probes[inModule]
        write.table(modProbes,paste(modules,".txt",sep=""),col.names = F, row.names = F,quote=F)
}


