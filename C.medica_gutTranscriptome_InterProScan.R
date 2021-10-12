root='~'
setwd(root)
IPRS_raw=read.table('InterProScan/Dcitr_OGSv3_acc_go_counts.txt',
                  sep="\t",fill = T,header=T,stringsAsFactors = F)
GOanna_raw=read.table('GOanna/Dcitr_OGSv3.0_beta_pep.desc.fa.goanna.inverteb.exponly.percentident70_qcov70.out_goanna_gaf_rep.tsv',
                  sep="\t",fill = T,header=F,stringsAsFactors = F)
GOanna_raw=unique(GOanna_raw)


workingDir="~"

## RData generated from C.medica_gutTranscriptome_WGCNA.R and C.medica_gutTranscriptome_DESeq2.R
setwd(workingDir)
lnames = load(file = "C._medica_DEGs.RData");
lnames
lnames = load(file = "C._medica_gut_signed.RData");
lnames

## GOanna formating
GOanna=by(GOanna_raw[,5],INDICES = GOanna_raw[,2],FUN =paste)
df=data.frame(names(GOanna))
for(i in 1:length(df[,1])){
  df[i,2]<-paste(GOanna[[df[i,1]]],collapse = ",")
}


## change OGSv3 ID to uniform format in InterProScan
for(i in 1:nrow(IPRS_raw)){
  for(j in 1:length(IPRS_raw$Accession_ID[i])){
    IPRS_raw$Accession_ID[i]=paste0(substr(IPRS_raw$Accession_ID[i],1,13),'.',
                                  substr(IPRS_raw$Accession_ID[i],14,14),'.',
                                  substr(IPRS_raw$Accession_ID[i],15,15))
  }
}

IPRS=IPRS_raw[,c(1,3,5,7)]
IPRS$GO=paste0(IPRS$BP_GO_IDs,IPRS$MF_GO_IDs,IPRS$CC_GO_IDs)
IPRS$GO=gsub(";","",IPRS$GO)
IPRS$GO=gsub("G",",G",IPRS$GO)
IPRS$GO=sub('.', '', IPRS$GO)
IPRS=IPRS[,c(1,5)]

## combine IPRS and GOanna
names(df)<-names(IPRS)
comb=rbind(IPRS,df)
GOMaps_ls=by(comb[,2],INDICES = comb[,1],FUN = paste)
df_2=data.frame(names(GOMaps_ls),stringsAsFactors = F)
for(i in 1:length(df_2[,1])){
  if(length(GOMaps_ls[[i]])>1){
    df_2[i,2]<-paste0(unique(unlist(strsplit(paste0(GOMaps_ls[[i]],collapse = ","),split = ","))),collapse = ",")
  }
  else{
    df_2[i,2]<-GOMaps_ls[[i]]
  }
}
write.table(df_2,file="Gene2GOanna.map",row.names = F,col.names = F,sep="\t",quote = F)

## GOslim
library(tidyr)
names(df_2)<-c("geneId","GOId")
df_3=separate_rows(df_2,GOId,sep = ",")
write.table(df_3,file="Gene2GO",row.names = F,col.names = F,sep="\t",quote = F)

for(color in names(table(moduleColors))){
  colorG=names(datExpr)[which(moduleColors==color)]
  colorGO=df_3[which(df_3$geneId%in%colorG),]
  write.table(colorGO,file=paste(color,"GO"),row.names = F,col.names = F,sep="\t",quote = F)
}

library(topGO)
geneID2GO <- readMappings(file="Gene2GOanna.map", sep = "\t", IDsep = ",")
geneNames <- names(geneID2GO)

myInterestingGenes=names(datExpr)[which(moduleColors=="blue")]
#myInterestingGenes=row.names(res_out)[which(res_out$log2FoldChange<0)]
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames

GOdata <- new("topGOdata",ontology="CC",allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GOdata, classic = resultFis,orderBy = "weight", ranksOf = "classic", topNodes = 30)
allRes$classic<-as.numeric(allRes$classic)
allRes <- allRes[allRes$classic < 0.05,] # filter terms for p<0.05

write.csv(allRes,file="C.medica_turquoise_GO_CC_top40.csv",row.names = F,quote = F)

png("Significant GO term_C.medica_turquoise_BP.png",
    width = 12, height = 9, units = "in", pointsize = 12,
    bg = "white", res = 500)
par(cex=1.5)
showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 5, useInfo = 'all')
dev.off()

# GO bubble plot

library(ggplot2)
library(forcats)

blue_CC=data.frame(allRes[,c(2,3,4,6)],"type"=rep("CC",nrow(allRes)))
blue_BP=data.frame(allRes[,c(2,4,6)],"type"=rep("BP",nrow(allRes)))
blue_MF=data.frame(allRes[,c(2,4,6)],"type"=rep("MF",nrow(allRes)))
mydat=blue_CC[seq(from=1,to=15),]

names(mydat)=c("Term","Genome.annotated","Module.annotated","p_value","type")
png("C.medica_blue.png",
    width = 11, height = 8, units = "in", pointsize = 12,
    bg = "white", res = 500)
ggplot(mydat, aes(y = Term, x = -log(p_value), size = Module.annotated/Genome.annotated)) + 
  geom_point(aes(color = p_value), alpha = 1.0)+
  theme(text = element_text(size=20))
dev.off()


