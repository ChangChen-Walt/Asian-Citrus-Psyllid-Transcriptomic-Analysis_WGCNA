root='~'
setwd(root)
Goanno=read.table('Dcitr_OGSv3.0_beta_pep.desc.fa.goanna.inverteb.exponly.percentident70_qcov70.out_goanna_gaf.tsv',
                  stringsAsFactors = F, header=F)

## seems have 1.2k replicates
Goanno=unique(Goanno)

## RData generated from C.medica_gutTranscriptome_WGCNA.R and C.medica_gutTranscriptome_DESeq2.R
lnames = load(file = "C._medica_gut_signed.RData");
lnames
lnames = load(file = "C._medica_DEGs.RData");
lnames

library(GO.db)

#totalColor=names(table(moduleColors))
totalColor=names(ls[-7])
probes = names(datExpr) 
module_info<-list()
for(modules in totalColor){
  inModule = is.finite(match(moduleColors, modules))
  modProbes = probes[inModule]
  module_info[[modules]]=modProbes
}

df=data.frame('DCid'=c(0),'GOid'=c(0))
extractGOid=function(module_info){
  GO_annol=list()
  for(i in 1:length(totalColor)){
    for(j in 1:length(module_info[[i]])){
      index=which(Goanno$V2%in%module_info[[totalColor[i]]][j])
      if(length(index)!=0){
        append=data.frame(cbind(rep(module_info[[i]][j],length(index)),Goanno[index,4]))
        colnames(append)=c('DCid','GOid')
        df=rbind(df,append) 
      }
    }
    df=unique(df)
    GO_annol=append(GO_annol,list(df))
    df=data.frame('DCid'=c(0),'GOid'=c(0))
  }
  return(GO_annol)
}
## list item names follow sequence of totalColor
GoannoL=extractGOid(module_info)

## module processing
GO_modlist=list()
GOannoIDls=list()
GOtotalIDls=list()


for(color in totalColor){
  i=which(totalColor==color)
  
  if(nrow(GoannoL[[i]])>2){
    Goanno_module=GoannoL[[i]][-1,] 
    annoted_ids=dim(table(Goanno_module[,1]))
    GOannoIDls[[color]]<-annoted_ids
    total_ids=length(module_info[[i]])
    GOtotalIDls[[color]]<-total_ids
    
    df_m=data.frame(sort(table(Goanno_module[,2]),decreasing = T))
## list module that has GOanno
    GO_modlist[[color]]<-df_m
    pdf(paste(color,'_module_GOanno.pdf'))
    p=barplot(df_m[,2],
              border = NA, axes = F,las=2,
              names.arg = df_m[,1],
              ylim=c(0,max(df_m[,2])+3),ylab="Go term Counts",cex.names=0.5)
    box(col = "grey62")
    title(main = paste(color,'module GO annotation results'), 
          sub = paste('annoted transcript ids =',annoted_ids,'/',total_ids))
    axis(side = 2, at = c(0,df_m[,2]), las = 1, col.axis = "grey62", col = "grey62", cex.axis = 0.8)
    dev.off()
    }
  else if(nrow(GoannoL[[i]])==2){
    Goanno_module=GoannoL[[i]][-1,] 
    annoted_ids=dim(table(Goanno_module[,1]))
    GOannoIDls[[color]]<-annoted_ids
    total_ids=length(module_info[[i]])
    GOtotalIDls[[color]]<-total_ids
    
    df_m=data.frame(table(Goanno_module[,2]))
    ## list module that has GOanno
    GO_modlist[[color]]<-df_m
    pdf(paste(color,'_module_GOanno.pdf'))
    p=barplot(df_m[,2],
              border = NA, axes = F,las=2,
              names.arg = df_m[,1],
              ylim=c(0,max(df_m[,2])+3),ylab="Go term Counts",cex.names=0.5)
    box(col = "grey62")
    title(main = paste(color,'module GO annotation results'), 
          sub = paste('annoted transcript ids =',annoted_ids,'/',total_ids))
    axis(side = 2, at = c(0,df_m[,2]), las = 1, col.axis = "grey62", col = "grey62", cex.axis = 0.8)
    dev.off()
  }
  
  else{
      GO_modlist[[color]]<-NA
      GOannoIDls[[color]]<-0
      
    } 
}

# export GOids for each module
for(i in 1:length(GO_modlist)){
  if(length(GO_modlist[[i]])!= 1){
    GO_modlist[[i]]$Term=c() 
    GO_modlist[[i]]$Aspect=c()
    module=GO_modlist[[i]][,1]
    for(GOid in module){
      GO_modlist[[i]]$Term[which(GOid==module)]<-Term(GOid)
      GO_modlist[[i]]$Aspect[which(GOid==module)]<-Ontology(GOid)
    }
    names(GO_modlist[[i]])=c("GOid","Freq","Term","Aspect")
    write.csv(GO_modlist[[i]],file = paste(names(GO_modlist[i]),'module GO annotation results.csv'),row.names = F)
    }
}

