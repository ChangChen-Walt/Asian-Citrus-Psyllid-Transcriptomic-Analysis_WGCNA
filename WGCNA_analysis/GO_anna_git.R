#clean
rm(list = ls())

root='~'
setwd(root)
Goanna=read.table('Dcitr_OGSv3.0_beta_pep.desc.fa.goanna.inverteb.exponly.percentident70_qcov70.out_goanna_gaf copy.tsv',
                  stringsAsFactors = F, header=F)

## seem has 1.2k replicates
Goanna=unique(Goanna)

lnames = load(file = "Total_gut_signed.RData");
lnames

# create a list contains each module's info
totalColor=names(table(moduleColors))
probes = names(datExpr) 
module_info<-list()
for(modules in totalColor){
  inModule = is.finite(match(moduleColors, modules))
  modProbes = probes[inModule]
  module_info[[modules]]=modProbes
}

df=data.frame('DCid'=c(0),'GOid'=c(0))
# Since the Goanna data follows a format: first mapped all the GO terms which are associated with one Gene, then move to map
# the next Gene. Hence, the number of GO terms one Gene obtaining can be measured by the length of the return index vector 
# using "which()" function.
extractGOid=function(module_info){
  GO_annal=list()
  for(i in 1:length(totalColor)){
    for(j in 1:length(module_info[[i]])){
      index=which(Goanna$V2%in%module_info[[totalColor[i]]][j])
      if(length(index)!=0){
        append=data.frame(cbind(rep(module_info[[i]][j],length(index)),Goanna[index,4]))
        colnames(append)=c('DCid','GOid')
        df=rbind(df,append) 
      }
    }
    df=unique(df)
    GO_annal=append(GO_annal,list(df))
    df=data.frame('DCid'=c(0),'GOid'=c(0))
  }
  return(GO_annal)
}
## list item names follow sequence of totalColor
GoannoL=extractGOid(module_info)

## module processing
intermod_intersec=c()
GO_modlist=list()
GOannaIDls=list()
GOtotalIDls=list()
module_color=c(as.character(names(table(moduleColors))))


for(color in module_color){
  i=which(totalColor==color)
  # this returns an overview of the Goanna results in each module
  if(nrow(GoannoL[[i]])!=1){
    Goanno_module=GoannoL[[i]][-1,] 
    annoted_ids=dim(table(Goanno_module[,1]))
    GOannaIDls[[color]]<-annoted_ids
    total_ids=length(module_info[[i]])
    GOtotalIDls[[color]]<-total_ids
    
    df_m=data.frame(sort(table(Goanno_module[,2]),decreasing = T))
## list module that has GOanna
    GO_modlist[[color]]<-df_m
    intermod_intersec=append(intermod_intersec,as.character(df_m[,1]))
    
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
    }else{
      GO_modlist[[color]]<-NA
      GOannaIDls[[color]]<-0
    }
}
## plot for whole intermodule
df_im=data.frame(sort(table(intermod_intersec),decreasing = T))
GO_modlist=append(GO_modlist,list('summary'=df_im))
## save modlist
save(GO_modlist,GOannaIDls,GOtotalIDls,file = "C._medica_gut_box1GO.RData")

pdf(paste('box2_module_GOanno.pdf'))
p=barplot(df_im[,2],
          border = NA, axes = F,las=2,
          names.arg = df_im[,1],
          ylim=c(0,max(df_im[,2])+3),ylab="intermodule Go term intersection",cex.names=0.5)
box(col = "grey62")
title(main = 'box2 module GO annotation results', 
      sub = paste('total go terms =',nrow(df_im)))
axis(side = 2, at = c(0,df_im[,2]), las = 1, col.axis = "grey62", col = "grey62", cex.axis = 0.8)
dev.off()
