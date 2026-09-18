args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(args[1]);out<-args[2]
suppressPackageStartupMessages(library(limma));suppressPackageStartupMessages(library(statmod))
stopifnot(as.character(getRversion())=='4.3.3',as.character(packageVersion('limma'))=='3.58.1',as.character(packageVersion('statmod'))=='1.5.0')
base<-as.matrix(read.delim(file.path(out,'protein_sample_scores.tsv'),row.names=1,check.names=FALSE));strict<-as.matrix(read.delim(file.path(out,'strict_matrix.tsv'),row.names=1,check.names=FALSE))
gmt<-strsplit(readLines(file.path(root,'data_raw/pathways/ReactomePathways.gmt')),'\t',fixed=TRUE);focus<-c('R-HSA-112315','R-HSA-112310','R-HSA-9907900','R-HSA-611105','R-HSA-72766','R-HSA-9612973','R-HSA-397014');targets<-c('SYP','SYT1','SNCA','NAA25','POMP','PSMA4','PSMB4','PSMC5','PSMD1','PSMD2','PSMD4','PSMG1');allfocus<-list();alltargets<-list();universe<-list()
writeTSV<-function(x,name)write.table(x,file.path(out,name),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
for(model in c('primary','strict_backbones','omit_NAA25_07')){
 x<-if(model=='strict_backbones')strict else base;if(model=='omit_NAA25_07')x<-x[,colnames(x)!='NAA25_07']
 g<-as.integer(grepl('NAA25',colnames(x)));design<-cbind(Intercept=1,Knockdown=g);elig<-rowSums(!is.na(x[,g==0]))>=4&rowSums(!is.na(x[,g==1]))>=4
 fit<-eBayes(lmFit(x[elig,],design),robust=TRUE,trend=FALSE);r<-topTable(fit,coef=2,number=Inf,sort.by='none',confint=TRUE);r$gene<-rownames(r);writeTSV(r,paste0('limma_',model,'.tsv'));z<-merge(data.frame(gene=targets),r,by='gene',all.x=TRUE,sort=FALSE);z$model<-model;alltargets[[model]]<-z
 y<-x[complete.cases(x),];idx<-lapply(gmt,function(z)which(rownames(y)%in%z[-c(1,2)]));names(idx)<-vapply(gmt,function(z)z[2],'');idx<-idx[lengths(idx)>=10&lengths(idx)<=500]
 universe[[model]]<-data.frame(model=model,tested_proteins=nrow(r),complete_proteins=nrow(y),pathways=length(idx))
 for(setting in if(model=='primary')c('fixed_0.01','estimated')else 'fixed_0.01'){
  cam<-camera(y,idx,design,contrast=2,inter.gene.cor=if(setting=='estimated')NA else .01);cam$pathway<-rownames(cam);if(!'Correlation'%in%names(cam))cam$Correlation<-.01;cam$model<-model;cam$setting<-setting;cam<-cam[,c('pathway','NGenes','Correlation','Direction','PValue','FDR','model','setting')];writeTSV(cam,paste0('camera_',model,'_',setting,'.tsv'));allfocus[[paste(model,setting)]]<-cam[cam$pathway%in%focus,]
 }
}
writeTSV(do.call(rbind,alltargets),'target_limma.tsv');writeTSV(do.call(rbind,allfocus),'focus_camera.tsv');writeTSV(do.call(rbind,universe),'model_universes.tsv');writeLines(trimws(capture.output(sessionInfo()),which='right'),file.path(out,'R_session.txt'))
print(do.call(rbind,universe))
