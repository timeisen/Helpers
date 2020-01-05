#test if normal
args<-commandArgs(trailingOnly=TRUE)

data<-read.table(args[1])
colnames(data)<-c("accession","identifier","state","tail_length")
gene<-unique(data$accession)


l<-sapply(1:length(gene),function(x){return(data[which(data$accession==gene[x]),4])})
names(l) <- gene
param_df<-unlist(sapply(l,function(x){if(length(x)<5000 & length(x)>100){return(shapiro.test(x)$p)}}))

write.table(param_df,file=args[2],quote=FALSE,col.names=FALSE)