#This script takes genes that are output from a model script and returns a file with a header that contains either the minimum residual or the mean residual and associated constants. 
library(plyr)
options("warn"=1)
print("min_residual or confidence_intervals or avg?")
args<-commandArgs(trailingOnly=TRUE)

m<-read.table(args[1])

min_residual <- function(x,m){
	m_temp <- m[which(m$accession %in% muniq[x]),]
	return(m_temp[which.min(m_temp$residual),])
}

avg <- function(x,m){
	m_temp <- m[which(m$accession %in% muniq[x]),]
	m_temp[,-1] = apply(m_temp[,-1],2,mean)
	return(m_temp)
}

confidence_intervals <- function(x,m){
	m_temp <- m[which(m$accession %in% muniq[x]),]
	accession <- as.character(m_temp[1,1])
	m_temp <- apply(m_temp[,-1],2,function(x){(1-min(x)/max(x))*100})
	m_temp <- c(accession,m_temp)
	return(m_temp)}


colnames(m)<-c("accession","st","a","k","b","residual")
muniq<-unique(m$accession)

if(args[3]=="min_residual"){
	data<-lapply(muniq,function(x){min_residual(x,m)})
	all_data <- ldply(data, data.frame)}

if(args[3]=="avg"){
	data<-lapply(muniq,function(x){min_residual(x,m)})
	all_data <- ldply(data, data.frame)}

if(args[3]=="confidence_intervals"){
	data<-unlist(lapply(muniq,function(x){confidence_intervals(x,m)}))
	all_data<-as.data.frame(matrix(data,ncol=6,byrow=TRUE))
	colnames(all_data)<-c("accession","st_conf","a_conf","k_conf","b_conf","residual_conf")}

all_data<-all_data[complete.cases(all_data),]

write.table(all_data,file=args[2],sep="\t",row.names=FALSE,quote=FALSE)

