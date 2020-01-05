#this is an expression computing script for the small RNA seq pipelines
args = commandArgs(trailingOnly=TRUE)

reads<-read.table(args[1])
lengths<-read.table(args[2])
star_file<-read.table(args[3],sep="\t",fill=TRUE)

uniq_read_number<-as.numeric(as.character(star_file[8,2]))

colnames(reads)<-c("tags","accession")
colnames(lengths)<-c("accession","length")

reads<-merge(reads,lengths,by="accession")
reads$rpkm<-reads$tags/(reads$length/1000)/(uniq_read_number/1E6)
reads$rpm<-reads$tags/(uniq_read_number/1E6)

write.table(reads,file=args[4],sep="\t",row.names=FALSE,quote=FALSE)


reads<-reads[which(reads$rpm>=10),]

write.table(reads,file=args[5],sep="\t",row.names=FALSE,quote=FALSE)
