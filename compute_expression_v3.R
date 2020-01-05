#this is an expression computing script for the small RNA seq pipelines
#updated 2019 07 17 to use a refflat for the chrom and strand info. 
args = commandArgs(trailingOnly=TRUE)

reads<-read.table(args[1])
lengths<-read.table(args[2])
star_file<-read.table(args[3],sep="\t",fill=TRUE)
refflat <- read.table(args[4],sep = "\t", head = FALSE, skip = 3)[,c(2,3,4)]
colnames(refflat) <- c("accession","chromosome","strand")
uniq_read_number<-as.numeric(as.character(star_file[8,2]))

colnames(reads)<-c("tags","accession")
colnames(lengths)<-c("accession","length")
reads <- merge(reads, refflat, by = "accession",all.x = TRUE)
reads<-merge(reads,lengths,by="accession")
reads$rpkm<-reads$tags/(reads$length/1000)/(uniq_read_number/1E6)
reads$rpm<-reads$tags/(uniq_read_number/1E6)


write.table(reads,file=args[5],sep="\t",row.names=FALSE,quote=FALSE)


reads<-reads[which(reads$rpm>=10),]

write.table(reads,file=args[6],sep="\t",row.names=FALSE,quote=FALSE)
