#FOR LUKAS
args <- commandArgs(trailingOnly=TRUE)
#convert gene names to symbol
#argument 1 !!in order!! is a file that contains an accession column. Argument 2 is the output file that will contain a gene symbol column. 

annot<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/symbol_to_accession.txt",head=TRUE)
input <- read.table(args[1],head=TRUE)
output <- write.table(merge(input,annot,by="accession"),file = args[2],quote=FALSE,sep="\t",row.names=FALSE)
