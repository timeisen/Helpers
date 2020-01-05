##Adds a header column to the single tails files.
## TJE 2019 07 17
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
InFile = read_tsv(args[1],col_names = FALSE)
colnames(InFile) = c("accession",paste0("TailLength",0:250))
OutFile = write_tsv(InFile, path = args[2])