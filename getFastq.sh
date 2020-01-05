#!/usr/bin/env bash

#Code originally from Justin via Jamie
#Modified by TJE to 
#2019 04 02

var=$1

bsub -q 18 -J "wget" wget "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/"${var:0:3}/${var:0:6}/$var/$var.sra
bsub -q 18 -J "fastq-dump" -w "ended(wget)" fastq-dump $var.sra
