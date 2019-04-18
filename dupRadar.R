library(dupRadar)
args = commandArgs(trailingOnly=TRUE)

bamDuprm <- args[1]

gtf <- '~/Documents/referentie/refgene/hg19refGene.sorted.gtf'

# '0' (unstranded), '1' (stranded) and '2' (reversely stranded)

stranded <- 0
paired   <- FALSE
threads  <- 4

# Duplication rate analysis
dm <- analyzeDuprates(bamDuprm,gtf,stranded,paired,threads)

write.table(dm, args[2], sep='\t')


