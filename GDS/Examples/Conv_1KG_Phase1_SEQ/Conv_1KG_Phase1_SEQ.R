#########################################################################
##
##  Convert 1KG Phase 1 Genotypes to Sequencing GDS Format
##
##  File: Conv_1KG_Phase1_SEQ.R
##  Output: ALL_chr_phase1_release_v3_20101123_snps_indels_svs_genotypes.seq.gds
##

library(SeqArray)


# file name
FTP_BASE <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp"
FTP_GENO_PATH <- "release/20110521"
FTP_FILE_TEMPLATE <- "ALL.chr%s.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"



#########################################################################
##  Download the gz files from 1KG Website

for (chr in c(1:22, "X"))
{
	fn <- sprintf(FTP_FILE_TEMPLATE, chr)
	vcf.fn <- paste(FTP_BASE, FTP_GENO_PATH, fn, sep="/")
	download.file(vcf.fn, fn)
}



#########################################################################
##  Create GDS file

vcf.fn <- sprintf(FTP_FILE_TEMPLATE, c(1:22, "X"))
gds.fn <- "ALL_chr_phase1_release_v3_20101123_snps_indels_svs_genotypes.seq.gds"

seqVCF2GDS(vcf.fn, gds.fn, compress.option=seqCompress.Option("ZIP_RA.max"))



#########################################################################
##  Show the dataset

(f <- openfn.gds(gds.fn))
closefn.gds(f)



#########################################################################
##  Session Info

sessionInfo()

q("no")
