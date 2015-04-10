#########################################################################
##
##  Convert 1KG Phase 1 Genotypes to SNP GDS Format
##
##  File: Conv_1KG_Phase1_SNP.R
##  Output: ALL_chr_phase1_release_v3_20101123_snps_indels_svs_genotypes.snp.gds
##

library(gdsfmt)
library(SNPRelate)


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
gds.fn <- "ALL_chr_phase1_release_v3_20101123_snps_indels_svs_genotypes.snp.gds"

snpgdsVCF2GDS(vcf.fn, gds.fn, method="copy.num.of.ref", snpfirstdim=FALSE,
	compress.annotation="ZIP_RA.max", compress.geno="ZIP_RA.max")



#########################################################################
##  Summarize the dataset

snpgdsSummary(gds.fn)



#########################################################################
##  Session Info

sessionInfo()

q("no")
