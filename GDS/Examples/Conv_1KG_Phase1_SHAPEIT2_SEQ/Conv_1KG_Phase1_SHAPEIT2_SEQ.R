#########################################################################
##
##  Convert 1KG Phase 1 Genotypes to Sequence GDS Format
##
##  File: Conv_1KG_Phase1_SHAPEIT2_SEQ.R
##  Output: 1KG_autosomes_SHAPEIT2_integrated_phase1_v3_20101123_snps_indels_svs_genotypes.seq.gds
##

library(SeqArray)


# File name
FTP_BASE <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp"
FTP_GENO_PATH <- "phase1/analysis_results/shapeit2_phased_haplotypes"
FTP_FILE_TEMPLATE <- "ALL.chr%d.SHAPEIT2_integrated_phase1_v3.20101123.snps_indels_svs.genotypes.all.vcf.gz"



#########################################################################
##  Download the gz files from 1KG Website

for (chr in 1:22)
{
	fn <- sprintf(FTP_FILE_TEMPLATE, chr)
	vcf.fn <- paste(FTP_BASE, FTP_GENO_PATH, fn, sep="/")
	download.file(vcf.fn, fn)
}


#########################################################################
##  Create GDS file

vcf.fn <- sprintf(FTP_FILE_TEMPLATE, 1:22)
gds.fn <- "1KG_autosomes_SHAPEIT2_integrated_phase1_v3_20101123_snps_indels_svs_genotypes.seq.gds"

seqVCF2GDS(vcf.fn, gds.fn, compress.option=seqCompress.Option("ZIP_RA.max"))



#########################################################################
##  Show the dataset

(f <- openfn.gds(gds.fn))
closefn.gds(f)



#########################################################################
##  Session Info

sessionInfo()

q("no")
