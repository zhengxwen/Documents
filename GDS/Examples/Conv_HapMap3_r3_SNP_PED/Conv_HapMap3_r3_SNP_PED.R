#########################################################################
##
##  Convert HapMap3_r3 Genotypes to SNP GDS Format from PED format
##
##  File: Conv_HapMap3_r3_SNP.R
##  Output: HapMap3_r3_b36_fwd_consensus_qc_poly.gds
##

library(gdsfmt)
library(SNPRelate)


# file name
FTP_BASE <- "ftp://ftp.ncbi.nlm.nih.gov/hapmap/"
FTP_GENO_PATH <- "genotypes/2010-05_phaseIII/plink_format"
FTP_PED_FILE <- "hapmap3_r3_b36_fwd.consensus.qc.poly.ped.gz"
FTP_MAP_FILE <- "hapmap3_r3_b36_fwd.consensus.qc.poly.map.gz"

FTP_SAMPLE_PATH <- "genotypes/2010-05_phaseIII"
FTP_SAMPLE_ANNOT_FILE <- "relationships_w_pops_041510.txt"



#########################################################################
##  Download the gz files from HapMap Website

ped.fn <- paste(FTP_BASE, FTP_GENO_PATH, FTP_PED_FILE, sep="/")
map.fn <- paste(FTP_BASE, FTP_GENO_PATH, FTP_MAP_FILE, sep="/")
samp.fn <- paste(FTP_BASE, FTP_SAMPLE_PATH, FTP_SAMPLE_ANNOT_FILE, sep="/")

# download ...
download.file(ped.fn, FTP_PED_FILE)
download.file(map.fn, FTP_MAP_FILE)
download.file(samp.fn, FTP_SAMPLE_ANNOT_FILE)


head(read.table(FTP_SAMPLE_ANNOT_FILE, header=TRUE, stringsAsFactors=FALSE))



#########################################################################
##  Create GDS file

gds.fn <- "HapMap3_r3_b36_fwd_consensus_qc_poly.gds"

snpgdsPED2GDS(FTP_PED_FILE, FTP_MAP_FILE, gds.fn, family=FALSE,
	snpfirstdim=FALSE, compress.annotation="ZIP_RA.max",
	compress.geno="ZIP_RA.max", verbose=TRUE)


# add sample.annotation

gfile <- openfn.gds(gds.fn, FALSE)

sample.id <- read.gdsn(index.gdsn(gfile, "sample.id"))

dat <- read.table(FTP_SAMPLE_ANNOT_FILE, header=TRUE, stringsAsFactors=FALSE)
names(dat) <- c("family", "IID", "father", "mother", "sex", "phenotype", "population")
dat$sex[dat$sex==1] <- "M"
dat$sex[dat$sex==2] <- "F"
dat <- dat[match(sample.id, dat$IID),
	c("family", "father", "mother", "sex", "population")]
table(dat$population, exclude=NULL)
add.gdsn(gfile, "sample.annot", dat, compress="ZIP_RA.max", closezip=TRUE)

# show
gfile

# Close the GDS file
closefn.gds(gfile)

cleanup.gds(gds.fn)



#########################################################################
##  Summarize the dataset

snpgdsSummary(gds.fn)



#########################################################################
##  Session Info

sessionInfo()

q("no")
