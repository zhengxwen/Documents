#########################################################################
##
##  Convert HapMap2_r24 Genotypes to GDS Format
##  File: Conv_HapMap2_r24.R
##  Output: HapMap2_r24_nr_b36_fwd.gds
##

library(gdsfmt)
library(SNPRelate)


# file name template
FILE_TEMPLATE <- "genotypes_chr%s_%s_r24_nr.b36_fwd.txt.gz"



#########################################################################
##  Download the gz files from HapMap Website

FTP_BASE <- "ftp://ftp.ncbi.nlm.nih.gov/hapmap/"
FTP_GENO_PATH <- "genotypes/2008-10_phaseII/fwd_strand/non-redundant"

for (chr.id in c(1:22, "X", "Y", "M"))
{
	for (pop in c("CEU", "JPT+CHB", "YRI"))
	{
		fn <- sprintf(FILE_TEMPLATE, chr.id, pop)
		ftp.fn <- file.path(FTP_BASE, FTP_GENO_PATH, fn)
		if (!file.exists(fn))
			download.file(ftp.fn, fn)
	}
}

# E.g., genotypes_chr1_CEU_r24_nr.b36_fwd.txt.gz:
#
# Col1: refSNP rs# identifier at the time of release (NB might merge 
#       with another rs# in the future)
# Col2: SNP alleles according to dbSNP
# Col3: chromosome that SNP maps to 
# Col4: chromosome position of SNP, in basepairs on reference sequence
# Col5: strand of reference sequence that SNP maps to
# Col6: version of reference sequence assembly (currently NCBI build36)
# Col7: HapMap genotyping center that produced the genotypes
# Col8: LSID for HapMap protocol used for genotyping
# Col9: LSID for HapMap assay used for genotyping
# Col10: LSID for panel of individuals genotyped
# Col11: QC-code, currently 'QC+' for all entries (for future use)
# Col12 and on: observed genotypes of samples, one per column, sample
#       identifiers in column headers (Coriell catalog numbers, example:
#       NA10847). Duplicate samples have .dup suffix.


# download the sample annotation file
FTP_SAMP_PATH <- "samples_individuals"
sample.fn <- "relationships_w_pops_121708.txt"
download.file(file.path(FTP_BASE, FTP_SAMP_PATH, sample.fn), sample.fn)

head(read.table(sample.fn, header=TRUE, stringsAsFactors=FALSE))



#########################################################################
##  Get HapMap individual IDs

SampID_File <- function(chr.id)
{
	ans <- c()
	for (pop in c("CEU", "JPT+CHB", "YRI"))
	{
		fn <- sprintf(FILE_TEMPLATE, "1", pop)
		tab <- read.table(fn, header=TRUE, comment.char="",
			stringsAsFactors=FALSE, nrows=2)
		ans <- c(ans, colnames(tab)[-(1:11)])
	}
	ans
}

sample.id <- SampID_File("1")

# check individual orders
for (chr.id in c(2:22, "X", "Y", "M"))
{
	s <- SampID_File(chr.id)
	if (!identical(sample.id, s))
		stop("Individual IDs are not consistent among chromosomes.")
	cat(sprintf("Chromosome %s: OK\n", chr.id))
}



#########################################################################
##  Create GDS file

# Create a new GDS file
gds.fn <- "HapMap2_r24_nr_b36_fwd.gds"
newfile <- createfn.gds(gds.fn)

# Add a format flag
put.attr.gdsn(newfile$root, "FileFormat", "SNP_ARRAY")

# Add variables
add.gdsn(newfile, "sample.id", sample.id, compress="ZIP_RA.max",
	closezip=TRUE)
var.snp <- add.gdsn(newfile, "snp.id", valdim=0L, storage="int",
	compress="ZIP_RA.max")
var.rsid <- add.gdsn(newfile, "snp.rs.id", valdim=0L, storage="string",
	compress="ZIP_RA.max")
var.chr <- add.gdsn(newfile, "snp.chromosome", valdim=0L, storage="string",
	compress="ZIP_RA.max")
var.pos <- add.gdsn(newfile, "snp.position", valdim=0L, storage="int",
	compress="ZIP_RA.max")
var.allele <- add.gdsn(newfile, "snp.allele", valdim=0L, storage="string",
	compress="ZIP_RA.max")

# Add genotypes
var.geno <- add.gdsn(newfile, "genotype", valdim=c(length(sample.id), 0),
	storage="bit2", compress="ZIP_RA.max")

# Indicate the SNP matrix is sample-by-snp
put.attr.gdsn(var.geno, "sample.order")

# initial var.
n.int <- 0L

cat("Starting ...\n")

# Write SNPs into the GDS file
for (chr.id in c(1:22, "X", "Y", "M"))
{
	v <- list()
	for (pop in c("CEU", "JPT+CHB", "YRI"))
	{
		fn <- sprintf(FILE_TEMPLATE, chr.id, pop)
		cat("Read:", fn)
		v[[pop]] <- read.table(fn, header=TRUE, comment.char="",
			stringsAsFactors=FALSE)
		cat(",", dim(v[[pop]]), "\n")
	}

	common.pos <- intersect(v[["CEU"]]$pos, v[["JPT+CHB"]]$pos)
	common.pos <- intersect(common.pos, v[["YRI"]]$pos)
	common.pos <- common.pos[order(common.pos)]
	cat("\tPosition Intersect:", length(common.pos), "\n")
	for (pop in c("CEU", "JPT+CHB", "YRI"))
		v[[pop]] <- v[[pop]][match(common.pos, v[[pop]]$pos), ]

	common.rs. <- intersect(v[["CEU"]]$rs., v[["JPT+CHB"]]$rs.)
	common.rs. <- intersect(common.rs., v[["YRI"]]$rs.)
	cat("\tRS ID Intersect:", length(common.rs.), "\n")
	for (pop in c("CEU", "JPT+CHB", "YRI"))
		v[[pop]] <- v[[pop]][match(common.rs., v[[pop]]$rs.), ]

	stopifnot(identical(v[["CEU"]]$SNPalleles, v[["JPT+CHB"]]$SNPalleles))
	stopifnot(identical(v[["CEU"]]$SNPalleles, v[["YRI"]]$SNPalleles))
	cat("\tSNP alleles:")
	print(table(v[["CEU"]]$SNPalleles))

	# detect and remove non-standard alleles
	allele <- v[["CEU"]]$SNPalleles
	s <- strsplit(allele, "/")
	flag <- sapply(s, function(x)
		(length(x)==2L) & all(x %in% c("A", "G", "C", "T")))
	if (sum(!flag) > 0)
	{
		for (pop in c("CEU", "JPT+CHB", "YRI"))
			v[[pop]] <- v[[pop]][flag, ]
		allele <- v[["CEU"]]$SNPalleles
		cat(sprintf("\tRemove non-standard alleles (%d):", sum(!flag)))
		print(table(v[["CEU"]]$SNPalleles))
	}

	# snp.id
	n <- dim(v$CEU)[1]
	append.gdsn(var.snp, seq_len(n) + n.int)
	n.int <- n.int + n

	# snp.rs.id
	append.gdsn(var.rsid, v[["CEU"]]$rs.)

	# snp.chromosome
	append.gdsn(var.chr, rep(chr.id, n))

	# snp.position
	append.gdsn(var.pos, v[["CEU"]]$pos)

	# snp.allele
	append.gdsn(var.allele, allele)

	s <- strsplit(allele, "/")
	g0  <- sapply(s, function(x) paste0(x[2],x[2]))
	g1a <- sapply(s, function(x) paste0(x[1],x[2]))
	g1b <- sapply(s, function(x) paste0(x[2],x[1]))
	g2  <- sapply(s, function(x) paste0(x[1],x[1]))
	gNN <- rep("NN", length(allele))

	# genotype
	m <- cbind(as.matrix(v[["CEU"]][, -c(1:11)]),
		as.matrix(v[["JPT+CHB"]][, -c(1:11)]),
		as.matrix(v[["YRI"]][, -c(1:11)]))
	f0  <- (m == g0)
	f1a <- (m == g1a)
	f1b <- (m == g1b)
	f2  <- (m == g2)
	fNN <- (m == gNN)
	if (!all(f0 | f1a | f1b | f2 | fNN))
		stop("Invalid SNP genotypes")

	g <- 2L*f2 + f1a + f1b
	g[fNN] <- 3L
	cat("\tGenotypes:")
	print(table(c(g)))

	append.gdsn(var.geno, t(g))
	cat("\n")
}


# Add sample.annotation

dat <- read.table(sample.fn, header=TRUE, stringsAsFactors=FALSE)
dat$sex[dat$sex==1] <- "M"
dat$sex[dat$sex==2] <- "F"
dat <- dat[match(sample.id, dat$IID),
	c("FID", "IID", "dad", "mom", "sex", "population")]
table(dat$population, exclude=NULL)
add.gdsn(newfile, "sample.annot", dat, compress="ZIP_RA.max", closezip=TRUE)

# Close the GDS file
closefn.gds(newfile)

cleanup.gds(gds.fn)



#########################################################################
##  Summarize the dataset

snpgdsSummary(gds.fn)



#########################################################################
##  Session Info

sessionInfo()

q("no")
