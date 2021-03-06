```r
#########################################################################
##
##  Convert HapMap3_r3 Genotypes to SNP GDS Format
##  File: Conv_HapMap3_r3.R
##  Output: HapMap3_r3_b36_fwd_consensus_qc_poly.gds
##

library(gdsfmt)
library(SNPRelate)
```

```
## SNPRelate -- supported by Streaming SIMD Extensions 2 (SSE2)
```

```r
# file name template
FILE_TEMPLATE <- "genotypes_chr%s_%s_phase3.3_consensus.b36_fwd.txt.gz"



#########################################################################
##  Download the gz files from HapMap Website

FTP_BASE <- "ftp://ftp.ncbi.nlm.nih.gov/hapmap/"
FTP_GENO_PATH <- "genotypes/2010-05_phaseIII/hapmap_format/consensus"
POP_LIST <- c("ASW", "CEU", "CHB", "CHD", "GIH", "JPT", "LWK", "MEX", "MKK", "TSI", "YRI")

for (chr.id in c(1:22, "X", "Y", "M"))
{
	for (pop in POP_LIST)
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
FTP_SAMP_PATH <- "genotypes/2010-05_phaseIII"
sample.fn <- "relationships_w_pops_041510.txt"
download.file(file.path(FTP_BASE, FTP_SAMP_PATH, sample.fn), sample.fn)

head(read.table(sample.fn, header=TRUE, stringsAsFactors=FALSE))
```

```
##    FID     IID     dad     mom sex pheno population
## 1 2427 NA19919 NA19908 NA19909   1     0        ASW
## 2 2431 NA19916       0       0   1     0        ASW
## 3 2424 NA19835       0       0   2     0        ASW
## 4 2469 NA20282       0       0   2     0        ASW
## 5 2368 NA19703       0       0   1     0        ASW
## 6 2425 NA19902 NA19900 NA19901   2     0        ASW
```

```r
#    FID     IID     dad     mom sex pheno population


#########################################################################
##  Get HapMap individual IDs

SampID_File <- function(chr.id)
{
	ans <- c()
	for (pop in POP_LIST)
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
```

```
## Chromosome 2: OK
## Chromosome 3: OK
## Chromosome 4: OK
## Chromosome 5: OK
## Chromosome 6: OK
## Chromosome 7: OK
## Chromosome 8: OK
## Chromosome 9: OK
## Chromosome 10: OK
## Chromosome 11: OK
## Chromosome 12: OK
## Chromosome 13: OK
## Chromosome 14: OK
## Chromosome 15: OK
## Chromosome 16: OK
## Chromosome 17: OK
## Chromosome 18: OK
## Chromosome 19: OK
## Chromosome 20: OK
## Chromosome 21: OK
## Chromosome 22: OK
## Chromosome X: OK
## Chromosome Y: OK
## Chromosome M: OK
```

```r
#########################################################################
##  Create GDS file

# Create a new GDS file
gds.fn <- "HapMap3_r3_b36_fwd_consensus_qc_poly.gds"
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
```

```
## Starting ...
```

```r
# Write SNPs into the GDS file
for (chr.id in c(1:22, "X", "Y"))
{
	# load data
	dat <- list()
	for (pop in POP_LIST)
	{
		fn <- sprintf(FILE_TEMPLATE, chr.id, pop)
		cat("Read:", fn)
		dat[[pop]] <- read.table(fn, header=TRUE, comment.char="",
			stringsAsFactors=FALSE)
		cat(",", dim(dat[[pop]]), "\n")

		# detect and remove non-standard alleles
		allele <- dat[[pop]]$alleles
		s <- strsplit(allele, "/")
		flag <- sapply(s, function(x)
			(length(x)==2L) & all(x %in% c("A", "G", "C", "T")))
		if (sum(!flag) > 0)
		{
			dat[[pop]] <- dat[[pop]][flag, ]
			allele <- dat[[pop]]$alleles
			cat(sprintf("\t%s: remove non-standard alleles (%d)\n",
				pop, sum(!flag)))
		}
	}

	# check position
	for (pop in POP_LIST)
	{
		if (any(duplicated(dat[[pop]]$pos)))
			stop("Positions in Population ", pop, " are not unique.")
	}
	common.pos <- dat[[POP_LIST[1]]]$pos
	for (pop in POP_LIST[-1])
		common.pos <- intersect(common.pos, dat[[pop]]$pos)
	common.pos <- common.pos[order(common.pos)]
	cat("\tPosition Intersect:", length(common.pos), "\n")
	for (pop in POP_LIST)
	{
		if (!identical(common.pos, dat[[pop]]$pos))
			dat[[pop]] <- dat[[pop]][match(common.pos, dat[[pop]]$pos), ]
	}

	# check RS ID
	for (pop in POP_LIST)
	{
		if (any(duplicated(dat[[pop]]$rs.)))
			stop("RS IDs in Population ", pop, " are not unique.")
	}
	common.rs. <- dat[[POP_LIST[1]]]$rs.
	for (pop in POP_LIST[-1])
		common.rs. <- intersect(common.rs., dat[[pop]]$rs.)
	cat("\tRS ID Intersect:", length(common.rs.), "\n")
	for (pop in POP_LIST)
	{
		if (!identical(common.rs., dat[[pop]]$rs.))
			dat[[pop]] <- dat[[pop]][match(common.rs., dat[[pop]]$rs.), ]
	}

	# check alleles
	allele <- dat[[POP_LIST[1]]]$alleles
	flag <- rep(TRUE, length(allele))
	for (pop in POP_LIST[-1])
	{
		flag <- flag & (allele == dat[[pop]]$alleles)
	}
	cat("\tAllele Intersect:", sum(flag), "\n")
	for (pop in POP_LIST)
		dat[[pop]] <- dat[[pop]][flag, ]


	# snp.id
	n <- dim(dat[[1]])[1]
	append.gdsn(var.snp, seq_len(n) + n.int)
	n.int <- n.int + n

	# snp.rs.id
	append.gdsn(var.rsid, dat[[1]]$rs.)

	# snp.chromosome
	append.gdsn(var.chr, rep(chr.id, n))

	# snp.position
	append.gdsn(var.pos, dat[[1]]$pos)

	# snp.allele
	allele <- dat[[1]]$alleles
	append.gdsn(var.allele, allele)

	s <- strsplit(allele, "/")
	g0  <- sapply(s, function(x) paste0(x[2],x[2]))
	g1a <- sapply(s, function(x) paste0(x[1],x[2]))
	g1b <- sapply(s, function(x) paste0(x[2],x[1]))
	g2  <- sapply(s, function(x) paste0(x[1],x[1]))
	gNN <- rep("NN", length(allele))

	# genotype
	m <- NULL
	for (pop in POP_LIST)
		m <- cbind(m, as.matrix(dat[[pop]][, -c(1:11)]))
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
```

```
## Read: genotypes_chr1_ASW_phase3.3_consensus.b36_fwd.txt.gz, 119487 98 
## 	ASW: remove non-standard alleles (58)
## Read: genotypes_chr1_CEU_phase3.3_consensus.b36_fwd.txt.gz, 119487 176 
## 	CEU: remove non-standard alleles (59)
## Read: genotypes_chr1_CHB_phase3.3_consensus.b36_fwd.txt.gz, 119487 148 
## 	CHB: remove non-standard alleles (62)
## Read: genotypes_chr1_CHD_phase3.3_consensus.b36_fwd.txt.gz, 119487 120 
## 	CHD: remove non-standard alleles (63)
## Read: genotypes_chr1_GIH_phase3.3_consensus.b36_fwd.txt.gz, 119487 112 
## 	GIH: remove non-standard alleles (61)
## Read: genotypes_chr1_JPT_phase3.3_consensus.b36_fwd.txt.gz, 119487 124 
## 	JPT: remove non-standard alleles (62)
## Read: genotypes_chr1_LWK_phase3.3_consensus.b36_fwd.txt.gz, 119487 121 
## 	LWK: remove non-standard alleles (62)
## Read: genotypes_chr1_MEX_phase3.3_consensus.b36_fwd.txt.gz, 119487 97 
## 	MEX: remove non-standard alleles (59)
## Read: genotypes_chr1_MKK_phase3.3_consensus.b36_fwd.txt.gz, 119487 195 
## 	MKK: remove non-standard alleles (57)
## Read: genotypes_chr1_TSI_phase3.3_consensus.b36_fwd.txt.gz, 119487 113 
## 	TSI: remove non-standard alleles (61)
## Read: genotypes_chr1_YRI_phase3.3_consensus.b36_fwd.txt.gz, 119487 214 
## 	YRI: remove non-standard alleles (61)
## 	Position Intersect: 119420 
## 	RS ID Intersect: 119420 
## 	Allele Intersect: 118856 
## 	Genotypes:
##        0        1        2        3 
## 60442882 44909950 60261024   427976 
## 
## Read: genotypes_chr2_ASW_phase3.3_consensus.b36_fwd.txt.gz, 119502 98 
## 	ASW: remove non-standard alleles (17)
## Read: genotypes_chr2_CEU_phase3.3_consensus.b36_fwd.txt.gz, 119502 176 
## 	CEU: remove non-standard alleles (17)
## Read: genotypes_chr2_CHB_phase3.3_consensus.b36_fwd.txt.gz, 119502 148 
## 	CHB: remove non-standard alleles (19)
## Read: genotypes_chr2_CHD_phase3.3_consensus.b36_fwd.txt.gz, 119502 120 
## 	CHD: remove non-standard alleles (19)
## Read: genotypes_chr2_GIH_phase3.3_consensus.b36_fwd.txt.gz, 119502 112 
## 	GIH: remove non-standard alleles (17)
## Read: genotypes_chr2_JPT_phase3.3_consensus.b36_fwd.txt.gz, 119502 124 
## 	JPT: remove non-standard alleles (19)
## Read: genotypes_chr2_LWK_phase3.3_consensus.b36_fwd.txt.gz, 119502 121 
## 	LWK: remove non-standard alleles (17)
## Read: genotypes_chr2_MEX_phase3.3_consensus.b36_fwd.txt.gz, 119502 97 
## 	MEX: remove non-standard alleles (17)
## Read: genotypes_chr2_MKK_phase3.3_consensus.b36_fwd.txt.gz, 119502 195 
## 	MKK: remove non-standard alleles (17)
## Read: genotypes_chr2_TSI_phase3.3_consensus.b36_fwd.txt.gz, 119502 113 
## 	TSI: remove non-standard alleles (17)
## Read: genotypes_chr2_YRI_phase3.3_consensus.b36_fwd.txt.gz, 119502 214 
## 	YRI: remove non-standard alleles (17)
## 	Position Intersect: 119481 
## 	RS ID Intersect: 119481 
## 	Allele Intersect: 119079 
## 	Genotypes:
##        0        1        2        3 
## 59892210 46082039 59931451   447663 
## 
## Read: genotypes_chr3_ASW_phase3.3_consensus.b36_fwd.txt.gz, 98971 98 
## 	ASW: remove non-standard alleles (12)
## Read: genotypes_chr3_CEU_phase3.3_consensus.b36_fwd.txt.gz, 98971 176 
## 	CEU: remove non-standard alleles (12)
## Read: genotypes_chr3_CHB_phase3.3_consensus.b36_fwd.txt.gz, 98971 148 
## 	CHB: remove non-standard alleles (13)
## Read: genotypes_chr3_CHD_phase3.3_consensus.b36_fwd.txt.gz, 98971 120 
## 	CHD: remove non-standard alleles (13)
## Read: genotypes_chr3_GIH_phase3.3_consensus.b36_fwd.txt.gz, 98971 112 
## 	GIH: remove non-standard alleles (12)
## Read: genotypes_chr3_JPT_phase3.3_consensus.b36_fwd.txt.gz, 98971 124 
## 	JPT: remove non-standard alleles (13)
## Read: genotypes_chr3_LWK_phase3.3_consensus.b36_fwd.txt.gz, 98971 121 
## 	LWK: remove non-standard alleles (13)
## Read: genotypes_chr3_MEX_phase3.3_consensus.b36_fwd.txt.gz, 98971 97 
## 	MEX: remove non-standard alleles (12)
## Read: genotypes_chr3_MKK_phase3.3_consensus.b36_fwd.txt.gz, 98971 195 
## 	MKK: remove non-standard alleles (12)
## Read: genotypes_chr3_TSI_phase3.3_consensus.b36_fwd.txt.gz, 98971 113 
## 	TSI: remove non-standard alleles (12)
## Read: genotypes_chr3_YRI_phase3.3_consensus.b36_fwd.txt.gz, 98971 214 
## 	YRI: remove non-standard alleles (13)
## 	Position Intersect: 98956 
## 	RS ID Intersect: 98956 
## 	Allele Intersect: 98618 
## 	Genotypes:
##        0        1        2        3 
## 49497197 38501693 49404923   365533 
## 
## Read: genotypes_chr4_ASW_phase3.3_consensus.b36_fwd.txt.gz, 88135 98 
## 	ASW: remove non-standard alleles (11)
## Read: genotypes_chr4_CEU_phase3.3_consensus.b36_fwd.txt.gz, 88135 176 
## 	CEU: remove non-standard alleles (12)
## Read: genotypes_chr4_CHB_phase3.3_consensus.b36_fwd.txt.gz, 88135 148 
## 	CHB: remove non-standard alleles (12)
## Read: genotypes_chr4_CHD_phase3.3_consensus.b36_fwd.txt.gz, 88135 120 
## 	CHD: remove non-standard alleles (12)
## Read: genotypes_chr4_GIH_phase3.3_consensus.b36_fwd.txt.gz, 88135 112 
## 	GIH: remove non-standard alleles (10)
## Read: genotypes_chr4_JPT_phase3.3_consensus.b36_fwd.txt.gz, 88135 124 
## 	JPT: remove non-standard alleles (12)
## Read: genotypes_chr4_LWK_phase3.3_consensus.b36_fwd.txt.gz, 88135 121 
## 	LWK: remove non-standard alleles (12)
## Read: genotypes_chr4_MEX_phase3.3_consensus.b36_fwd.txt.gz, 88135 97 
## 	MEX: remove non-standard alleles (9)
## Read: genotypes_chr4_MKK_phase3.3_consensus.b36_fwd.txt.gz, 88135 195 
## 	MKK: remove non-standard alleles (12)
## Read: genotypes_chr4_TSI_phase3.3_consensus.b36_fwd.txt.gz, 88135 113 
## 	TSI: remove non-standard alleles (11)
## Read: genotypes_chr4_YRI_phase3.3_consensus.b36_fwd.txt.gz, 88135 214 
## 	YRI: remove non-standard alleles (13)
## 	Position Intersect: 88120 
## 	RS ID Intersect: 88120 
## 	Allele Intersect: 87878 
## 	Genotypes:
##        0        1        2        3 
## 43959495 34287124 44173370   345577 
## 
## Read: genotypes_chr5_ASW_phase3.3_consensus.b36_fwd.txt.gz, 90368 98 
## 	ASW: remove non-standard alleles (15)
## Read: genotypes_chr5_CEU_phase3.3_consensus.b36_fwd.txt.gz, 90368 176 
## 	CEU: remove non-standard alleles (16)
## Read: genotypes_chr5_CHB_phase3.3_consensus.b36_fwd.txt.gz, 90368 148 
## 	CHB: remove non-standard alleles (17)
## Read: genotypes_chr5_CHD_phase3.3_consensus.b36_fwd.txt.gz, 90368 120 
## 	CHD: remove non-standard alleles (22)
## Read: genotypes_chr5_GIH_phase3.3_consensus.b36_fwd.txt.gz, 90368 112 
## 	GIH: remove non-standard alleles (16)
## Read: genotypes_chr5_JPT_phase3.3_consensus.b36_fwd.txt.gz, 90368 124 
## 	JPT: remove non-standard alleles (22)
## Read: genotypes_chr5_LWK_phase3.3_consensus.b36_fwd.txt.gz, 90368 121 
## 	LWK: remove non-standard alleles (16)
## Read: genotypes_chr5_MEX_phase3.3_consensus.b36_fwd.txt.gz, 90368 97 
## 	MEX: remove non-standard alleles (15)
## Read: genotypes_chr5_MKK_phase3.3_consensus.b36_fwd.txt.gz, 90368 195 
## 	MKK: remove non-standard alleles (14)
## Read: genotypes_chr5_TSI_phase3.3_consensus.b36_fwd.txt.gz, 90368 113 
## 	TSI: remove non-standard alleles (17)
## Read: genotypes_chr5_YRI_phase3.3_consensus.b36_fwd.txt.gz, 90368 214 
## 	YRI: remove non-standard alleles (20)
## 	Position Intersect: 90340 
## 	RS ID Intersect: 90340 
## 	Allele Intersect: 90054 
## 	Genotypes:
##        0        1        2        3 
## 45226484 35223436 45021267   334251 
## 
## Read: genotypes_chr6_ASW_phase3.3_consensus.b36_fwd.txt.gz, 93671 98 
## 	ASW: remove non-standard alleles (38)
## Read: genotypes_chr6_CEU_phase3.3_consensus.b36_fwd.txt.gz, 93671 176 
## 	CEU: remove non-standard alleles (38)
## Read: genotypes_chr6_CHB_phase3.3_consensus.b36_fwd.txt.gz, 93671 148 
## 	CHB: remove non-standard alleles (37)
## Read: genotypes_chr6_CHD_phase3.3_consensus.b36_fwd.txt.gz, 93671 120 
## 	CHD: remove non-standard alleles (38)
## Read: genotypes_chr6_GIH_phase3.3_consensus.b36_fwd.txt.gz, 93671 112 
## 	GIH: remove non-standard alleles (38)
## Read: genotypes_chr6_JPT_phase3.3_consensus.b36_fwd.txt.gz, 93671 124 
## 	JPT: remove non-standard alleles (38)
## Read: genotypes_chr6_LWK_phase3.3_consensus.b36_fwd.txt.gz, 93671 121 
## 	LWK: remove non-standard alleles (39)
## Read: genotypes_chr6_MEX_phase3.3_consensus.b36_fwd.txt.gz, 93671 97 
## 	MEX: remove non-standard alleles (38)
## Read: genotypes_chr6_MKK_phase3.3_consensus.b36_fwd.txt.gz, 93671 195 
## 	MKK: remove non-standard alleles (37)
## Read: genotypes_chr6_TSI_phase3.3_consensus.b36_fwd.txt.gz, 93671 113 
## 	TSI: remove non-standard alleles (38)
## Read: genotypes_chr6_YRI_phase3.3_consensus.b36_fwd.txt.gz, 93671 214 
## 	YRI: remove non-standard alleles (40)
## 	Position Intersect: 93631 
## 	RS ID Intersect: 93631 
## 	Allele Intersect: 93325 
## 	Genotypes:
##        0        1        2        3 
## 46891470 36115726 47019442   348387 
## 
## Read: genotypes_chr7_ASW_phase3.3_consensus.b36_fwd.txt.gz, 77377 98 
## 	ASW: remove non-standard alleles (21)
## Read: genotypes_chr7_CEU_phase3.3_consensus.b36_fwd.txt.gz, 77377 176 
## 	CEU: remove non-standard alleles (22)
## Read: genotypes_chr7_CHB_phase3.3_consensus.b36_fwd.txt.gz, 77377 148 
## 	CHB: remove non-standard alleles (24)
## Read: genotypes_chr7_CHD_phase3.3_consensus.b36_fwd.txt.gz, 77377 120 
## 	CHD: remove non-standard alleles (25)
## Read: genotypes_chr7_GIH_phase3.3_consensus.b36_fwd.txt.gz, 77377 112 
## 	GIH: remove non-standard alleles (22)
## Read: genotypes_chr7_JPT_phase3.3_consensus.b36_fwd.txt.gz, 77377 124 
## 	JPT: remove non-standard alleles (24)
## Read: genotypes_chr7_LWK_phase3.3_consensus.b36_fwd.txt.gz, 77377 121 
## 	LWK: remove non-standard alleles (21)
## Read: genotypes_chr7_MEX_phase3.3_consensus.b36_fwd.txt.gz, 77377 97 
## 	MEX: remove non-standard alleles (22)
## Read: genotypes_chr7_MKK_phase3.3_consensus.b36_fwd.txt.gz, 77377 195 
## 	MKK: remove non-standard alleles (21)
## Read: genotypes_chr7_TSI_phase3.3_consensus.b36_fwd.txt.gz, 77377 113 
## 	TSI: remove non-standard alleles (22)
## Read: genotypes_chr7_YRI_phase3.3_consensus.b36_fwd.txt.gz, 77377 214 
## 	YRI: remove non-standard alleles (22)
## 	Position Intersect: 77350 
## 	RS ID Intersect: 77350 
## 	Allele Intersect: 77075 
## 	Genotypes:
##        0        1        2        3 
## 38588384 30211729 38577810   295852 
## 
## Read: genotypes_chr8_ASW_phase3.3_consensus.b36_fwd.txt.gz, 77111 98 
## 	ASW: remove non-standard alleles (11)
## Read: genotypes_chr8_CEU_phase3.3_consensus.b36_fwd.txt.gz, 77111 176 
## 	CEU: remove non-standard alleles (11)
## Read: genotypes_chr8_CHB_phase3.3_consensus.b36_fwd.txt.gz, 77111 148 
## 	CHB: remove non-standard alleles (11)
## Read: genotypes_chr8_CHD_phase3.3_consensus.b36_fwd.txt.gz, 77111 120 
## 	CHD: remove non-standard alleles (11)
## Read: genotypes_chr8_GIH_phase3.3_consensus.b36_fwd.txt.gz, 77111 112 
## 	GIH: remove non-standard alleles (10)
## Read: genotypes_chr8_JPT_phase3.3_consensus.b36_fwd.txt.gz, 77111 124 
## 	JPT: remove non-standard alleles (11)
## Read: genotypes_chr8_LWK_phase3.3_consensus.b36_fwd.txt.gz, 77111 121 
## 	LWK: remove non-standard alleles (10)
## Read: genotypes_chr8_MEX_phase3.3_consensus.b36_fwd.txt.gz, 77111 97 
## 	MEX: remove non-standard alleles (11)
## Read: genotypes_chr8_MKK_phase3.3_consensus.b36_fwd.txt.gz, 77111 195 
## 	MKK: remove non-standard alleles (11)
## Read: genotypes_chr8_TSI_phase3.3_consensus.b36_fwd.txt.gz, 77111 113 
## 	TSI: remove non-standard alleles (11)
## Read: genotypes_chr8_YRI_phase3.3_consensus.b36_fwd.txt.gz, 77111 214 
## 	YRI: remove non-standard alleles (12)
## 	Position Intersect: 77098 
## 	RS ID Intersect: 77098 
## 	Allele Intersect: 76870 
## 	Genotypes:
##        0        1        2        3 
## 38548510 30123495 38420492   294893 
## 
## Read: genotypes_chr9_ASW_phase3.3_consensus.b36_fwd.txt.gz, 65251 98 
## 	ASW: remove non-standard alleles (11)
## Read: genotypes_chr9_CEU_phase3.3_consensus.b36_fwd.txt.gz, 65251 176 
## 	CEU: remove non-standard alleles (12)
## Read: genotypes_chr9_CHB_phase3.3_consensus.b36_fwd.txt.gz, 65251 148 
## 	CHB: remove non-standard alleles (12)
## Read: genotypes_chr9_CHD_phase3.3_consensus.b36_fwd.txt.gz, 65251 120 
## 	CHD: remove non-standard alleles (13)
## Read: genotypes_chr9_GIH_phase3.3_consensus.b36_fwd.txt.gz, 65251 112 
## 	GIH: remove non-standard alleles (13)
## Read: genotypes_chr9_JPT_phase3.3_consensus.b36_fwd.txt.gz, 65251 124 
## 	JPT: remove non-standard alleles (12)
## Read: genotypes_chr9_LWK_phase3.3_consensus.b36_fwd.txt.gz, 65251 121 
## 	LWK: remove non-standard alleles (10)
## Read: genotypes_chr9_MEX_phase3.3_consensus.b36_fwd.txt.gz, 65251 97 
## 	MEX: remove non-standard alleles (11)
## Read: genotypes_chr9_MKK_phase3.3_consensus.b36_fwd.txt.gz, 65251 195 
## 	MKK: remove non-standard alleles (11)
## Read: genotypes_chr9_TSI_phase3.3_consensus.b36_fwd.txt.gz, 65251 113 
## 	TSI: remove non-standard alleles (11)
## Read: genotypes_chr9_YRI_phase3.3_consensus.b36_fwd.txt.gz, 65251 214 
## 	YRI: remove non-standard alleles (10)
## 	Position Intersect: 65238 
## 	RS ID Intersect: 65238 
## 	Allele Intersect: 65050 
## 	Genotypes:
##        0        1        2        3 
## 32651941 25429693 32543150   250066 
## 
## Read: genotypes_chr10_ASW_phase3.3_consensus.b36_fwd.txt.gz, 75616 98 
## 	ASW: remove non-standard alleles (16)
## Read: genotypes_chr10_CEU_phase3.3_consensus.b36_fwd.txt.gz, 75616 176 
## 	CEU: remove non-standard alleles (15)
## Read: genotypes_chr10_CHB_phase3.3_consensus.b36_fwd.txt.gz, 75616 148 
## 	CHB: remove non-standard alleles (15)
## Read: genotypes_chr10_CHD_phase3.3_consensus.b36_fwd.txt.gz, 75616 120 
## 	CHD: remove non-standard alleles (15)
## Read: genotypes_chr10_GIH_phase3.3_consensus.b36_fwd.txt.gz, 75616 112 
## 	GIH: remove non-standard alleles (15)
## Read: genotypes_chr10_JPT_phase3.3_consensus.b36_fwd.txt.gz, 75616 124 
## 	JPT: remove non-standard alleles (15)
## Read: genotypes_chr10_LWK_phase3.3_consensus.b36_fwd.txt.gz, 75616 121 
## 	LWK: remove non-standard alleles (14)
## Read: genotypes_chr10_MEX_phase3.3_consensus.b36_fwd.txt.gz, 75616 97 
## 	MEX: remove non-standard alleles (15)
## Read: genotypes_chr10_MKK_phase3.3_consensus.b36_fwd.txt.gz, 75616 195 
## 	MKK: remove non-standard alleles (14)
## Read: genotypes_chr10_TSI_phase3.3_consensus.b36_fwd.txt.gz, 75616 113 
## 	TSI: remove non-standard alleles (14)
## Read: genotypes_chr10_YRI_phase3.3_consensus.b36_fwd.txt.gz, 75616 214 
## 	YRI: remove non-standard alleles (16)
## 	Position Intersect: 75599 
## 	RS ID Intersect: 75599 
## 	Allele Intersect: 75363 
## 	Genotypes:
##        0        1        2        3 
## 37962984 28974637 38075204   269286 
## 
## Read: genotypes_chr11_ASW_phase3.3_consensus.b36_fwd.txt.gz, 72993 98 
## 	ASW: remove non-standard alleles (21)
## Read: genotypes_chr11_CEU_phase3.3_consensus.b36_fwd.txt.gz, 72993 176 
## 	CEU: remove non-standard alleles (21)
## Read: genotypes_chr11_CHB_phase3.3_consensus.b36_fwd.txt.gz, 72993 148 
## 	CHB: remove non-standard alleles (22)
## Read: genotypes_chr11_CHD_phase3.3_consensus.b36_fwd.txt.gz, 72993 120 
## 	CHD: remove non-standard alleles (25)
## Read: genotypes_chr11_GIH_phase3.3_consensus.b36_fwd.txt.gz, 72993 112 
## 	GIH: remove non-standard alleles (22)
## Read: genotypes_chr11_JPT_phase3.3_consensus.b36_fwd.txt.gz, 72993 124 
## 	JPT: remove non-standard alleles (23)
## Read: genotypes_chr11_LWK_phase3.3_consensus.b36_fwd.txt.gz, 72993 121 
## 	LWK: remove non-standard alleles (21)
## Read: genotypes_chr11_MEX_phase3.3_consensus.b36_fwd.txt.gz, 72993 97 
## 	MEX: remove non-standard alleles (21)
## Read: genotypes_chr11_MKK_phase3.3_consensus.b36_fwd.txt.gz, 72993 195 
## 	MKK: remove non-standard alleles (21)
## Read: genotypes_chr11_TSI_phase3.3_consensus.b36_fwd.txt.gz, 72993 113 
## 	TSI: remove non-standard alleles (22)
## Read: genotypes_chr11_YRI_phase3.3_consensus.b36_fwd.txt.gz, 72993 214 
## 	YRI: remove non-standard alleles (20)
## 	Position Intersect: 72967 
## 	RS ID Intersect: 72967 
## 	Allele Intersect: 72682 
## 	Genotypes:
##        0        1        2        3 
## 36569339 28230781 36463872   272762 
## 
## Read: genotypes_chr12_ASW_phase3.3_consensus.b36_fwd.txt.gz, 70482 98 
## 	ASW: remove non-standard alleles (32)
## Read: genotypes_chr12_CEU_phase3.3_consensus.b36_fwd.txt.gz, 70482 176 
## 	CEU: remove non-standard alleles (33)
## Read: genotypes_chr12_CHB_phase3.3_consensus.b36_fwd.txt.gz, 70482 148 
## 	CHB: remove non-standard alleles (33)
## Read: genotypes_chr12_CHD_phase3.3_consensus.b36_fwd.txt.gz, 70482 120 
## 	CHD: remove non-standard alleles (34)
## Read: genotypes_chr12_GIH_phase3.3_consensus.b36_fwd.txt.gz, 70482 112 
## 	GIH: remove non-standard alleles (33)
## Read: genotypes_chr12_JPT_phase3.3_consensus.b36_fwd.txt.gz, 70482 124 
## 	JPT: remove non-standard alleles (33)
## Read: genotypes_chr12_LWK_phase3.3_consensus.b36_fwd.txt.gz, 70482 121 
## 	LWK: remove non-standard alleles (31)
## Read: genotypes_chr12_MEX_phase3.3_consensus.b36_fwd.txt.gz, 70482 97 
## 	MEX: remove non-standard alleles (32)
## Read: genotypes_chr12_MKK_phase3.3_consensus.b36_fwd.txt.gz, 70482 195 
## 	MKK: remove non-standard alleles (32)
## Read: genotypes_chr12_TSI_phase3.3_consensus.b36_fwd.txt.gz, 70482 113 
## 	TSI: remove non-standard alleles (35)
## Read: genotypes_chr12_YRI_phase3.3_consensus.b36_fwd.txt.gz, 70482 214 
## 	YRI: remove non-standard alleles (31)
## 	Position Intersect: 70447 
## 	RS ID Intersect: 70447 
## 	Allele Intersect: 70182 
## 	Genotypes:
##        0        1        2        3 
## 35482435 26709344 35534641   317834 
## 
## Read: genotypes_chr13_ASW_phase3.3_consensus.b36_fwd.txt.gz, 53293 98 
## 	ASW: remove non-standard alleles (8)
## Read: genotypes_chr13_CEU_phase3.3_consensus.b36_fwd.txt.gz, 53293 176 
## 	CEU: remove non-standard alleles (9)
## Read: genotypes_chr13_CHB_phase3.3_consensus.b36_fwd.txt.gz, 53293 148 
## 	CHB: remove non-standard alleles (10)
## Read: genotypes_chr13_CHD_phase3.3_consensus.b36_fwd.txt.gz, 53293 120 
## 	CHD: remove non-standard alleles (10)
## Read: genotypes_chr13_GIH_phase3.3_consensus.b36_fwd.txt.gz, 53293 112 
## 	GIH: remove non-standard alleles (9)
## Read: genotypes_chr13_JPT_phase3.3_consensus.b36_fwd.txt.gz, 53293 124 
## 	JPT: remove non-standard alleles (10)
## Read: genotypes_chr13_LWK_phase3.3_consensus.b36_fwd.txt.gz, 53293 121 
## 	LWK: remove non-standard alleles (8)
## Read: genotypes_chr13_MEX_phase3.3_consensus.b36_fwd.txt.gz, 53293 97 
## 	MEX: remove non-standard alleles (9)
## Read: genotypes_chr13_MKK_phase3.3_consensus.b36_fwd.txt.gz, 53293 195 
## 	MKK: remove non-standard alleles (9)
## Read: genotypes_chr13_TSI_phase3.3_consensus.b36_fwd.txt.gz, 53293 113 
## 	TSI: remove non-standard alleles (9)
## Read: genotypes_chr13_YRI_phase3.3_consensus.b36_fwd.txt.gz, 53293 214 
## 	YRI: remove non-standard alleles (9)
## 	Position Intersect: 53282 
## 	RS ID Intersect: 53282 
## 	Allele Intersect: 53163 
## 	Genotypes:
##        0        1        2        3 
## 26696920 20645349 26720531   205911 
## 
## Read: genotypes_chr14_ASW_phase3.3_consensus.b36_fwd.txt.gz, 46655 98 
## 	ASW: remove non-standard alleles (10)
## Read: genotypes_chr14_CEU_phase3.3_consensus.b36_fwd.txt.gz, 46655 176 
## 	CEU: remove non-standard alleles (10)
## Read: genotypes_chr14_CHB_phase3.3_consensus.b36_fwd.txt.gz, 46655 148 
## 	CHB: remove non-standard alleles (10)
## Read: genotypes_chr14_CHD_phase3.3_consensus.b36_fwd.txt.gz, 46655 120 
## 	CHD: remove non-standard alleles (10)
## Read: genotypes_chr14_GIH_phase3.3_consensus.b36_fwd.txt.gz, 46655 112 
## 	GIH: remove non-standard alleles (10)
## Read: genotypes_chr14_JPT_phase3.3_consensus.b36_fwd.txt.gz, 46655 124 
## 	JPT: remove non-standard alleles (9)
## Read: genotypes_chr14_LWK_phase3.3_consensus.b36_fwd.txt.gz, 46655 121 
## 	LWK: remove non-standard alleles (10)
## Read: genotypes_chr14_MEX_phase3.3_consensus.b36_fwd.txt.gz, 46655 97 
## 	MEX: remove non-standard alleles (11)
## Read: genotypes_chr14_MKK_phase3.3_consensus.b36_fwd.txt.gz, 46655 195 
## 	MKK: remove non-standard alleles (10)
## Read: genotypes_chr14_TSI_phase3.3_consensus.b36_fwd.txt.gz, 46655 113 
## 	TSI: remove non-standard alleles (10)
## Read: genotypes_chr14_YRI_phase3.3_consensus.b36_fwd.txt.gz, 46655 214 
## 	YRI: remove non-standard alleles (11)
## 	Position Intersect: 46644 
## 	RS ID Intersect: 46644 
## 	Allele Intersect: 46505 
## 	Genotypes:
##        0        1        2        3 
## 23406373 18033206 23348298   179608 
## 
## Read: genotypes_chr15_ASW_phase3.3_consensus.b36_fwd.txt.gz, 43309 98 
## 	ASW: remove non-standard alleles (14)
## Read: genotypes_chr15_CEU_phase3.3_consensus.b36_fwd.txt.gz, 43309 176 
## 	CEU: remove non-standard alleles (14)
## Read: genotypes_chr15_CHB_phase3.3_consensus.b36_fwd.txt.gz, 43309 148 
## 	CHB: remove non-standard alleles (14)
## Read: genotypes_chr15_CHD_phase3.3_consensus.b36_fwd.txt.gz, 43309 120 
## 	CHD: remove non-standard alleles (14)
## Read: genotypes_chr15_GIH_phase3.3_consensus.b36_fwd.txt.gz, 43309 112 
## 	GIH: remove non-standard alleles (14)
## Read: genotypes_chr15_JPT_phase3.3_consensus.b36_fwd.txt.gz, 43309 124 
## 	JPT: remove non-standard alleles (14)
## Read: genotypes_chr15_LWK_phase3.3_consensus.b36_fwd.txt.gz, 43309 121 
## 	LWK: remove non-standard alleles (15)
## Read: genotypes_chr15_MEX_phase3.3_consensus.b36_fwd.txt.gz, 43309 97 
## 	MEX: remove non-standard alleles (14)
## Read: genotypes_chr15_MKK_phase3.3_consensus.b36_fwd.txt.gz, 43309 195 
## 	MKK: remove non-standard alleles (15)
## Read: genotypes_chr15_TSI_phase3.3_consensus.b36_fwd.txt.gz, 43309 113 
## 	TSI: remove non-standard alleles (14)
## Read: genotypes_chr15_YRI_phase3.3_consensus.b36_fwd.txt.gz, 43309 214 
## 	YRI: remove non-standard alleles (16)
## 	Position Intersect: 43293 
## 	RS ID Intersect: 43293 
## 	Allele Intersect: 43132 
## 	Genotypes:
##        0        1        2        3 
## 21626823 16832371 21630291   165919 
## 
## Read: genotypes_chr16_ASW_phase3.3_consensus.b36_fwd.txt.gz, 45778 98 
## 	ASW: remove non-standard alleles (12)
## Read: genotypes_chr16_CEU_phase3.3_consensus.b36_fwd.txt.gz, 45778 176 
## 	CEU: remove non-standard alleles (16)
## Read: genotypes_chr16_CHB_phase3.3_consensus.b36_fwd.txt.gz, 45778 148 
## 	CHB: remove non-standard alleles (16)
## Read: genotypes_chr16_CHD_phase3.3_consensus.b36_fwd.txt.gz, 45778 120 
## 	CHD: remove non-standard alleles (16)
## Read: genotypes_chr16_GIH_phase3.3_consensus.b36_fwd.txt.gz, 45778 112 
## 	GIH: remove non-standard alleles (17)
## Read: genotypes_chr16_JPT_phase3.3_consensus.b36_fwd.txt.gz, 45778 124 
## 	JPT: remove non-standard alleles (16)
## Read: genotypes_chr16_LWK_phase3.3_consensus.b36_fwd.txt.gz, 45778 121 
## 	LWK: remove non-standard alleles (12)
## Read: genotypes_chr16_MEX_phase3.3_consensus.b36_fwd.txt.gz, 45778 97 
## 	MEX: remove non-standard alleles (15)
## Read: genotypes_chr16_MKK_phase3.3_consensus.b36_fwd.txt.gz, 45778 195 
## 	MKK: remove non-standard alleles (15)
## Read: genotypes_chr16_TSI_phase3.3_consensus.b36_fwd.txt.gz, 45778 113 
## 	TSI: remove non-standard alleles (14)
## Read: genotypes_chr16_YRI_phase3.3_consensus.b36_fwd.txt.gz, 45778 214 
## 	YRI: remove non-standard alleles (16)
## 	Position Intersect: 45757 
## 	RS ID Intersect: 45757 
## 	Allele Intersect: 45536 
## 	Genotypes:
##        0        1        2        3 
## 22818095 17451055 23191387   153255 
## 
## Read: genotypes_chr17_ASW_phase3.3_consensus.b36_fwd.txt.gz, 39329 98 
## 	ASW: remove non-standard alleles (17)
## Read: genotypes_chr17_CEU_phase3.3_consensus.b36_fwd.txt.gz, 39329 176 
## 	CEU: remove non-standard alleles (16)
## Read: genotypes_chr17_CHB_phase3.3_consensus.b36_fwd.txt.gz, 39329 148 
## 	CHB: remove non-standard alleles (17)
## Read: genotypes_chr17_CHD_phase3.3_consensus.b36_fwd.txt.gz, 39329 120 
## 	CHD: remove non-standard alleles (17)
## Read: genotypes_chr17_GIH_phase3.3_consensus.b36_fwd.txt.gz, 39329 112 
## 	GIH: remove non-standard alleles (17)
## Read: genotypes_chr17_JPT_phase3.3_consensus.b36_fwd.txt.gz, 39329 124 
## 	JPT: remove non-standard alleles (17)
## Read: genotypes_chr17_LWK_phase3.3_consensus.b36_fwd.txt.gz, 39329 121 
## 	LWK: remove non-standard alleles (17)
## Read: genotypes_chr17_MEX_phase3.3_consensus.b36_fwd.txt.gz, 39329 97 
## 	MEX: remove non-standard alleles (17)
## Read: genotypes_chr17_MKK_phase3.3_consensus.b36_fwd.txt.gz, 39329 195 
## 	MKK: remove non-standard alleles (17)
## Read: genotypes_chr17_TSI_phase3.3_consensus.b36_fwd.txt.gz, 39329 113 
## 	TSI: remove non-standard alleles (16)
## Read: genotypes_chr17_YRI_phase3.3_consensus.b36_fwd.txt.gz, 39329 214 
## 	YRI: remove non-standard alleles (18)
## 	Position Intersect: 39311 
## 	RS ID Intersect: 39311 
## 	Allele Intersect: 39041 
## 	Genotypes:
##        0        1        2        3 
## 19821704 14983086 19602627   132860 
## 
## Read: genotypes_chr18_ASW_phase3.3_consensus.b36_fwd.txt.gz, 41942 98 
## 	ASW: remove non-standard alleles (2)
## Read: genotypes_chr18_CEU_phase3.3_consensus.b36_fwd.txt.gz, 41942 176 
## 	CEU: remove non-standard alleles (2)
## Read: genotypes_chr18_CHB_phase3.3_consensus.b36_fwd.txt.gz, 41942 148 
## 	CHB: remove non-standard alleles (1)
## Read: genotypes_chr18_CHD_phase3.3_consensus.b36_fwd.txt.gz, 41942 120 
## 	CHD: remove non-standard alleles (2)
## Read: genotypes_chr18_GIH_phase3.3_consensus.b36_fwd.txt.gz, 41942 112 
## 	GIH: remove non-standard alleles (2)
## Read: genotypes_chr18_JPT_phase3.3_consensus.b36_fwd.txt.gz, 41942 124 
## 	JPT: remove non-standard alleles (2)
## Read: genotypes_chr18_LWK_phase3.3_consensus.b36_fwd.txt.gz, 41942 121 
## 	LWK: remove non-standard alleles (1)
## Read: genotypes_chr18_MEX_phase3.3_consensus.b36_fwd.txt.gz, 41942 97 
## 	MEX: remove non-standard alleles (1)
## Read: genotypes_chr18_MKK_phase3.3_consensus.b36_fwd.txt.gz, 41942 195 
## 	MKK: remove non-standard alleles (2)
## Read: genotypes_chr18_TSI_phase3.3_consensus.b36_fwd.txt.gz, 41942 113 
## 	TSI: remove non-standard alleles (2)
## Read: genotypes_chr18_YRI_phase3.3_consensus.b36_fwd.txt.gz, 41942 214 
## 	YRI: remove non-standard alleles (2)
## 	Position Intersect: 41940 
## 	RS ID Intersect: 41940 
## 	Allele Intersect: 41832 
## 	Genotypes:
##        0        1        2        3 
## 20806034 16418478 21056093   158699 
## 
## Read: genotypes_chr19_ASW_phase3.3_consensus.b36_fwd.txt.gz, 26953 98 
## 	ASW: remove non-standard alleles (19)
## Read: genotypes_chr19_CEU_phase3.3_consensus.b36_fwd.txt.gz, 26953 176 
## 	CEU: remove non-standard alleles (21)
## Read: genotypes_chr19_CHB_phase3.3_consensus.b36_fwd.txt.gz, 26953 148 
## 	CHB: remove non-standard alleles (23)
## Read: genotypes_chr19_CHD_phase3.3_consensus.b36_fwd.txt.gz, 26953 120 
## 	CHD: remove non-standard alleles (24)
## Read: genotypes_chr19_GIH_phase3.3_consensus.b36_fwd.txt.gz, 26953 112 
## 	GIH: remove non-standard alleles (21)
## Read: genotypes_chr19_JPT_phase3.3_consensus.b36_fwd.txt.gz, 26953 124 
## 	JPT: remove non-standard alleles (21)
## Read: genotypes_chr19_LWK_phase3.3_consensus.b36_fwd.txt.gz, 26953 121 
## 	LWK: remove non-standard alleles (21)
## Read: genotypes_chr19_MEX_phase3.3_consensus.b36_fwd.txt.gz, 26953 97 
## 	MEX: remove non-standard alleles (19)
## Read: genotypes_chr19_MKK_phase3.3_consensus.b36_fwd.txt.gz, 26953 195 
## 	MKK: remove non-standard alleles (19)
## Read: genotypes_chr19_TSI_phase3.3_consensus.b36_fwd.txt.gz, 26953 113 
## 	TSI: remove non-standard alleles (20)
## Read: genotypes_chr19_YRI_phase3.3_consensus.b36_fwd.txt.gz, 26953 214 
## 	YRI: remove non-standard alleles (21)
## 	Position Intersect: 26929 
## 	RS ID Intersect: 26929 
## 	Allele Intersect: 26713 
## 	Genotypes:
##        0        1        2        3 
## 13462585 10268154 13488174    99148 
## 
## Read: genotypes_chr20_ASW_phase3.3_consensus.b36_fwd.txt.gz, 37159 98 
## 	ASW: remove non-standard alleles (2)
## Read: genotypes_chr20_CEU_phase3.3_consensus.b36_fwd.txt.gz, 37159 176 
## 	CEU: remove non-standard alleles (2)
## Read: genotypes_chr20_CHB_phase3.3_consensus.b36_fwd.txt.gz, 37159 148 
## 	CHB: remove non-standard alleles (2)
## Read: genotypes_chr20_CHD_phase3.3_consensus.b36_fwd.txt.gz, 37159 120 
## 	CHD: remove non-standard alleles (2)
## Read: genotypes_chr20_GIH_phase3.3_consensus.b36_fwd.txt.gz, 37159 112 
## 	GIH: remove non-standard alleles (2)
## Read: genotypes_chr20_JPT_phase3.3_consensus.b36_fwd.txt.gz, 37159 124 
## 	JPT: remove non-standard alleles (2)
## Read: genotypes_chr20_LWK_phase3.3_consensus.b36_fwd.txt.gz, 37159 121 
## 	LWK: remove non-standard alleles (2)
## Read: genotypes_chr20_MEX_phase3.3_consensus.b36_fwd.txt.gz, 37159 97 
## 	MEX: remove non-standard alleles (2)
## Read: genotypes_chr20_MKK_phase3.3_consensus.b36_fwd.txt.gz, 37159 195 
## 	MKK: remove non-standard alleles (2)
## Read: genotypes_chr20_TSI_phase3.3_consensus.b36_fwd.txt.gz, 37159 113 
## 	TSI: remove non-standard alleles (2)
## Read: genotypes_chr20_YRI_phase3.3_consensus.b36_fwd.txt.gz, 37159 214 
## 	YRI: remove non-standard alleles (2)
## 	Position Intersect: 37157 
## 	RS ID Intersect: 37157 
## 	Allele Intersect: 37009 
## 	Genotypes:
##        0        1        2        3 
## 18511454 14330205 18730023   129891 
## 
## Read: genotypes_chr21_ASW_phase3.3_consensus.b36_fwd.txt.gz, 19802 98 
## 	ASW: remove non-standard alleles (3)
## Read: genotypes_chr21_CEU_phase3.3_consensus.b36_fwd.txt.gz, 19802 176 
## 	CEU: remove non-standard alleles (3)
## Read: genotypes_chr21_CHB_phase3.3_consensus.b36_fwd.txt.gz, 19802 148 
## 	CHB: remove non-standard alleles (3)
## Read: genotypes_chr21_CHD_phase3.3_consensus.b36_fwd.txt.gz, 19802 120 
## 	CHD: remove non-standard alleles (3)
## Read: genotypes_chr21_GIH_phase3.3_consensus.b36_fwd.txt.gz, 19802 112 
## 	GIH: remove non-standard alleles (4)
## Read: genotypes_chr21_JPT_phase3.3_consensus.b36_fwd.txt.gz, 19802 124 
## 	JPT: remove non-standard alleles (4)
## Read: genotypes_chr21_LWK_phase3.3_consensus.b36_fwd.txt.gz, 19802 121 
## 	LWK: remove non-standard alleles (3)
## Read: genotypes_chr21_MEX_phase3.3_consensus.b36_fwd.txt.gz, 19802 97 
## 	MEX: remove non-standard alleles (3)
## Read: genotypes_chr21_MKK_phase3.3_consensus.b36_fwd.txt.gz, 19802 195 
## 	MKK: remove non-standard alleles (3)
## Read: genotypes_chr21_TSI_phase3.3_consensus.b36_fwd.txt.gz, 19802 113 
## 	TSI: remove non-standard alleles (3)
## Read: genotypes_chr21_YRI_phase3.3_consensus.b36_fwd.txt.gz, 19802 214 
## 	YRI: remove non-standard alleles (3)
## 	Position Intersect: 19798 
## 	RS ID Intersect: 19798 
## 	Allele Intersect: 19756 
## 	Genotypes:
##       0       1       2       3 
## 9828188 7909846 9783068   78030 
## 
## Read: genotypes_chr22_ASW_phase3.3_consensus.b36_fwd.txt.gz, 20649 98 
## 	ASW: remove non-standard alleles (13)
## Read: genotypes_chr22_CEU_phase3.3_consensus.b36_fwd.txt.gz, 20649 176 
## 	CEU: remove non-standard alleles (16)
## Read: genotypes_chr22_CHB_phase3.3_consensus.b36_fwd.txt.gz, 20649 148 
## 	CHB: remove non-standard alleles (15)
## Read: genotypes_chr22_CHD_phase3.3_consensus.b36_fwd.txt.gz, 20649 120 
## 	CHD: remove non-standard alleles (15)
## Read: genotypes_chr22_GIH_phase3.3_consensus.b36_fwd.txt.gz, 20649 112 
## 	GIH: remove non-standard alleles (15)
## Read: genotypes_chr22_JPT_phase3.3_consensus.b36_fwd.txt.gz, 20649 124 
## 	JPT: remove non-standard alleles (15)
## Read: genotypes_chr22_LWK_phase3.3_consensus.b36_fwd.txt.gz, 20649 121 
## 	LWK: remove non-standard alleles (14)
## Read: genotypes_chr22_MEX_phase3.3_consensus.b36_fwd.txt.gz, 20649 97 
## 	MEX: remove non-standard alleles (15)
## Read: genotypes_chr22_MKK_phase3.3_consensus.b36_fwd.txt.gz, 20649 195 
## 	MKK: remove non-standard alleles (14)
## Read: genotypes_chr22_TSI_phase3.3_consensus.b36_fwd.txt.gz, 20649 113 
## 	TSI: remove non-standard alleles (15)
## Read: genotypes_chr22_YRI_phase3.3_consensus.b36_fwd.txt.gz, 20649 214 
## 	YRI: remove non-standard alleles (14)
## 	Position Intersect: 20631 
## 	RS ID Intersect: 20631 
## 	Allele Intersect: 20543 
## 	Genotypes:
##        0        1        2        3 
## 10457440  7777322 10392874    70935 
## 
## Read: genotypes_chrX_ASW_phase3.3_consensus.b36_fwd.txt.gz, 34064 98 
## 	ASW: remove non-standard alleles (3)
## Read: genotypes_chrX_CEU_phase3.3_consensus.b36_fwd.txt.gz, 34064 176 
## 	CEU: remove non-standard alleles (7)
## Read: genotypes_chrX_CHB_phase3.3_consensus.b36_fwd.txt.gz, 34064 148 
## 	CHB: remove non-standard alleles (7)
## Read: genotypes_chrX_CHD_phase3.3_consensus.b36_fwd.txt.gz, 34064 120 
## 	CHD: remove non-standard alleles (7)
## Read: genotypes_chrX_GIH_phase3.3_consensus.b36_fwd.txt.gz, 34064 112 
## 	GIH: remove non-standard alleles (7)
## Read: genotypes_chrX_JPT_phase3.3_consensus.b36_fwd.txt.gz, 34064 124 
## 	JPT: remove non-standard alleles (7)
## Read: genotypes_chrX_LWK_phase3.3_consensus.b36_fwd.txt.gz, 34064 121 
## 	LWK: remove non-standard alleles (3)
## Read: genotypes_chrX_MEX_phase3.3_consensus.b36_fwd.txt.gz, 34064 97 
## 	MEX: remove non-standard alleles (5)
## Read: genotypes_chrX_MKK_phase3.3_consensus.b36_fwd.txt.gz, 34064 195 
## 	MKK: remove non-standard alleles (4)
## Read: genotypes_chrX_TSI_phase3.3_consensus.b36_fwd.txt.gz, 34064 113 
## 	TSI: remove non-standard alleles (7)
## Read: genotypes_chrX_YRI_phase3.3_consensus.b36_fwd.txt.gz, 34064 214 
## 	YRI: remove non-standard alleles (4)
## 	Position Intersect: 34056 
## 	RS ID Intersect: 34056 
## 	Allele Intersect: 33783 
## 	Genotypes:
##        0        1        2        3 
## 20427718  5898143 20532646   336344 
## 
## Read: genotypes_chrY_ASW_phase3.3_consensus.b36_fwd.txt.gz, 490 98 
## Read: genotypes_chrY_CEU_phase3.3_consensus.b36_fwd.txt.gz, 490 176 
## Read: genotypes_chrY_CHB_phase3.3_consensus.b36_fwd.txt.gz, 490 148 
## Read: genotypes_chrY_CHD_phase3.3_consensus.b36_fwd.txt.gz, 490 120 
## Read: genotypes_chrY_GIH_phase3.3_consensus.b36_fwd.txt.gz, 490 112 
## Read: genotypes_chrY_JPT_phase3.3_consensus.b36_fwd.txt.gz, 490 124 
## Read: genotypes_chrY_LWK_phase3.3_consensus.b36_fwd.txt.gz, 490 121 
## Read: genotypes_chrY_MEX_phase3.3_consensus.b36_fwd.txt.gz, 490 97 
## Read: genotypes_chrY_MKK_phase3.3_consensus.b36_fwd.txt.gz, 490 195 
## Read: genotypes_chrY_TSI_phase3.3_consensus.b36_fwd.txt.gz, 490 113 
## Read: genotypes_chrY_YRI_phase3.3_consensus.b36_fwd.txt.gz, 490 214 
## 	Position Intersect: 490 
## 	RS ID Intersect: 490 
## 	Allele Intersect: 432 
## 	Genotypes:
##      0      1      2      3 
## 194659 208667 194797   5381
```

```r
# Add sample.annotation

dat <- read.table(sample.fn, header=TRUE, stringsAsFactors=FALSE)
names(dat) <- c("family", "IID", "father", "mother", "sex", "phenotype", "population")
dat$sex[dat$sex==1] <- "M"
dat$sex[dat$sex==2] <- "F"
dat <- dat[match(sample.id, dat$IID),
	c("family", "father", "mother", "sex", "population")]
table(dat$population, exclude=NULL)
```

```
## 
##  ASW  CEU  CHB  CHD  GIH  JPT  LWK  MEX  MKK  TSI  YRI <NA> 
##   87  165  137  109  101  113  110   86  184  102  203    0
```

```r
add.gdsn(newfile, "sample.annot", dat, compress="ZIP_RA.max", closezip=TRUE)

# show
newfile
```

```
## File: /home/postdoc/zhengx/my/CreateGDS/HapMap3/HapMap3_r3_b36_fwd_consensus_qc_poly.gds
## +    [  ] *
## |--+ sample.id   { VStr8 1397 ZIP_RA(27.43%) }
## |--+ snp.id   { Int32 1452477 ZIP_RA(34.38%) }
## |--+ snp.rs.id   { VStr8 1452477 ZIP_RA(34.11%) }
## |--+ snp.chromosome   { VStr8 1452477 ZIP_RA(0.00%) }
## |--+ snp.position   { Int32 1452477 ZIP_RA(80.41%) }
## |--+ snp.allele   { VStr8 1452477 ZIP_RA(9.07%) }
## |--+ genotype   { Bit2 1397x1452477 ZIP_RA(54.28%) } *
## |--+ sample.annot   [ data.frame ] *
## |  |--+ family   { VStr8 1397 ZIP_RA(26.06%) }
## |  |--+ father   { VStr8 1397 ZIP_RA(15.64%) }
## |  |--+ mother   { VStr8 1397 ZIP_RA(15.97%) }
## |  |--+ sex   { VStr8 1397 ZIP_RA(11.56%) }
## |  |--+ population   { VStr8 1397 ZIP_RA(1.65%) }
```

```r
# Close the GDS file
closefn.gds(newfile)

cleanup.gds(gds.fn)
```

```
## Clean up the fragments of GDS file:
## 	open the file "HapMap3_r3_b36_fwd_consensus_qc_poly.gds" (size: 287748146).
## 	# of fragments in total: 165.
## 	save it to "HapMap3_r3_b36_fwd_consensus_qc_poly.gds.tmp".
## 	rename "HapMap3_r3_b36_fwd_consensus_qc_poly.gds.tmp" (size: 287746478).
## 	# of fragments in total: 26.
```

```r
#########################################################################
##  Summarize the dataset

snpgdsSummary(gds.fn)
```

```
## The file name: /home/postdoc/zhengx/my/CreateGDS/HapMap3/HapMap3_r3_b36_fwd_consensus_qc_poly.gds 
## The total number of samples: 1397 
## The total number of SNPs: 1452477 
## SNP genotypes are stored in SNP-major mode (Sample X SNP).
```

```r
#########################################################################
##  Session Info

sessionInfo()
```

```
## R version 3.1.2 (2014-10-31)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
## [1] SNPRelate_1.0.1 gdsfmt_1.2.2    knitr_1.9      
## 
## loaded via a namespace (and not attached):
## [1] evaluate_0.5.5 formatR_1.0    methods_3.1.2  stringr_0.6.2 
## [5] tools_3.1.2
```
