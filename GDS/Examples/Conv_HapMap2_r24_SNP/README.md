```r
#########################################################################
##
##  Convert HapMap2_r24 Genotypes to SNP GDS Format
##
##  File: Conv_HapMap2_r24_SNP.R
##  Output: HapMap2_r24_nr_b36_fwd.gds
##

library(gdsfmt)
library(SNPRelate)
```

```
## SNPRelate -- supported by Streaming SIMD Extensions 2 (SSE2)
```

```r
# file name template
FILE_TEMPLATE <- "genotypes_chr%s_%s_r24_nr.b36_fwd.txt.gz"



#########################################################################
##  Download the gz files from HapMap Website

FTP_BASE <- "ftp://ftp.ncbi.nlm.nih.gov/hapmap/"
FTP_GENO_PATH <- "genotypes/2008-10_phaseII/fwd_strand/non-redundant"
POP_LIST <- c("CEU", "JPT+CHB", "YRI")

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
FTP_SAMP_PATH <- "samples_individuals"
sample.fn <- "relationships_w_pops_121708.txt"
download.file(file.path(FTP_BASE, FTP_SAMP_PATH, sample.fn), sample.fn)

head(read.table(sample.fn, header=TRUE, stringsAsFactors=FALSE))
```

```
##    FID     IID     dad     mom sex pheno population
## 1 2357 NA19625       0       0   2     0        ASW
## 2 2367 NA19702 NA19700 NA19701   1     0        ASW
## 3 2367 NA19700       0       0   1     0        ASW
## 4 2367 NA19701       0       0   2     0        ASW
## 5 2368 NA19705 NA19703 NA19704   1     0        ASW
## 6 2368 NA19703       0       0   1     0        ASW
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
```

```
## Starting ...
```

```r
# Write SNPs into the GDS file
for (chr.id in c(1:22, "X", "Y", "M"))
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
		allele <- dat[[pop]]$SNPalleles
		s <- strsplit(allele, "/")
		flag <- sapply(s, function(x)
			(length(x)==2L) & all(x %in% c("A", "G", "C", "T")))
		if (sum(!flag) > 0)
		{
			dat[[pop]] <- dat[[pop]][flag, ]
			allele <- dat[[pop]]$SNPalleles
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
		flag <- duplicated(dat[[pop]]$rs.)
		if (any(flag))
		{
			dat[[pop]] <- dat[[pop]][!flag, ]
			cat("\tRS IDs in Population", pop, "are not unique.\n")
		}
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
	allele <- dat[[POP_LIST[1]]]$SNPalleles
	flag <- rep(TRUE, length(allele))
	for (pop in POP_LIST[-1])
	{
		flag <- flag & (allele == dat[[pop]]$SNPalleles)
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
	allele <- dat[[1]]$SNPalleles
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
## Read: genotypes_chr1_CEU_r24_nr.b36_fwd.txt.gz, 307691 101 
## Read: genotypes_chr1_JPT+CHB_r24_nr.b36_fwd.txt.gz, 311854 101 
## Read: genotypes_chr1_YRI_r24_nr.b36_fwd.txt.gz, 305929 101 
## 	Position Intersect: 293961 
## 	RS ID Intersect: 293961 
## 	Allele Intersect: 293961 
## 	Genotypes:
##        0        1        2        3 
## 31353638 15375804 31292006  1348022 
## 
## Read: genotypes_chr2_CEU_r24_nr.b36_fwd.txt.gz, 326231 101 
## Read: genotypes_chr2_JPT+CHB_r24_nr.b36_fwd.txt.gz, 327180 101 
## Read: genotypes_chr2_YRI_r24_nr.b36_fwd.txt.gz, 318602 101 
## 	Position Intersect: 311138 
## 	RS ID Intersect: 311138 
## 	Allele Intersect: 311138 
## 	Genotypes:
##        0        1        2        3 
## 32613209 17657338 32574767  1161946 
## 
## Read: genotypes_chr3_CEU_r24_nr.b36_fwd.txt.gz, 255391 101 
## Read: genotypes_chr3_JPT+CHB_r24_nr.b36_fwd.txt.gz, 255618 101 
## Read: genotypes_chr3_YRI_r24_nr.b36_fwd.txt.gz, 250155 101 
## 	Position Intersect: 240048 
## 	RS ID Intersect: 240048 
## 	Allele Intersect: 240048 
## 	Genotypes:
##        0        1        2        3 
## 25043393 13818481 24903744  1047342 
## 
## Read: genotypes_chr4_CEU_r24_nr.b36_fwd.txt.gz, 244849 101 
## Read: genotypes_chr4_JPT+CHB_r24_nr.b36_fwd.txt.gz, 245102 101 
## Read: genotypes_chr4_YRI_r24_nr.b36_fwd.txt.gz, 238922 101 
## 	Position Intersect: 229888 
## 	RS ID Intersect: 229888 
## 	Allele Intersect: 229888 
## 	Genotypes:
##        0        1        2        3 
## 24073953 12824498 24003056  1168253 
## 
## Read: genotypes_chr5_CEU_r24_nr.b36_fwd.txt.gz, 247632 101 
## Read: genotypes_chr5_JPT+CHB_r24_nr.b36_fwd.txt.gz, 248154 101 
## Read: genotypes_chr5_YRI_r24_nr.b36_fwd.txt.gz, 242186 101 
## 	Position Intersect: 234786 
## 	RS ID Intersect: 234786 
## 	Allele Intersect: 234786 
## 	Genotypes:
##        0        1        2        3 
## 24562687 13303184 24506196  1020153 
## 
## Read: genotypes_chr6_CEU_r24_nr.b36_fwd.txt.gz, 268348 101 
## Read: genotypes_chr6_JPT+CHB_r24_nr.b36_fwd.txt.gz, 272814 101 
## Read: genotypes_chr6_YRI_r24_nr.b36_fwd.txt.gz, 265955 101 
## 	Position Intersect: 256301 
## 	RS ID Intersect: 256301 
## 	Allele Intersect: 256301 
## 	Genotypes:
##        0        1        2        3 
## 27007124 14169112 26961502  1063532 
## 
## Read: genotypes_chr7_CEU_r24_nr.b36_fwd.txt.gz, 213023 101 
## Read: genotypes_chr7_JPT+CHB_r24_nr.b36_fwd.txt.gz, 213891 101 
## Read: genotypes_chr7_YRI_r24_nr.b36_fwd.txt.gz, 208708 101 
## 	Position Intersect: 200508 
## 	RS ID Intersect: 200508 
## 	Allele Intersect: 200508 
## 	Genotypes:
##        0        1        2        3 
## 21029169 11198998 20962212   946781 
## 
## Read: genotypes_chr8_CEU_r24_nr.b36_fwd.txt.gz, 213095 101 
## Read: genotypes_chr8_JPT+CHB_r24_nr.b36_fwd.txt.gz, 216811 101 
## Read: genotypes_chr8_YRI_r24_nr.b36_fwd.txt.gz, 212014 101 
## 	Position Intersect: 202138 
## 	RS ID Intersect: 202138 
## 	Allele Intersect: 202138 
## 	Genotypes:
##        0        1        2        3 
## 21114323 11756971 20954509   751457 
## 
## Read: genotypes_chr9_CEU_r24_nr.b36_fwd.txt.gz, 181445 101 
## Read: genotypes_chr9_JPT+CHB_r24_nr.b36_fwd.txt.gz, 183433 101 
## Read: genotypes_chr9_YRI_r24_nr.b36_fwd.txt.gz, 180147 101 
## 	Position Intersect: 171793 
## 	RS ID Intersect: 171793 
## 	Allele Intersect: 171793 
## 	Genotypes:
##        0        1        2        3 
## 18076595  9589741 18028566   689208 
## 
## Read: genotypes_chr10_CEU_r24_nr.b36_fwd.txt.gz, 209342 101 
## Read: genotypes_chr10_JPT+CHB_r24_nr.b36_fwd.txt.gz, 211862 101 
## Read: genotypes_chr10_YRI_r24_nr.b36_fwd.txt.gz, 204146 101 
## 	Position Intersect: 196093 
## 	RS ID Intersect: 196090 
## 	Allele Intersect: 196090 
## 	Genotypes:
##        0        1        2        3 
## 20848586 10424637 20737507   933570 
## 
## Read: genotypes_chr11_CEU_r24_nr.b36_fwd.txt.gz, 204228 101 
## Read: genotypes_chr11_JPT+CHB_r24_nr.b36_fwd.txt.gz, 205538 101 
## Read: genotypes_chr11_YRI_r24_nr.b36_fwd.txt.gz, 195110 101 
## 	Position Intersect: 189512 
## 	RS ID Intersect: 189511 
## 	Allele Intersect: 189511 
## 	Genotypes:
##        0        1        2        3 
## 20224947  9898324 20139379   905320 
## 
## Read: genotypes_chr12_CEU_r24_nr.b36_fwd.txt.gz, 191979 101 
## Read: genotypes_chr12_JPT+CHB_r24_nr.b36_fwd.txt.gz, 193071 101 
## Read: genotypes_chr12_YRI_r24_nr.b36_fwd.txt.gz, 187294 101 
## 	Position Intersect: 175487 
## 	RS ID Intersect: 175487 
## 	Allele Intersect: 175487 
## 	Genotypes:
##        0        1        2        3 
## 18711609  9171176 18627721   870984 
## 
## Read: genotypes_chr13_CEU_r24_nr.b36_fwd.txt.gz, 155905 101 
## Read: genotypes_chr13_JPT+CHB_r24_nr.b36_fwd.txt.gz, 158406 101 
## Read: genotypes_chr13_YRI_r24_nr.b36_fwd.txt.gz, 152674 101 
## 	Position Intersect: 146362 
## 	RS ID Intersect: 146359 
## 	Allele Intersect: 146359 
## 	Genotypes:
##        0        1        2        3 
## 15406596  7939618 15508206   662510 
## 
## Read: genotypes_chr14_CEU_r24_nr.b36_fwd.txt.gz, 123071 101 
## Read: genotypes_chr14_JPT+CHB_r24_nr.b36_fwd.txt.gz, 123764 101 
## Read: genotypes_chr14_YRI_r24_nr.b36_fwd.txt.gz, 118518 101 
## 	Position Intersect: 114921 
## 	RS ID Intersect: 114921 
## 	Allele Intersect: 114921 
## 	Genotypes:
##        0        1        2        3 
## 12002642  6454026 12117533   454469 
## 
## Read: genotypes_chr15_CEU_r24_nr.b36_fwd.txt.gz, 106814 101 
## Read: genotypes_chr15_JPT+CHB_r24_nr.b36_fwd.txt.gz, 107363 101 
## Read: genotypes_chr15_YRI_r24_nr.b36_fwd.txt.gz, 102431 101 
## 	Position Intersect: 99440 
## 	RS ID Intersect: 99440 
## 	Allele Intersect: 99440 
## 	Genotypes:
##        0        1        2        3 
## 10396930  5650682 10373528   427660 
## 
## Read: genotypes_chr16_CEU_r24_nr.b36_fwd.txt.gz, 109692 101 
## Read: genotypes_chr16_JPT+CHB_r24_nr.b36_fwd.txt.gz, 109734 101 
## Read: genotypes_chr16_YRI_r24_nr.b36_fwd.txt.gz, 104530 101 
## 	Position Intersect: 101193 
## 	RS ID Intersect: 101193 
## 	Allele Intersect: 101193 
## 	Genotypes:
##        0        1        2        3 
## 10574295  5484275 10756766   506774 
## 
## Read: genotypes_chr17_CEU_r24_nr.b36_fwd.txt.gz, 89701 101 
## Read: genotypes_chr17_JPT+CHB_r24_nr.b36_fwd.txt.gz, 89576 101 
## Read: genotypes_chr17_YRI_r24_nr.b36_fwd.txt.gz, 85541 101 
## 	Position Intersect: 82941 
## 	RS ID Intersect: 82941 
## 	Allele Intersect: 82941 
## 	Genotypes:
##       0       1       2       3 
## 8792301 4504001 8703914  393854 
## 
## Read: genotypes_chr18_CEU_r24_nr.b36_fwd.txt.gz, 119118 101 
## Read: genotypes_chr18_JPT+CHB_r24_nr.b36_fwd.txt.gz, 120025 101 
## Read: genotypes_chr18_YRI_r24_nr.b36_fwd.txt.gz, 115768 101 
## 	Position Intersect: 110859 
## 	RS ID Intersect: 110859 
## 	Allele Intersect: 110859 
## 	Genotypes:
##        0        1        2        3 
## 11728352  5962780 11762914   477884 
## 
## Read: genotypes_chr19_CEU_r24_nr.b36_fwd.txt.gz, 56607 101 
## Read: genotypes_chr19_JPT+CHB_r24_nr.b36_fwd.txt.gz, 56687 101 
## Read: genotypes_chr19_YRI_r24_nr.b36_fwd.txt.gz, 53766 101 
## 	Position Intersect: 51739 
## 	RS ID Intersect: 51739 
## 	Allele Intersect: 51739 
## 	Genotypes:
##       0       1       2       3 
## 5416329 2877886 5435145  240170 
## 
## Read: genotypes_chr20_CEU_r24_nr.b36_fwd.txt.gz, 119921 101 
## Read: genotypes_chr20_JPT+CHB_r24_nr.b36_fwd.txt.gz, 119989 101 
## Read: genotypes_chr20_YRI_r24_nr.b36_fwd.txt.gz, 115921 101 
## 	Position Intersect: 112857 
## 	RS ID Intersect: 112857 
## 	Allele Intersect: 112857 
## 	Genotypes:
##        0        1        2        3 
## 12511508  4877837 12678425   403620 
## 
## Read: genotypes_chr21_CEU_r24_nr.b36_fwd.txt.gz, 50165 101 
## Read: genotypes_chr21_JPT+CHB_r24_nr.b36_fwd.txt.gz, 51900 101 
## Read: genotypes_chr21_YRI_r24_nr.b36_fwd.txt.gz, 49154 101 
## 	Position Intersect: 45882 
## 	RS ID Intersect: 45882 
## 	Allele Intersect: 45882 
## 	Genotypes:
##       0       1       2       3 
## 4789336 2653822 4751773  193209 
## 
## Read: genotypes_chr22_CEU_r24_nr.b36_fwd.txt.gz, 54786 101 
## Read: genotypes_chr22_JPT+CHB_r24_nr.b36_fwd.txt.gz, 56716 101 
## Read: genotypes_chr22_YRI_r24_nr.b36_fwd.txt.gz, 54840 101 
## 	Position Intersect: 51386 
## 	RS ID Intersect: 51386 
## 	Allele Intersect: 51386 
## 	Genotypes:
##       0       1       2       3 
## 5606161 2516052 5549585  202422 
## 
## Read: genotypes_chrX_CEU_r24_nr.b36_fwd.txt.gz, 119769 101 
## 	CEU: remove non-standard alleles (26)
## Read: genotypes_chrX_JPT+CHB_r24_nr.b36_fwd.txt.gz, 120666 101 
## 	JPT+CHB: remove non-standard alleles (20)
## Read: genotypes_chrX_YRI_r24_nr.b36_fwd.txt.gz, 118973 101 
## 	YRI: remove non-standard alleles (15)
## 	Position Intersect: 115219 
## 	RS IDs in Population CEU are not unique.
## 	RS IDs in Population JPT+CHB are not unique.
## 	RS IDs in Population YRI are not unique.
## 	RS ID Intersect: 115212 
## 	Allele Intersect: 115212 
## 	Genotypes:
##        0        1        2        3 
## 14021749  2748295 13873765   463431 
## 
## Read: genotypes_chrY_CEU_r24_nr.b36_fwd.txt.gz, 672 101 
## Read: genotypes_chrY_JPT+CHB_r24_nr.b36_fwd.txt.gz, 667 101 
## Read: genotypes_chrY_YRI_r24_nr.b36_fwd.txt.gz, 668 101 
## 	Position Intersect: 651 
## 	RS ID Intersect: 651 
## 	Allele Intersect: 651 
## 	Genotypes:
##     0     1     2     3 
## 51694 31808 54835 37433 
## 
## Read: genotypes_chrM_CEU_r24_nr.b36_fwd.txt.gz, 212 101 
## Read: genotypes_chrM_JPT+CHB_r24_nr.b36_fwd.txt.gz, 206 101 
## Read: genotypes_chrM_YRI_r24_nr.b36_fwd.txt.gz, 211 101 
## 	Position Intersect: 203 
## 	RS ID Intersect: 203 
## 	Allele Intersect: 203 
## 	Genotypes:
##     0     1     2     3 
## 30195    23 23990   602
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
##  CEU  CHB  JPT  YRI <NA> 
##   90   45   45   90    0
```

```r
add.gdsn(newfile, "sample.annot", dat, compress="ZIP_RA.max", closezip=TRUE)

# show
newfile
```

```
## File: /home/postdoc/zhengx/my/CreateGDS/HapMap2/HapMap2_r24_nr_b36_fwd.gds
## +    [  ] *
## |--+ sample.id   { VStr8 270 ZIP_RA(26.06%) }
## |--+ snp.id   { Int32 3735292 ZIP_RA(34.48%) }
## |--+ snp.rs.id   { VStr8 3735292 ZIP_RA(32.93%) }
## |--+ snp.chromosome   { VStr8 3735292 ZIP_RA(0.09%) }
## |--+ snp.position   { Int32 3735292 ZIP_RA(73.35%) }
## |--+ snp.allele   { VStr8 3735292 ZIP_RA(10.36%) }
## |--+ genotype   { Bit2 270x3735292 ZIP_RA(39.10%) } *
## |--+ sample.annot   [ data.frame ] *
## |  |--+ family   { VStr8 270 ZIP_RA(23.46%) }
## |  |--+ father   { VStr8 270 ZIP_RA(23.89%) }
## |  |--+ mother   { VStr8 270 ZIP_RA(23.22%) }
## |  |--+ sex   { VStr8 270 ZIP_RA(18.89%) }
## |  |--+ population   { VStr8 270 ZIP_RA(4.63%) }
```

```r
# Close the GDS file
closefn.gds(newfile)

cleanup.gds(gds.fn)
```

```
## Clean up the fragments of GDS file:
## 	open the file "HapMap2_r24_nr_b36_fwd.gds" (size: 128987494).
## 	# of fragments in total: 175.
## 	save it to "HapMap2_r24_nr_b36_fwd.gds.tmp".
## 	rename "HapMap2_r24_nr_b36_fwd.gds.tmp" (size: 128985706).
## 	# of fragments in total: 26.
```

```r
#########################################################################
##  Summarize the dataset

snpgdsSummary(gds.fn)
```

```
## The file name: /home/postdoc/zhengx/my/CreateGDS/HapMap2/HapMap2_r24_nr_b36_fwd.gds 
## The total number of samples: 270 
## The total number of SNPs: 3735292 
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
