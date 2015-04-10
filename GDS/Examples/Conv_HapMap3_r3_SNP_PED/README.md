
```r
#########################################################################
##
##  Convert HapMap3_r3 Genotypes to SNP GDS Format from PED format
##
##  File: Conv_HapMap3_r3_SNP.R
##  Output: HapMap3_r3_b36_fwd_consensus_qc_poly.gds
##

library(gdsfmt)
library(SNPRelate)
```

```
## SNPRelate -- supported by Streaming SIMD Extensions 2 (SSE2)
```

```r
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
#########################################################################
##  Create GDS file

gds.fn <- "HapMap3_r3_b36_fwd_consensus_qc_poly.gds"

snpgdsPED2GDS(FTP_PED_FILE, FTP_MAP_FILE, gds.fn, family=FALSE,
	snpfirstdim=FALSE, compress.annotation="ZIP_RA.max",
	compress.geno="ZIP_RA.max", verbose=TRUE)
```

```
## PLINK PED/MAP to GDS Format:
## Import 1457897 variants from 'hapmap3_r3_b36_fwd.consensus.qc.poly.map.gz'
## Chromosome:
##      1      2      3      4      5      6      7      8      9     10 
## 119487 119502  98971  88135  90368  93671  77377  77111  65251  75616 
##     11     12     13     14     15     16     17     18     19     20 
##  72993  70482  53293  46655  43309  45778  39329  41942  26953  37159 
##     21     22     23     25 
##  19802  20649  33574    490 
## Reading 'hapmap3_r3_b36_fwd.consensus.qc.poly.ped.gz'
## Output: 'HapMap3_r3_b36_fwd_consensus_qc_poly.gds'
## Import 1397 samples
## Transpose the genotypic matrix ...
## Done.
## Optimize the access efficiency ...
## Clean up the fragments of GDS file:
## 	open the file "HapMap3_r3_b36_fwd_consensus_qc_poly.gds" (size: 776935439).
## 	# of fragments in total: 34.
## 	save it to "HapMap3_r3_b36_fwd_consensus_qc_poly.gds.tmp".
## 	rename "HapMap3_r3_b36_fwd_consensus_qc_poly.gds.tmp" (size: 267767746).
## 	# of fragments in total: 15.
```

```r
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
```

```
## 
##  ASW  CEU  CHB  CHD  GIH  JPT  LWK  MEX  MKK  TSI  YRI <NA> 
##   87  165  137  109  101  113  110   86  184  102  203    0
```

```r
add.gdsn(gfile, "sample.annot", dat, compress="ZIP_RA.max", closezip=TRUE)

# show
gfile
```

```
## File: /home/postdoc/zhengx/tmp/HapMap3/HapMap3_r3_b36_fwd_consensus_qc_poly.gds
## +    [  ] *
## |--+ sample.id   { VStr8 1397 ZIP_RA(28.15%) }
## |--+ snp.id   { Int32 1457897 ZIP_RA(34.58%) }
## |--+ snp.rs.id   { VStr8 1457897 ZIP_RA(34.28%) }
## |--+ snp.position   { Int32 1457897 ZIP_RA(80.98%) }
## |--+ snp.chromosome   { Int32 1457897 ZIP_RA(0.10%) } *
## |--+ snp.allele   { VStr8 1457897 ZIP_RA(14.07%) }
## |--+ genotype   { Bit2 1397x1457897 ZIP_RA(50.10%) } *
## |--+ sample.annot   [ data.frame ] *
## |  |--+ family   { VStr8 1397 ZIP_RA(29.66%) }
## |  |--+ father   { VStr8 1397 ZIP_RA(15.87%) }
## |  |--+ mother   { VStr8 1397 ZIP_RA(16.27%) }
## |  |--+ sex   { VStr8 1397 ZIP_RA(12.13%) }
## |  |--+ population   { VStr8 1397 ZIP_RA(1.65%) }
```

```r
# Close the GDS file
closefn.gds(gfile)

cleanup.gds(gds.fn)
```

```
## Clean up the fragments of GDS file:
## 	open the file "HapMap3_r3_b36_fwd_consensus_qc_poly.gds" (size: 267773373).
## 	# of fragments in total: 38.
## 	save it to "HapMap3_r3_b36_fwd_consensus_qc_poly.gds.tmp".
## 	rename "HapMap3_r3_b36_fwd_consensus_qc_poly.gds.tmp" (size: 267773229).
## 	# of fragments in total: 26.
```

```r
#########################################################################
##  Summarize the dataset

snpgdsSummary(gds.fn)
```

```
## Some of 'snp.allele' are not standard (e.g., I/D).
```

```
## The file name: /home/postdoc/zhengx/tmp/HapMap3/HapMap3_r3_b36_fwd_consensus_qc_poly.gds 
## The total number of samples: 1397 
## The total number of SNPs: 1457897 
## SNP genotypes are stored in SNP-major mode (Sample X SNP).
## The number of valid samples: 1397 
## The number of biallelic unique SNPs: 1457554
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
## [1] SNPRelate_1.1.11 gdsfmt_1.3.10    knitr_1.9       
## 
## loaded via a namespace (and not attached):
## [1] evaluate_0.5.5 formatR_1.0    methods_3.1.2  stringr_0.6.2 
## [5] tools_3.1.2
```

