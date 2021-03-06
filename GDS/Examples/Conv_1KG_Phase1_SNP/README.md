
```r
#########################################################################
##
##  Convert 1KG Phase 1 Genotypes to SNP GDS Format
##
##  File: Conv_1KG_Phase1_SNP.R
##  Output: ALL_chr_phase1_release_v3_20101123_snps_indels_svs_genotypes.snp.gds
##

library(gdsfmt)
library(SNPRelate)
```

```
## SNPRelate -- supported by Streaming SIMD Extensions 2 (SSE2)
```

```r
# file name
FTP_BASE <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp"
FTP_GENO_PATH <- "release/20110521"
FTP_FILE_TEMPLATE <- "ALL.chr%s.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"


#########################################################################
##  Create GDS file

vcf.fn <- sprintf(FTP_FILE_TEMPLATE, c(1:22, "X"))
gds.fn <- "ALL_chr_phase1_release_v3_20101123_snps_indels_svs_genotypes.snp.gds"

snpgdsVCF2GDS(vcf.fn, gds.fn, method="copy.num.of.ref", snpfirstdim=FALSE,
	compress.annotation="ZIP_RA.max", compress.geno="ZIP_RA.max")
```

```
## VCF Format --> SNP GDS Format
## Method: dosage (0,1,2) of reference allele for all variant sites
## Number of samples: 1092
## Parsing "ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 3007196 variants.
## + genotype   { Bit2 1092x3007196 ZIP_RA(7.78%) } *
## Parsing "ALL.chr2.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 3307592 variants.
## + genotype   { Bit2 1092x6314788 ZIP_RA(7.75%) } *
## Parsing "ALL.chr3.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 2763454 variants.
## + genotype   { Bit2 1092x9078242 ZIP_RA(7.74%) } *
## Parsing "ALL.chr4.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 2736765 variants.
## + genotype   { Bit2 1092x11815007 ZIP_RA(7.74%) } *
## Parsing "ALL.chr5.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 2530217 variants.
## + genotype   { Bit2 1092x14345224 ZIP_RA(7.72%) } *
## Parsing "ALL.chr6.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 2424425 variants.
## + genotype   { Bit2 1092x16769649 ZIP_RA(7.76%) } *
## Parsing "ALL.chr7.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 2215231 variants.
## + genotype   { Bit2 1092x18984880 ZIP_RA(7.81%) } *
## Parsing "ALL.chr8.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 2183839 variants.
## + genotype   { Bit2 1092x21168719 ZIP_RA(7.80%) } *
## Parsing "ALL.chr9.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 1652388 variants.
## + genotype   { Bit2 1092x22821107 ZIP_RA(7.84%) } *
## Parsing "ALL.chr10.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 1882663 variants.
## + genotype   { Bit2 1092x24703770 ZIP_RA(7.85%) } *
## Parsing "ALL.chr11.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 1894908 variants.
## + genotype   { Bit2 1092x26598678 ZIP_RA(7.86%) } *
## Parsing "ALL.chr12.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 1828006 variants.
## + genotype   { Bit2 1092x28426684 ZIP_RA(7.87%) } *
## Parsing "ALL.chr13.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 1373000 variants.
## + genotype   { Bit2 1092x29799684 ZIP_RA(7.87%) } *
## Parsing "ALL.chr14.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 1258254 variants.
## + genotype   { Bit2 1092x31057938 ZIP_RA(7.88%) } *
## Parsing "ALL.chr15.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 1130554 variants.
## + genotype   { Bit2 1092x32188492 ZIP_RA(7.90%) } *
## Parsing "ALL.chr16.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 1210619 variants.
## + genotype   { Bit2 1092x33399111 ZIP_RA(7.93%) } *
## Parsing "ALL.chr17.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 1046733 variants.
## + genotype   { Bit2 1092x34445844 ZIP_RA(7.94%) } *
## Parsing "ALL.chr18.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 1088820 variants.
## + genotype   { Bit2 1092x35534664 ZIP_RA(7.95%) } *
## Parsing "ALL.chr19.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 816115 variants.
## + genotype   { Bit2 1092x36350779 ZIP_RA(7.98%) } *
## Parsing "ALL.chr20.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 855166 variants.
## + genotype   { Bit2 1092x37205945 ZIP_RA(7.99%) } *
## Parsing "ALL.chr21.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 518965 variants.
## + genotype   { Bit2 1092x37724910 ZIP_RA(8.00%) } *
## Parsing "ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 494328 variants.
## + genotype   { Bit2 1092x38219238 ZIP_RA(8.01%) } *
## Parsing "ALL.chrX.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" ...
## 	import 1487477 variants.
## + genotype   { Bit2 1092x39706715 ZIP_RA(7.98%) } *
## Optimize the access efficiency ...
## Clean up the fragments of GDS file:
## 	open the file "ALL_chr_phase1_release_v3_20101123_snps_indels_svs_genotypes.snp.gds" (size: 1216036222).
## 	# of fragments in total: 13706.
## 	save it to "ALL_chr_phase1_release_v3_20101123_snps_indels_svs_genotypes.snp.gds.tmp".
## 	rename "ALL_chr_phase1_release_v3_20101123_snps_indels_svs_genotypes.snp.gds.tmp" (size: 1215871990).
## 	# of fragments in total: 20.
```

```r
#########################################################################
##  Summarize the dataset

snpgdsSummary(gds.fn)
```

```
## Some of 'snp.allele' are not standard (e.g., TC/T).
```

```
## The file name: ZhengX/1KG/Phase1/ALL_chr_phase1_release_v3_20101123_snps_indels_svs_genotypes.snp.gds 
## The total number of samples: 1092 
## The total number of SNPs: 39706715 
## SNP genotypes are stored in SNP-major mode (Sample X SNP).
## The number of valid samples: 1092 
## The number of biallelic unique SNPs: 38248779
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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] SNPRelate_1.1.11 gdsfmt_1.3.10    knitr_1.9       
## 
## loaded via a namespace (and not attached):
## [1] evaluate_0.5.5 formatR_1.0    stringr_0.6.2  tools_3.1.2
```

