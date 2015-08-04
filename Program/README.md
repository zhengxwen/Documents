# Program

## viewgds.R

`viewgds` is a shell script written in R, to view the contents of a GDS file. The R packages `gdsfmt`, `getopt` and `optparse` should be installed before running `viewgds`, and the package `crayon` is optional.

```R
install.packages("getopt", repos="http://cran.r-project.org")
install.packages("optparse", repos="http://cran.r-project.org")
install.packages("crayon", repos="http://cran.r-project.org")

source("http://bioconductor.org/biocLite.R")
biocLite("gdsfmt")
```

```sh
## download
wget --no-check-certificate https://raw.githubusercontent.com/zhengxwen/Documents/master/Program/viewgds.R -O viewgds
## or (on mac)
curl -L https://raw.githubusercontent.com/zhengxwen/Documents/master/Program/viewgds.R -o viewgds

## make it executable
chmod +x viewgds
```
