# Program

## viewgds.R

`viewgds` is a shell script written in R, to view the contents of a GDS file. The R packages `gdsfmt`, `getopt` and `optparse` should be installed before running `viewgds`, and the package `crayon` is optional.

In the R environment,
```R
install.packages("getopt", repos="http://cran.r-project.org")
install.packages("optparse", repos="http://cran.r-project.org")
install.packages("crayon", repos="http://cran.r-project.org")

source("http://bioconductor.org/biocLite.R")
biocLite("gdsfmt")
```

Command line,
```sh
## download
wget --no-check-certificate https://raw.githubusercontent.com/zhengxwen/Documents/master/Program/viewgds.R -O viewgds
## or (on mac)
curl -L https://raw.githubusercontent.com/zhengxwen/Documents/master/Program/viewgds.R -o viewgds

## make it executable
chmod +x viewgds
```

You might have to modify the first line of `viewgds` to locate the program `Rscript`:
```R
#! /usr/bin/Rscript --vanilla
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("gdsfmt"))
```
E.g., change `/usr/bin/Rscript` to `/usr/local/bin/Rscript`.

```
Usage: viewgds [options] file1 [file2] [file3]

Options:
	-a, --all
		Include hidden GDS node(s)

	-n NODE, --node=NODE
		Specify a GDS node

	-e FILENAME, --export=FILENAME
		Export a GDS node to a text file (-e "" for standard output)

	-i FILENAME, --import=FILENAME
		Import a text file and create a GDS node (-i "" for standard input)

	-c FORMAT, --create=FORMAT
		Create a GDS node with the format TYPE:DIM:COMPRESSION (e.g., -n NAME -c "int:4,0:ZIP_RA.max")

	--delete
		Delete the GDS node (e.g., -n NAME --delete)

	--newfile
		Create a new GDS file

	--show-attr
		Show the attribute(s)

	--quiet
		No screen output

	--clean
		Clean up the fragments of GDS file

	-h, --help
		Show this help message and exit
```
