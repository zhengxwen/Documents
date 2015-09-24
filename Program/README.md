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

Installation with command line,
```sh
echo '#!' `which Rscript` '--vanilla' > viewgds
curl -L https://raw.githubusercontent.com/zhengxwen/Documents/master/Program/viewgds.R >> viewgds
chmod +x viewgds

## Or
echo '#!' `which Rscript` '--vanilla' > viewgds
wget -qO- --no-check-certificate https://raw.githubusercontent.com/zhengxwen/Documents/master/Program/viewgds.R >> viewgds
chmod +x viewgds
```

```
Usage: viewgds [options] file1 [file2] [file3]

Options:
	-a, --all
		Include hidden GDS node(s)
	-n NODE, --node=NODE
		Specify a GDS node (e.g., -n VAR1[,VAR2,VAR3...]
	-f FUNCTION, --fun=FUNCTION
		Specify a function for data processing, e.g., --fun "table(x)" for tabulation
	-e FILENAME, --export=FILENAME
		Export a GDS node to a text file (-e "" for standard output)
	--efun=FUNCTION
		Specify a function for exporting data
		E.g., --efun "function(s) paste('prefix', s)" for data preparation
	-i FILENAME, --import=FILENAME
		Import a text file and create a GDS node "-c" (-i "" for standard input)
	--icol=COLUMN
		Specify the columns for importing data (multiple variables splited by semicolon)
		E.g., --icol "-(1:4)" for excluding the first four columns
	--ifun=FUNCTION
		Specify the functions for importing data (multiple functions splited by pound)
		E.g., --ifun "function(x) x+1 # function(x) 4-x" for data preparation
	--iskip=NUMBER
		Specify the number of line skipped when importing data
	--inmax=NUMBER
		Specify the maximun number of line when importing data
	-s CHARACTER, --separator=CHARACTER
		Set the field separator character
	-c FORMAT, --create=FORMAT
		Create a GDS node with the format TYPE:DIM:COMPRESSION[;TYPE2:DIM2:COMPRESSION2;...]
		E.g., -n NAME -c "int:4,0:ZIP_RA.max"
	--delete
		Delete the GDS node (e.g., -n NAME --delete)
	--newfile
		Create a new GDS file
	--show-attr
		Show the attribute(s)
	--attr-set=VALUE
		Set the attribute(s) (e.g., --attr-set "x=1,y=1:4")
	--attr-del=VALUE
		Delete the attribute(s) (e.g., --attr-del "x,y")
	--system
		Show the system configuration
	-q, --quiet
		No screen output
	--clean
		Clean up the fragments of GDS file
	-h, --help
		Show this help message and exit
```
