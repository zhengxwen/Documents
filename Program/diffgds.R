#! /usr/local/bin/Rscript --vanilla
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("gdsfmt"))
suppressPackageStartupMessages(library("parallel"))
crayon.flag <- requireNamespace("crayon", quietly=TRUE)


####  define the options and parse  ####

option_list <- list(
	make_option(c("-a", "--all"), action="store_true", default=FALSE,
		help="Include hidden GDS node(s)"),
	make_option(c("-n", "--node"), action="store", type="character",
		help="Specify a GDS node (e.g., -n VAR1[,VAR2,VAR3...]"),
	make_option(c("-q", "--quiet"), action="store_true", default=FALSE,
		help="No screen output"),
	make_option("--clean", action="store_true", default=FALSE,
		help="Clean up the fragments of GDS file")
)
parser <- OptionParser(usage="%prog [options] file1 file2",
	option_list=option_list)

arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
files <- arguments$args



####  define the internal functions  ####

ERROR <- if (crayon.flag) crayon::red else `c`

common_nodes <- character()
file1_nodes <- character()
file2_nodes <- character()


scan_gdsn <- function(node1, node2)
{
	n1 <- ls.gdsn(node1, include.hidden=TRUE)
	n2 <- ls.gdsn(node1, include.hidden=TRUE)

	nn <- intersect(n1, n2)
	file1_nodes <<- c(file1_nodes, setdiff(n1, nn))
	file2_nodes <<- c(file2_nodes, setdiff(n2, nn))

	for (n in nn)
	{
		i1 <- index.gdsn(node1, n)
		i2 <- index.gdsn(node2, n)
		d1 <- objdesp.gdsn(i1)$dim
		d2 <- objdesp.gdsn(i2)$dim
		fullname <- name.gdsn(i1, TRUE)

		if (!is.null(d1) & !is.null(d2))
		{
			cat("checking ", fullname, " ...", sep="")

			flag <- identical(d1, d2)
			if (flag)
			{
				if (.Platform$OS.type == "windows")
				{
					m1 <- digest.gdsn(i1, algo="md5", action="Robject")
					m2 <- digest.gdsn(i2, algo="md5", action="Robject")
				} else {
					v <- unlist(mclapply(list(i1, i2), FUN=function(x) {
						digest.gdsn(x, algo="md5", action="Robject")
					}, mc.cores=2L))
					m1 <- v[1L]; m2 <- v[2L]
				}
				flag <- identical(m1, m2)
				cat(ifelse(flag, "\rdim, md5 [OK]",
					paste0("\rdim [OK], md5 [", ERROR("error"), "]")))
			} else {
				cat("\rdim [", ERROR("error"), "]", sep="")
			}

			cat("  ")
			if (flag)
				print(i1, all=TRUE)
			else
				cat(fullname, "\n", sep="")
		}

		scan_gdsn(i1, i2)
	}
}



####  run the main part  ####

main <- function()
{
	if (length(files) != 2L)
		stop("No input file, see \"diffgds -h\".")	

	# open the files
	f1 <- openfn.gds(files[1L], allow.fork=TRUE)
	f2 <- openfn.gds(files[2L], allow.fork=TRUE)

	# scan files
	cat(">>>> common\n")
	scan_gdsn(f1$root, f2$root)

	cat("<<<< ", files[1L], "\n", sep="")
    if (length(file1_nodes))
    	cat(paste0(file1_nodes, "\n"))

	cat("<<<< ", files[2L], "\n", sep="")
	if (length(file2_nodes))
		cat(paste0(file2_nodes, "\n"))

	# close the files
	closefn.gds(f1)
	closefn.gds(f2)

	invisible()
}


res <- try(main())
v <- warnings()
if (!is.null(v)) cat(v, file=stderr(), sep="\n")

# quit
q(status = ifelse(inherits(res, "try-error"), 1L, 0L))
