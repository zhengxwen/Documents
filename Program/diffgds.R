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
		help="Clean up the fragments of GDS file"),
	make_option(c("-v", "--version", action="store_true", default=FALSE,
		help="Show version")
)
parser <- OptionParser(usage="%prog [options] file1 file2",
	option_list=option_list)

arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
files <- arguments$args



####  define the internal functions  ####

hash <- "sha512"

ERROR <- if (crayon.flag) crayon::red else `c`
INVERSE <- if (crayon.flag) crayon::inverse else `c`
OK <- if (crayon.flag) crayon::blue("OK") else "OK"

common_nodes <- character()
file1_nodes <- character()
file2_nodes <- character()


scan_gdsn <- function(node1, node2, path)
{
	n1 <- ls.gdsn(node1, include.hidden=TRUE)
	n2 <- ls.gdsn(node2, include.hidden=TRUE)
	nn <- intersect(n1, n2)

	s <- setdiff(n1, nn)
	if (length(s) > 0L)
		file1_nodes <<- c(file1_nodes, paste0(path, s))
	s <- setdiff(n2, nn)
	if (length(s) > 0L)
		file2_nodes <<- c(file2_nodes, paste0(path, s))

	for (n in nn)
	{
		i1 <- index.gdsn(node1, n)
		i2 <- index.gdsn(node2, n)
		d1 <- objdesp.gdsn(i1)$dim
		d2 <- objdesp.gdsn(i2)$dim
		fullname <- name.gdsn(i1, TRUE)

		if (!is.null(d1) & !is.null(d2))
		{
			cat("checking ", hash, ", ", fullname, " ...", sep="")

			flag <- identical(d1, d2)
			if (flag)
			{
				if (.Platform$OS.type == "windows")
				{
					m1 <- digest.gdsn(i1, algo=hash, action="Robject")
					m2 <- digest.gdsn(i2, algo=hash, action="Robject")
				} else {
					v <- unlist(mclapply(list(i1, i2), FUN=function(x) {
						digest.gdsn(x, algo=hash, action="Robject")
					}, mc.cores=2L))
					m1 <- v[1L]; m2 <- v[2L]
				}
				flag <- identical(m1, m2)
				cat(ifelse(flag,
					paste0("\rdim, ", hash, " [", OK, "]"),
					paste0("\rdim [OK], ", hash, " [", ERROR("error"), "]")))
			} else {
				cat("\rdim [", ERROR("error"), "]", sep="")
			}

			cat("  ")
			if (flag)
				print(i1, all=TRUE)
			else
				cat(fullname, "\n", sep="")
		}

		scan_gdsn(i1, i2, paste0(path, n, "/"))
	}
}



####  run the main part  ####

main <- function()
{
	if (opt$version)
	{
		cat("diffgds 0.9.0\n")
		return(invisible())
	}

	if (length(files) != 2L)
		stop("No input file, see \"diffgds -h\".")	

	# open the files
	f1 <- openfn.gds(files[1L], allow.fork=TRUE)
	f2 <- openfn.gds(files[2L], allow.fork=TRUE)

	# scan files
	cat(INVERSE(">>>> common"), "\n", sep="")
	scan_gdsn(f1$root, f2$root, "")

	cat(INVERSE(paste0("<<<< ", basename(files[1L]))), "\n", sep="")
	for (nm in file1_nodes)
    	print(index.gdsn(f1, nm), expand=FALSE)

	cat(INVERSE(paste0("<<<< ", basename(files[2L]))), "\n", sep="")
	for (nm in file2_nodes)
    	print(index.gdsn(f2, nm), expand=FALSE)

	# close the files
	closefn.gds(f1)
	closefn.gds(f2)

	invisible()
}


res <- try(main())

# quit
q(status = ifelse(inherits(res, "try-error"), 1L, 0L))
