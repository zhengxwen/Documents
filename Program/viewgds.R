#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("gdsfmt"))
suppressPackageStartupMessages(library("compiler"))
crayon.flag <- requireNamespace("crayon", quietly=TRUE)


####  define the options and parse  ####

option_list <- list(
	make_option(c("-a", "--all"), action="store_true", default=FALSE,
		help="Include hidden GDS node(s)"),
	make_option(c("-n", "--node"), action="store", type="character",
		help="Specify a GDS node (e.g., -n VAR1[,VAR2,VAR3...]"),
	make_option(c("-f", "--fun"), type="character",
		help="Specify a function for data processing, e.g., --fun \"table(x)\" for tabulation",
		metavar="function"),
	make_option(c("-e", "--export"), type="character",
		help="Export a GDS node to a text file (-e \"\" for standard output)",
		metavar="filename"),
	make_option(c("--efun"), type="character",
		help="Specify a function for exporting data\n\t\tE.g., --efun \"function(s) paste('prefix', s)\" for data preparation",
		metavar="function"),
	make_option(c("-i", "--import"), type="character",
		help="Import a text file and create a GDS node \"-c\" (-i \"\" for standard input)",
		metavar="filename"),
	make_option(c("--icol"), type="character",
		help="Specify the columns for importing data (multiple variables splited by semicolon)\n\t\tE.g., --icol \"-(1:4)\" for excluding the first four columns",
		metavar="column"),
	make_option(c("--ifun"), type="character",
		help="Specify the functions for importing data (multiple functions splited by pound)\n\t\tE.g., --ifun \"function(x) x+1 # function(x) 4-x\" for data preparation",
		metavar="function"),
	make_option(c("--iskip"), type="integer",
		help="Specify the number of line skipped when importing data",
		metavar="number"),
	make_option(c("--inmax"), type="integer",
		help="Specify the maximun number of line when importing data",
		metavar="number"),
	make_option("--separator", type="character",
		help="Set the field separator character", metavar="character"),
	make_option(c("-c", "--create"), action="store", type="character",
		help="Create a GDS node with the format TYPE:DIM:COMPRESSION[;TYPE2:DIM2:COMPRESSION2;...]\n\t\tE.g., -n NAME -c \"int:4,0:ZIP_RA.max\"",
		metavar="format"),
	make_option("--digest", type="character",
		help="Create hash function digest for R object ignoring compression (ALGO can be md5, sha1, sha256, sha384 or sha512)",
		metavar="algo"),
	make_option("--hash", type="character",
		help="Create hash function digest for raw data field (ALGO can be md5, sha1, sha256, sha384 or sha512)",
		metavar="algo"),
	make_option("--summary", action="store_true", default=FALSE,
		help="Summarize a GDS node"),
	make_option("--delete", action="store_true", default=FALSE,
		help="Delete the GDS node (e.g., -n NAME --delete)"),
	make_option("--newfile", action="store_true", default=FALSE,
		help="Create a new GDS file"),
	make_option(c("-s", "--show-attr"), action="store_true", default=FALSE,
		help="Show the attribute(s)"),
	make_option("--attr-set", type="character",
		help="Set the attribute(s) (e.g., --attr-set \"x=1,y=1:4\")",
		metavar="value"),
	make_option("--attr-del", type="character",
		help="Delete the attribute(s) (e.g., --attr-del \"x,y\")",
		metavar="value"),
	make_option("--system", action="store_true", default=FALSE,
		help="Show the system configuration"),
	make_option(c("-q", "--quiet"), action="store_true", default=FALSE,
		help="No screen output"),
	make_option("--clean", action="store_true", default=FALSE,
		help="Clean up the fragments of GDS file"),
	make_option(c("-v", "--version"), action="store_true", default=FALSE,
		help="Show version")
)
parser <- OptionParser(usage="%prog [options] file1 [file2] [file3]",
	option_list=option_list)

arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
files <- arguments$args



####  define the internal functions  ####

# default
nprev <- 6L
BLURRED <- if (crayon.flag) crayon::blurred else `c`
INVERSE <- if (crayon.flag) crayon::inverse else `c`
UNDERLINE<- if (crayon.flag) crayon::underline else `c`



last.text <- function(s, len)
{
	substring(s, nchar(s)-len+1L, nchar(s))
}


## set/delete the attributes
do_attribute <- function(gdsfile, namelist)
{
	if (!is.null(opt$`attr-set`))
	{
		if (length(namelist) != 1L)
			stop("Please specify only one GDS node for the attribute assignment.")
		n <- index.gdsn(gdsfile, namelist[1L])
		args <- eval(parse(text=paste("list(", opt$`attr-set`, ")")))
		nm <- names(args)
		for (i in seq_len(length(nm)))
			put.attr.gdsn(n, nm[i], args[[i]])
	}

	if (!is.null(opt$`attr-del`))
	{
		if (length(namelist) != 1L)
			stop("Please specify only one GDS node for the attribute assignment.")
		n <- index.gdsn(gdsfile, namelist[1L])
		nm <- unlist(strsplit(opt$`attr-del`, ",", fixed=TRUE))
		delete.attr.gdsn(n, nm)
	}
}

## clean up the GDS file(s)
do_clean <- function()
{
	if (opt$clean)
	{
		for (fn in files)
			cleanup.gds(fn, verbose=!opt$quiet)
	}
}



####  run the main part  ####

main <- function()
{
	if (opt$version)
	{
		cat("viewgds_0.9.1\n")
		cat("gdsfmt_", packageVersion("gdsfmt"), "\n", sep="")
		return(invisible())
	}

	if (opt$system)
	{
		v <- c(package.version=as.character(packageVersion("gdsfmt")),
			system.gds())
		print(v)
	}

	if (opt$newfile)
	{
		for (fn in files)
		{
			f <- createfn.gds(fn)
			if (!opt$quiet)
			{
				cat(INVERSE("New file:"), "")
				cat(f$filename)
			}
			closefn.gds(f)
		}
		if (!opt$quiet) cat("\n")
	}


	if (!is.null(opt$node))
	{
		if (length(files) != 1L)
		{
			if (length(files) == 0L)
				stop("No input GDS file.")
			cat(files, sep="\n")
			stop("Please specify only one GDS file.")
		}

		# open the existing GDS file
		readonly <- is.null(opt$create) & !opt$delete &
			is.null(opt$`attr-set`) & is.null(opt$`attr-del`)
		gfile <- openfn.gds(files, readonly=readonly)
		on.exit({ closefn.gds(gfile) })

		# get GDS variable(s)
		nodenames <- unlist(strsplit(opt$node, ",", fixed=TRUE))
		if (length(nodenames) <= 0L) nodenames <- ""

		if (!is.null(opt$create))
		{
			if (opt$delete)
				stop("\"-c FORMAT --delete\" should not be specified simultaneously.")

			createfmt <- trimws(unlist(strsplit(opt$create, ";", fixed=TRUE)))
			if (length(createfmt) != length(nodenames))
				stop("'-c' should have the same number of elements as '-n'.")

			nodelist <- list()
			for (i in seq_len(length(nodenames)))
			{
				nm <- unlist(strsplit(nodenames[i], "/", fixed=TRUE))
				if (length(nm) > 1L)
					folder <- index.gdsn(gfile, index=nm[-length(nm)])
				else
					folder <- gfile$root

				fmt <- trimws(unlist(strsplit(createfmt[i], ":", fixed=TRUE)))
				if (length(fmt) <= 0L)
					stop("Please specify the format (e.g., -n NAME -c \"int:4,0:ZIP_RA.max\").")

				if (length(fmt) < 2L) fmt[2L] <- "0"
				if (trimws(fmt[2L]) == "") fmt[2L] <- "0"
				dm <- scan(text=fmt[2L], what=integer(), quiet=TRUE, sep=",")

				if (length(fmt) < 3L) fmt[3L] <- ""
				if (length(fmt) > 4L)
					stop("Maximum number of cells in -c FORMAT is 4.")
				if (length(fmt) < 4L) fmt[4L] <- ""

				args <- list(node=folder, name=nm[length(nm)], storage=fmt[1L],
					valdim=dm, compress=fmt[3L], replace=TRUE)
				args <- c(args, eval(parse(text=paste("list(", fmt[4L], ")"))))

				n <- do.call(add.gdsn, args)
				if (!opt$quiet & is.null(opt$import))
					print(n, attribute=TRUE, attribute.trim=FALSE)
				nodelist[[i]] <- n
			}

			do_attribute(gfile, nodenames)

			# import
			if (!is.null(opt$import))
			{
				if (!opt$quiet)
				{
					for (i in seq_len(length(nodelist)))
					{
						cat(INVERSE("GDS node:"), "")
						print(nodelist[[i]], attribute=TRUE, attribute.trim=FALSE)
					}
					cat(INVERSE("Import:"), " ", opt$import, "\n", sep="")
				}

				ic <- vector("list", length(nodenames))
				if (!is.null(opt$icol))
				{
					s <- trimws(unlist(strsplit(opt$icol, ";", fixed=TRUE)))
					if (length(s) != length(nodenames))
					{
						stop("'--icol' should have the same number of elements as ",
							"'-n' (", length(nodenames), ").")
					}
					for (i in seq_len((length(s))))
					{
						if (s[i] != "")
							ic[[i]] <- eval(parse(text=paste("c(", s[i], ")")))
					}
				}

				ifun <- vector("list", length(nodenames))
				for (i in seq_len(length(ifun)))
					ifun[[i]] <- `c`
				if (!is.null(opt$ifun))
				{
					s <- trimws(unlist(strsplit(opt$ifun, "#", fixed=TRUE)))
					if (length(s) > length(nodenames))
					{
						stop("'--ifun' should have the same number of elements as ",
							"'-n' (", length(nodenames), ").")
					}
					if (length(s) < length(nodenames))
					{
						s <- c(s, rep("", length(nodenames)-length(s)))
					}
					for (i in seq_len((length(s))))
					{
						if (s[i] != "")
						{
							ifun[[i]] <- eval(parse(text=s[i]))
							if (!is.function(ifun[[i]]))
								stop("--ifun \"", s[i], "\" should be a function.")
						}
					}
				}

				if (opt$import == "") opt$import <- "stdin"
				infile <- file(opt$import)
				open(infile)
				on.exit({ close(infile) }, add=TRUE)
				nc <- NA_integer_
				k <- 1L

				if (!is.null(opt$iskip))
				{
					i <- opt$iskip
					if (i > 0L)
					{
						while(length(readLines(f, n=1L)) > 0L)
						{
							k <- k + 1L
							i <- i - 1L
							if (i <= 0L) break
						}
					}
				}

				if (!is.null(opt$inmax))
				{
					if (opt$inmax <= 0L)
						stop("--inmax should be > 0.")
					kmax <- k + opt$inmax
				} else
					kmax <- 2147483647L  # 2^31 - 1

				sep <- ""
				if (!is.null(opt$separator))
				{
					sep <- eval(parse(text=paste('"', opt$separator, '"',
						sep="")))
				}

				while(length(s <- readLines(infile, n=1L)) > 0L)
				{
					ss <- scan(text=s, what="", quiet=TRUE, sep=sep,
						strip.white=TRUE)
					if (is.na(nc))
					{
						nc <- length(ss)
					} else {
						if (nc != length(ss))
						{
							stop(sprintf("Line %d should have %d columns.",
								k, length(ss)))
						}
					}

					for (i in seq_len(length(nodelist)))
					{
						f <- ifun[[i]]
						v <- ic[[i]]
						if (is.null(v))
							append.gdsn(nodelist[[i]], f(ss))
						else
							append.gdsn(nodelist[[i]], f(ss[v]))
					}

					k <- k + 1L
					if (k >= kmax) break
				}

				for (i in seq_len(length(nodelist)))
					readmode.gdsn(nodelist[[i]])

				on.exit({ closefn.gds(gfile) })
				close(infile)
			}

			if (!opt$quiet)
			{
				for (i in seq_len(length(nodelist)))
				{
					cat(INVERSE("New GDS node:"), "")
					print(nodelist[[i]], attribute=TRUE, attribute.trim=FALSE)
				}
			}

			on.exit()
			closefn.gds(gfile)
			do_clean()

			q(status=0L)

		} else if (opt$delete)
		{
			for (nm in nodenames)
			{
				n <- index.gdsn(gfile, nm)
				if (!opt$quiet)
				{
					cat(INVERSE("Delete GDS node:"), "")
					print(n, all=opt$all, attribute=TRUE, attribute.trim=FALSE)
				}
				delete.gdsn(n)
			}

			on.exit()
			closefn.gds(gfile)
			do_clean()

			q(status=0L)
		}

		do_attribute(gfile, nodenames)

		for (nm in nodenames)
		{
			node <- index.gdsn(gfile, nm)
			if (!opt$quiet)
			{
				cat(INVERSE("GDS node:"), "")
				print(node, all=opt$all, attribute=TRUE, attribute.trim=FALSE)

				dp <- objdesp.gdsn(node)
				if (dp$is.array & !is.null(dp$dim))
				{
					cat(INVERSE("Preview:\n"))
					show(node)
				}

				if (opt$summary)
				{
					cat(INVERSE("Summary:\n"))
					v <- summarize.gdsn(node)
					if (!is.null(v$decimal))
					{
						d <- v$decimal
						v$decimal <- NULL
						cat(paste(names(v), format(v, justify="none"), sep=": "), sep="\n")
						cat("decimal:", paste(paste(
							sQuote(names(d)), d, sprintf("(%.1f%%)", d/sum(d)*100)),
							collapse="; "))
						cat("\n")
					} else
						cat(paste(names(v), format(v, justify="none"), sep=": "), sep="\n")
				}
			}
			if (!is.null(opt$fun))
			{
				x <- read.gdsn(node)
				fun <- eval(parse(text=opt$fun))
				if (is.function(fun))
					fun(x)
				else
					print(fun)
			}
			if (!is.null(opt$digest))
			{
				if (!opt$quiet)
					cat(INVERSE(paste(opt$digest, "digest of R object:")), "\n", sep="")
				cat(digest.gdsn(node, opt$hash, action="Robject"), "\n", sep="")
			}
			if (!is.null(opt$hash))
			{
				if (!opt$quiet)
					cat(INVERSE(paste(opt$hash, "digest:")), "\n", sep="")
				cat(digest.gdsn(node, opt$digest), "\n", sep="")
			}
		}

		if (!is.null(opt$export))
		{
			if (length(nodenames) > 1L)
				stop("Only a GDS node is allowed for the export.")

			efun <- NULL
			if (!is.null(opt$efun))
			{
				s <- trimws(opt$efun)
				if (s != "")
				{
					efun <- eval(parse(text=s))
					if (!is.function(efun))
						stop("--efun should be a function.")
				}
			}

			sep <- " "
			if (!is.null(opt$separator))
			{
				sep <- eval(parse(text=paste('"', opt$separator, '"', sep="")))
			}

			dp <- objdesp.gdsn(node)
			if (dp$is.array & !is.null(dp$dim))
			{
				if (!opt$quiet)
					cat(INVERSE("Export:"), " ", opt$export, "\n", sep="")

				fn <- opt$export
				if (fn != "")
				{
					if (last.text(fn, 3L) == ".gz")
						fout <- gzfile(fn, "wb")
					else if (last.text(fn, 3L) == ".xz")
						fout <- xzfile(fn, "wb")
					else
						fout <- file(fn, "wt")
				} else
					fout <- ""

				if (length(dp$dim) == 1L)
				{
					if (is.null(efun))
					{
						apply.gdsn(node, 1L, FUN=cat, as.is="none", file=fout,
							fill=4194304L)
					} else {
						apply.gdsn(node, 1L, FUN=function(x) {
							cat(efun(x), file=fout, fill=4194304L) },
							as.is="none")
					}
				} else if (length(dp$dim) >= 2L)
				{
					if (is.null(efun))
					{
						apply.gdsn(node, length(dp$dim), FUN=cat, as.is="none",
							file=fout, fill=4194304L, sep=sep)
					} else {
						apply.gdsn(node, length(dp$dim), FUN=function(x) {
							cat(efun(x), file=fout, fill=4194304L, sep=sep)
						}, as.is="none")
					}
				}

				if (fn != "") close(fout)
			}
		}

		on.exit()
		closefn.gds(gfile)
		do_clean()

		q(status=0L)
	}


	if (length(files) <= 0L)
	{
		stop("No input file, see \"viewgds -h\".")	
	}

	if (!opt$newfile)
	{
		for (fn in files)
		{
			# open the GDS file
			f <- openfn.gds(fn)
			if (!opt$quiet)
			{
				print(f, all=opt$all, attribute=opt$`show-attr`,
					attribute.trim=FALSE)
			}
			closefn.gds(f)
		}
		if (!opt$quiet) cat("\n")
	}

	# clean up the GDS file(s)
	do_clean()
}


main <- cmpfun(main)
res <- try(main())

# quit
q(status = ifelse(inherits(res, "try-error"), 1L, 0L))
