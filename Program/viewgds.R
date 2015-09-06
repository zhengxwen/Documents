#! /usr/bin/Rscript --vanilla
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
	make_option(c("-e", "--export"), type="character",
		help="Export a GDS node to a text file (-e \"\" for standard output)",
		metavar="filename"),
	make_option(c("-i", "--import"), type="character",
		help="Import a text file and create a GDS node \"-c\" (-i \"\" for standard input)",
		metavar="filename"),
	make_option(c("--icol"), type="character",
		help="Specify the columns for importing data, allowing multiple variables splited by semicolon\n\t\tE.g., --icol \"-(1:4)\" for excluding the first four columns",
		metavar="column"),
	make_option(c("-c", "--create"), action="store", type="character",
		help="Create a GDS node with the format TYPE:DIM:COMPRESSION[;TYPE2:DIM2:COMPRESSION2;...]\n\t\tE.g., -n NAME -c \"int:4,0:ZIP_RA.max\"",
		metavar="format"),
	make_option("--delete", action="store_true", default=FALSE,
		help="Delete the GDS node (e.g., -n NAME --delete)"),
	make_option("--newfile", action="store_true", default=FALSE,
		help="Create a new GDS file"),
	make_option("--show-attr", action="store_true", default=FALSE,
		help="Show the attribute(s)"),
	make_option("--system", action="store_true", default=FALSE,
		help="Show the system configuration"),
	make_option("--quiet", action="store_true", default=FALSE,
		help="No screen output"),
	make_option("--clean", action="store_true", default=FALSE,
		help="Clean up the fragments of GDS file")
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


# read and drop upper dimensions
read <- function(node, start, count)
{
	v <- read.gdsn(node, start=start, count=count, simplify="none")
	if (!is.null(dm <- dim(v)))
	{
		if (length(dm) > 2L)
			dim(v) <- dm[c(1L,2L)]
	}
	v
}


# view 1-dim array
view.dim1 <- function(dm, node)
{
	if (dm <= 0L)
	{
		s <- ""
	} else if (dm <= nprev*2L)
	{
		s <- format(read.gdsn(node))
	} else {
		s <- format(c(read(node, 1L, nprev), read(node, dm-nprev+1L, nprev)))
		s <- s[c(1:nprev, NA, seq(nprev+1L, length(s)))]
		s[nprev+1L] <- BLURRED("...")
	}
	cat(s, sep="\n")
}

# view 2-dim array
view.dim2 <- function(dm, node, st=NULL)
{
	if (any(dm[1L] <= 0L, dm[2L] <= 0L))
	{
		cat("\n")
		return(invisible())
	}

	cn <- rep(1L, length(st))
	if (dm[1L] <= nprev*2L)
	{
		if (dm[2L] <= nprev*2L)
		{
			v <- read(node, c(1L,1L,st), c(-1L,-1L,cn))
		} else {
			v <- cbind(
				read(node, c(1L,1L,st), c(-1L,nprev,cn)),
				read(node, c(1L,dm[2L]-nprev+1L,st), c(-1L,nprev,cn))
			)
		}
	} else {
		if (dm[2L] <= nprev*2L)
		{
			v <- rbind(
				read(node, c(1L,1L,st), c(nprev,-1L,cn)),
				read(node, c(dm[1L]-nprev+1L,1L,st), c(nprev,-1L,cn))
			)
		} else {
			v1 <- cbind(
				read(node, c(1L,1L,st), c(nprev,nprev,cn)),
				read(node, c(1L,dm[2L]-nprev+1L,st), c(nprev,nprev,cn))
			)
			v2 <- cbind(
				read(node, c(dm[1L]-nprev+1L,1L,st), c(nprev,nprev,cn)),
				read(node, c(dm[1L]-nprev+1L,dm[2L]-nprev+1L,st), c(nprev,nprev,cn))
			)
			v <- rbind(v1, v2)
		}
	}

	s <- format(v)
	if (dm[2L] > nprev*2L)
	{
		s <- s[, c(1:nprev, NA, seq(nprev+1L,ncol(s))), drop=FALSE]
		s[, nprev+1L] <- BLURRED("..")
	}
	if (dm[1L] > nprev*2L)
	{
		s <- s[c(1:nprev, NA, seq(nprev+1L,nrow(s))), , drop=FALSE]
		s[nprev+1L, ] <- ""
		s[nprev+1L, 1L] <- BLURRED("......")
	}

	write.table(s, col.names=FALSE, row.names=FALSE, quote=FALSE)
	invisible()
}

# view >2-dim array
view.dim <- function(i, st, dm, node)
{
	if (i > length(dm))
	{
		cat(UNDERLINE(sprintf("[,,%s]:\n", paste(st, collapse=","))))
		view.dim2(dm, node, st)
	} else {
		for (j in seq_len(min(dm[i], nprev)))
		{
			st2 <- c(st, j)
			view.dim(i + 1L, st2, dm, node)
		}
	}
}


last.text <- function(s, len)
{
	substring(s, nchar(s)-len+1L, nchar(s))
}



####  run the main part  ####

main <- function()
{
	if (opt$system)
	{
		v <- c(package.version = as.character(packageVersion("gdsfmt")),
			system.gds())
		print(v)
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
		readonly <- is.null(opt$create) & !opt$delete
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
					valdim=dm, compress=fmt[3L])
				args <- c(args, eval(parse(text=paste("list(", fmt[4L], ")"))))

				n <- do.call(add.gdsn, args)
				if (!opt$quiet & is.null(opt$import))
					print(n, attribute=TRUE, attribute.trim=FALSE)
				nodelist[[i]] <- n
			}

			# import
			if (!is.null(opt$import))
			{
				ic <- vector("list", length(nodenames))
				if (!is.null(opt$icol))
				{
					s <- trimws(unlist(strsplit(opt$icol, ";", fixed=TRUE)))
					if (length(s) != length(nodenames))
						stop("'--icol' should have the same number of elements as '-n'.")
					for (i in seq_len((length(s))))
					{
						if (s[i] != "")
							ic[[i]] <- eval(parse(text=paste("c(", s[i], ")")))
					}
				}

				if (opt$import == "") opt$import <- "stdin"
				f <- file(opt$import)
				open(f)
				while(length(s <- readLines(f, n=1L)) > 0L)
				{
					ss <- scan(text=s, what="", quiet=TRUE, strip.white=TRUE)
					for (i in seq_len(length(nodelist)))
					{
						v <- ic[[i]]
						if (is.null(v))
							append.gdsn(nodelist[[i]], ss)
						else
							append.gdsn(nodelist[[i]], ss[v])
					}
				}
				close(f)

				for (i in seq_len(length(nodelist)))
					readmode.gdsn(nodelist[[i]])

				if (!opt$quiet)
				{
					for (i in seq_len(length(nodelist)))
					{
						cat(INVERSE("New GDS node:"), "")
						print(nodelist[[i]], attribute=TRUE, attribute.trim=FALSE)
					}
				}
			}

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
			
			q(status=0L)
		}

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
					if (length(dp$dim) == 1L)
					{
						view.dim1(dp$dim, node)
					} else if (length(dp$dim) == 2L)
					{
						view.dim2(dp$dim, node)
					} else {
						view.dim(3L, NULL, dp$dim, node)
					}
				}
			}
		}

		if (!is.null(opt$export))
		{
			if (length(nodenames) > 1L)
				stop("Only a GDS node is allowed for the export.")

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
					apply.gdsn(node, 1L, FUN=cat, as.is="none", file=fout, fill=TRUE)
				} else if (length(dp$dim) == 2L)
				{
					apply.gdsn(node, 1L, FUN=cat, as.is="none", file=fout, fill=65536L)
				}

				if (fn != "") close(fout)
			}
		}

		q(status=0L)
	}


	if (length(files) <= 0L)
	{
		stop("No input file, see \"viewgds -h\".")	
	}

	for (fn in files)
	{
		if (opt$newfile)
		{
			f <- createfn.gds(fn)
			if (!opt$quiet)
			{
				cat(INVERSE("New file:"), "")
				cat(f$filename)
			}
			closefn.gds(f)
		} else {
			# open the GDS file
			f <- openfn.gds(fn)
			if (!opt$quiet)
			{
				print(f, all=opt$all, attribute=opt$`show-attr`,
					attribute.trim=FALSE)
			}
			closefn.gds(f)

			# clean up the GDS file
			if (opt$clean)
				cleanup.gds(fn, verbose=!opt$quiet)
		}

		if (!opt$quiet) cat("\n")
	}
}


main <- cmpfun(main)
res <- try(main())
v <- warnings()
if (!is.null(v)) cat(v, sep="\n")

# quit
q(status = ifelse(inherits(res, "try-error"), 1L, 0L))
