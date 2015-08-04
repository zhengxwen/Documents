#! /usr/bin/Rscript --vanilla
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("gdsfmt"))
crayon.flag <- requireNamespace("crayon", quietly=TRUE)

option_list <- list(
	make_option(c("-a", "--all"), action="store_true", default=FALSE,
		help="Include hidden GDS node(s)"),
	make_option(c("-n", "--node"), action="store", type="character",
		help="Specify a GDS node"),
	make_option("--quiet", action="store_true", default=FALSE,
		help="No screen output")
)
parser <- OptionParser(usage="%prog [options] file1 [file2] [file3]",
	option_list=option_list)

arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
files <- arguments$args


# default
nprev <- 6L
BLURRED <- if (crayon.flag) crayon::blurred else `c`
INVERSE <- if (crayon.flag) crayon::inverse else `c`
UNDERLINE<- if (crayon.flag) crayon::underline else `c`


# view 1-dim array
view.dim1 <- function(dm, node)
{
	if (dm <= nprev*2L)
	{
		format(read.gdsn(node))
	} else {
		c(
			format(read.gdsn(node, start=1L, count=nprev)),
			BLURRED("..."),
			format(read.gdsn(node, start=dm-nprev+1L, count=nprev))
		)
	}
}

# view 2-dim array
view.dim2 <- function(dm, node, st=NULL)
{
	cn <- rep(1L, length(st))

	if (dm[1L] <= nprev*2L)
	{
		if (dm[2L] <= nprev*2L)
		{
			v <- read.gdsn(node, start=c(1L,1L,st), count=c(-1L,-1L,cn))
		} else {
			v <- cbind(
				read.gdsn(node, start=c(1L,1L,st), count=c(-1L,nprev,cn)),
				read.gdsn(node, start=c(1L,dm[2L]-nprev+1L,st), count=c(-1L,nprev,cn))
			)
		}
	} else {
		if (dm[2L] <= nprev*2L)
		{
			v <- cbind(
				read.gdsn(node, start=c(1L,1L,st), count=c(nprev,-1L,cn)),
				read.gdsn(node, start=c(dm[1L]-nprev+1L,1L,st), count=c(nprev,-1L,cn))
			)
		} else {
			v1 <- cbind(
				read.gdsn(node, start=c(1L,1L,st), count=c(nprev,nprev,cn)),
				read.gdsn(node, start=c(1L,dm[2L]-nprev+1L,st), count=c(nprev,nprev,cn))
			)
			v2 <- cbind(
				read.gdsn(node, start=c(dm[1L]-nprev+1L,1L,st),
					count=c(nprev,nprev,cn)),
				read.gdsn(node, start=c(dm[1L]-nprev+1L,dm[2L]-nprev+1L,st),
					count=c(nprev,nprev,cn))
			)
			v <- rbind(v1, v2)
		}
	}

	s <- format(v)
	if (dm[2L] > nprev*2L)
	{
		s <- s[, c(1:nprev, NA, seq(nprev+1L,ncol(s)))]
		s[, nprev+1L] <- BLURRED("..")
	}
	if (dm[1L] > nprev*2L)
	{
		s <- s[c(1:nprev, NA, seq(nprev+1L,nrow(s))), ]
		s[nprev+1L, ] <- ""
		s[nprev+1L, 1L] <- BLURRED("......")
	}

	write.table(s, col.names=FALSE, row.names=FALSE, quote=FALSE)	
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


# run
res <- try({
	if (!is.null(opt$node))
	{
		if (length(files) != 1L)
			stop("Please specify only one GDS file.")

		gfile <- openfn.gds(files)

		node <- index.gdsn(gfile, opt$node)
		if (!opt$quiet)
		{
			cat(INVERSE("GDS node:"), "")
			print(node, all=opt$all)

			dp <- objdesp.gdsn(node)
			if (dp$is.array & !is.null(dp$dim))
			{
				cat(INVERSE("Preview:\n"))
				if (length(dp$dim) == 1L)
				{
					cat(view.dim1(dp$dim, node), sep="\n")
				} else if (length(dp$dim) == 2L)
				{
					view.dim2(dp$dim, node)
				} else {
					view.dim(3L, NULL, dp$dim, node)
				}
			}
		}

		closefn.gds(gfile)
		q(status=0L)
	}


	for (fn in files)
	{
		gfile <- openfn.gds(fn)
		if (!opt$quiet)
			print(gfile, all=opt$all)
		closefn.gds(gfile)
		cat("\n")
	}
})

# quit
q(status = if (inherits(res, "try-error")) 1L else 0L)
