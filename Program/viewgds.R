#! /usr/bin/Rscript --vanilla
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("gdsfmt"))

option_list <- list(
	make_option(c("-a", "--all"), action="store_true", default=FALSE,
		help="include hidden GDS node(s)")
)
parser <- OptionParser(usage="%prog [options] file1 [file2] [file3]",
	option_list=option_list)

arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
files <- arguments$args


res <- try({
	for (fn in files)
	{
		gfile <- openfn.gds(fn)
		print(gfile, all=opt$all)
		closefn.gds(gfile)
		cat("\n")
	}
})

q("no", status = if (inherits(res, "try-error")) 1L else 0L)
