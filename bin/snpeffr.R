#! /usr/bin/env Rscript

# load docopt
tt <- require(docopt, quietly = TRUE)

if(!tt) {
  print("Installing docopt package first!\n")
  install.packages("docopt", repos = "https://cloud.r-project.org")
  library(docopt, warn.conflicts = FALSE)
}

# configuration for docopt
"
Parse snpeff output to summary of mutations at positions of interest

Usage:
    snpeffr.r [--fpath=FPATH] [--pos=POS] [--genes=GENES] [-exc=EXCL] [-out=OUT]

Options:
    -v, --version           Show version.
    -f FPATH --fpath=FPATH  path to input vcf file from snpeff
    -p POS --pos=POS        format as named comma separated list
                            of positions (no spaces), i.e. see default
                            [default: fks1_hs1=221637:221663,fks1_hs2=223782:223805,fks1_hs3=221805:221807]
    -g GENES --genes=GENES  a list of comma separated gene names (no spaces) [default: CAB11_002014]
    -e EXCL --exc=EXCL      a quoted regular expression for effects to exclude [default: 'synonymous_variant']
    -o OUT --out=OUT        csv or gz file to save output to [default: out.csv]
    -h, --help              show this help text
" -> doc

argus <- docopt(doc, version = 'v1.0\n')

# check arguments
genes <- trimws(unlist(strsplit(argus$genes, split = ",")))
if(length(genes) < 1 || !is.character(genes)) {

  stop("List of genes passed not formatted correctly. These should be
       formatted as: my_gene_1, my_gene_2")

}

pos <- eval(str2lang(paste0("list(", argus$pos, ")")))
if(!is.list(pos) || length(pos) < 1) {
  stop("List of positions passed is not formatted correctly. These should be
       formatted as: gene_name1=posfirst:poslast, gene_name2=posfirst:poslast")
}

# deal with quoting of regular expression
argus$exc <- substr(argus$exc, 2, nchar(argus$exc) - 1)

# run parser
out <- snpeffr::snpeffr(vcf_path = argus$fpath,
                        positions = pos,
                        genes = genes,
                        exclude_effects = argus$exc)

data.table::fwrite(out, argus$out)
