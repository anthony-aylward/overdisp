#!/usr/bin/Rscript
#===============================================================================
# overdisp.R
#===============================================================================

"Estimate overdispersion parameters from allele count data

Usage: overdisp.R [<file>]

" -> doc




# Libraries ====================================================================

library(npbin)

suppressMessages(library(chenimbalance))
suppressMessages(library(docopt))




# Functions ====================================================================

load_data_from_stdin_or_file <- function(file_option) {
  if (is.null(file_option)) {
    f <- file("stdin", "r")
    data_frame <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
    close(f)
    data_frame
  } else {
    read.table(file_option, header = TRUE, stringsAsFactors = FALSE)
  }
}

main <- function(opt) {
  counts <- load_data_from_stdin_or_file(opt[["file"]])
  print(head(counts))
}




# Execute ======================================================================

opt <- docopt(doc)
main(opt)
