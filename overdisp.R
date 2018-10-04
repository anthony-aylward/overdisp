#!/usr/bin/Rscript
#===============================================================================
# overdisp.R
#===============================================================================

"Estimate overdispersion parameters from allele count data

Usage: overdisp [options] [<file>]

-h, --help             show this help message
--min-coverage <int>   minimum coverage level [default: 10]
--nbreaks <int>        number of breaks for the NPBin spline [default: 11]
--ncores <int>         number of cores to use [default: 1]
--order <int>          spline order for NPBin [default: 4]
" -> doc




# Libraries ====================================================================

library(docopt)
library(npbin)

suppressMessages(library(chenimbalance))



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

estimate_null_parameters <- function(
  data_frame,
  minimum_coverage = 10,
  n_breaks = 11,
  spline_order = 4,
  n_cores = 1
) {
  data_table <- convert_to_data_table(
    npbin_preprocess_counts(data_frame),
    minimum_coverage = minimum_coverage
  )

  n <- nrow(data_table)
  breaks <- seq(0, 1, length.out = n_breaks)
  pi_init <- initialize_weights(data_table, n_breaks, spline_order)

  overall_model_estimate <- emBinBspl(
    data_table[, xm],
    data_table[, m],
    breaks = breaks,
    k = spline_order,
    pi.init = pi_init,
    ncores = n_cores,
    err.max = 1e-3,
    iter.max = 200
  )
  null_model_estimate <- estNull(
    data_table[, xm],
    data_table[, m],
    overall_model_estimate,
    init = NULL,
    iter.max = 200,
    ncores = n_cores,
    ub = rep(log(1e4), 2),
    err.max = 1e-4
  )
  null_model_estimate[["coef.null"]]
}

main <- function(opt) {
  counts_frame <- load_data_from_stdin_or_file(opt[["file"]])
  counts_frame <- counts_frame[["coverage"]] >= opt[["min-coverage"]]
  total <- counts_frame[["coverage"]]
  allelic_ratio <- counts_frame[["ref_count"]] / total
  lsse_parameters <- alleledb_beta_binomial(
    total[!is.na(total)],
    allelic_ratio[!is.na(total)],
    r_by = 0.025
  )
  null_parameters <- estimate_null_parameters(
    counts_frame,
    minimum_coverage = opt[["min-coverage"]],
    n_breaks = opt[["nbreaks"]],
    spline_order = opt[["order"]],
    n_cores = opt[["ncores"]]
  )
  print(null_parameters)
}

parse_options <- function(doc) {
  opt <- docopt(doc)
  for (option in c("min-coverage", "nbreaks", "order", "ncores")) {
    opt[[option]] <- as.integer(opt[[option]])
  }
  opt
}




# Execute ======================================================================

opt <- parse_options(doc)
main(opt)
