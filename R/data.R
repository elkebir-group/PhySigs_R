#' Sample set of active mutational signatures.
#'
#' A sample character vector containing subset of active mutational signatures to include.
#'
#' @format A character vector
"S"

#' Sample tumor phylogeny.
#'
#' A sample graphNEL tree where the nodes correspond to clones.
#'
#' @format A graphNEL object with 9 nodes:
#' \describe{
#'   \item{nodes}{clone, with node ID label}
#' }
"T_tree"

#' Sample raw feature matrix
#'
#' A sample dataset containing the trinucleotide counts for each node (clone) in
#' tumor phylogeny T.
#'
#' @format A data frame with 9 rows and 96 variables:
#' \describe{A\[C>A\]A: trinucleotide context count, C to A substitution flanked by A and A}
"P"

#' Sample normalized feature matrix
#'
#' A sample dataset containing the normalized trinucleotide counts for each node (clone) in
#' tumor phylogeny T. These counts were normalized using 'genome' option from
#' deconstructSigs::getTriContextFraction() and then scaled back up to counts.
#'
#' @format A data frame with 9 rows and 96 variables:
#' \describe{E.g., A\[C>A\]A: normalized trinucleotide context count, C to A substitution flanked by A and A}
"P_Norm"

#' Sample list of exposure data frames
#'
#' A sample dataset containing a list mutational signature exposure dataframes for all possible numbers
#' of exposure shifts (i.e., k ranges from 0 to 8).
#'
#' @format A list of data frames, each with 2 rows (active signatures) and x variables (number of clone clusters):
#' \describe{E.g., B;C: Cluster of comprised of nodes B and C}
"E_list"
