% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{P_Norm}
\alias{P_Norm}
\title{Sample normalized feature matrix}
\format{
A data frame with 9 rows and 96 variables:
\describe{E.g., A[C>A]A: normalized trinucleotide context count, C to A substitution flanked by A and A}
}
\usage{
P_Norm
}
\description{
A sample dataset containing the normalized trinucleotide counts for each node (clone) in
tumor phylogeny T. These counts were normalized using 'genome' option from
deconstructSigs::getTriContextFraction() and then scaled back up to counts.
}
\keyword{datasets}
