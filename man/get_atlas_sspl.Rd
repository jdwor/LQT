% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_atlas_sspl.R
\name{get_atlas_sspl}
\alias{get_atlas_sspl}
\title{Compute the structural shortest path length (SSPL) for the SC atlas}
\usage{
get_atlas_sspl(cfg)
}
\arguments{
\item{cfg}{a pre-made cfg structure (as list object).}
}
\value{
An .RData file with the suffix _sspl_matrix.RData, which contains the raw SSPL
matrix computed from the binarized atlas SC matrix.
}
\description{
This function uses the atlas SC matrix to compute an atlas SSPL matrix
containing SSPLs between each pair of brain regions. It assumes that you have already
run get_parcel_atlas to obtain the atlas SC matrix.
}
