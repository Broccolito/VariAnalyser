% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VariantAnnotator.R
\name{VariantAnnotator}
\alias{VariantAnnotator}
\title{Annotate the GWAS results}
\usage{
VariantAnnotator(
  variant_file = "data/example/FHS_EA_MRS_5e8_snplist.txt",
  chr = "CHR",
  pos = "POS",
  ref = "NEA",
  alt = "EA",
  delimiter = " ",
  write_csv = TRUE
)
}
\arguments{
\item{variant_file}{The name of the text file containing GWAS statistics}

\item{chr}{The name of the chromosome number column in the text file}

\item{pos}{The name of the chromosome position column in the text file}

\item{ref}{The name of the reference allele column in the text file}

\item{alt}{The name of the alternative allele column in the text file}

\item{delimiter}{The delimiter used in the text file}

\item{write_csv}{Whether to write the result to a csv file, else skip saving}

\item{pvalue}{The name of the pvalue column in the text file}
}
\value{
A data frame of the annotated results
}
\description{
Annotate the GWAS results
}
\examples{
VariantAnnotator("data/example/FHS_EA_MRS_merged.txt")
}
