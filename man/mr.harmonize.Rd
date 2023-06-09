% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mr.harmonize.R
\name{mr.harmonize}
\alias{mr.harmonize}
\title{Perform mr harmonization}
\usage{
mr.harmonize(By, Bx, By.se, Bx.se, outcome.pval, exposure.pval,
  outcome.effect_allele, exposure.effect_allele, exposure.other_allele,
  outcome.SNP, exposure.SNP)
}
\arguments{
\item{By}{Vector of genetic effects  from outcome}

\item{Bx}{Vector of genetic effects from exposure}

\item{By.se}{Standard errors of genetic effects from outcome}

\item{Bx.se}{Standard errors of genetic effects from exposure}

\item{outcome.pval}{p-value from outcome}

\item{exposure.pval}{p-value from exposure}

\item{outcome.effect_allele}{vector containing the effect_alleles from outcome}

\item{exposure.effect_allele}{vector containing the effect_alleles from exposure}

\item{exposure.other_allele}{vector containing the effect_alleles from exposure}

\item{outcome.SNP}{vector containing rs_ids from outcome}

\item{exposure.SNP}{vector containing SNPs from exposure}

\item{outcome.other_allele}{vector containing the effect_alleles from outcome}
}
\value{
harmonized data
}
\description{
allelic harmonization between outcome and exposure. The intersect is taken from both outcome and exposure. The alleles and effect sizes are aligned. Mismatch alles are removed. And duplicates are removed.
}
\examples{
mr.harmonize()

}
\author{
Matthijs Knigge
}
\keyword{harmonize}
