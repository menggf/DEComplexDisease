% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/module.screen.R
\docType{methods}
\name{module.screen}
\alias{module.screen}
\title{Screen the feature patients or genes in predicted modules}
\usage{
module.screen(res.module, feature.patients = NULL, feature.genes = NULL,
  show.mods = NULL, show.n = 4, method = c("ratio", "fisher.test")[1],
  cores = 1)
}
\arguments{
\item{res.module}{a 'seed.module' or 'cluster.module' object returned by \code{\link{seed.module}} or \code{\link{cluster.module}}}

\item{feature.patients}{the patients to screen}

\item{feature.genes}{the genes to screen}

\item{show.mods}{the modules to display}

\item{show.n}{the number of modules to display}

\item{method}{the method to find the most associated modules}

\item{cores}{the thread number}
}
\value{
A plot for gene or patient overlaps with the feature genes or patients.
}
\description{
Screen feature patients or genes given by users among the predicted modules
}
\details{
This function is used to find the modules associated with the 'feature.patients' or 'feature genes'.

In current version, two methods can be used: 'ratio' and 'fisher.test'. `ratio` is to rank the modules based on ratio between observed overlaps and the expected overlaps that estimated using all the samples.  `fisher.test` is to use Fisher's exact test to check the significance of association.
}
\examples{
# screen the modules for feature patients.
module.screen(seed.mod, feature.patients=sample(colnames(deg),15))
}
\author{
Guofeng Meng
}
