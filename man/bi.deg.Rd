% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bi.deg.R
\name{bi.deg}
\alias{bi.deg}
\title{Transformation an expression matrix to binary differential expression matrix}
\usage{
bi.deg(exp, cl, method = c("edger", "deseq2", "normalized")[1],
  cutoff = 0.05, cores = 1)
}
\arguments{
\item{exp}{a matrix or data frame for expression data. The expression value can be counts or normalized expression data}

\item{cl}{a vector of 0 and 1. It has equal length with the column number of exp. 1 indicates the corresponding samples are patients and 0 is control or normal}

\item{method}{defines the methods applied for DE analysis. The possible value is 'edger', 'deseq2', 'normalized'. 'edger' or 'deseq2' is used for RNA-seq count data; 'normalized' is used for normalized RNA-seq or microarray data}

\item{cutoff}{the p-value cutoff for DEGs}

\item{cores}{the thread number}
}
\value{
A deg class object with value of 1, 0 and -1.
}
\description{
Transform the RNA-seq counts or normalized expression matrix into binary differential expression matrix of -1, 0 and 1, which indicates the down-regulation, no change and up-regulation.
}
\details{
For each sample in 'exp', 'cl' defines the patients and normal. The normal samples are used to construct the expression references with negative binomial distribution (e.g. method='edger' or method='deseq2') or a normal distribution (method='normalized').


When counts data are used, the DEG analysis is performed using the functions implemented by `DESeq2` or `edgeR`. The dispersion and mu values are estimated.
}
\examples{
deg <- bi.deg(exp,cl=cl, method='edger', cutoff=0.05) # exp is the RNA-seq counts matrix
}
\author{
Guofeng Meng
}
