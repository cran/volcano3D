% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deseq_polar.R
\name{deseq_polar}
\alias{deseq_polar}
\title{Convert DESeq2 objects to a volcano3d object}
\usage{
deseq_polar(
  object,
  objectLRT,
  contrast = NULL,
  data = NULL,
  pcutoff = 0.05,
  padj.method = "BH",
  filter_pairwise = TRUE,
  ...
)
}
\arguments{
\item{object}{An object of class 'DESeqDataSet' with the full design formula.
The function \code{DESeq} needs to have been run.}

\item{objectLRT}{An object of class 'DESeqDataSet' with the reduced design
formula. The function \code{DESeq} needs to have been run on this object with
argument \code{test="LRT"}.}

\item{contrast}{Character value specifying column within the metadata stored
in the DESeq2 dataset objects is the outcome variable. This column must
contain a factor with 3 levels. If not set, the function will select the
last term in the design formula of \code{object} as per DESeq2 convention.}

\item{data}{Optional matrix containing gene expression data. If not supplied,
the function will pull the expression data from within the DESeq2 object
using the DESeq2 function \code{assay()}. NOTE: for consistency with gene
expression datasets, genes are in rows.}

\item{pcutoff}{Cut-off for p-value significance}

\item{padj.method}{Can be any method available in \code{p.adjust} or \code{"qvalue"}.
The option \code{"none"} is a pass-through.}

\item{filter_pairwise}{Logical whether adjusted p-value pairwise statistical
tests are only conducted on genes which reach significant adjusted p-value
cut-off on the group likelihood ratio test}

\item{...}{Optional arguments passed to \code{\link{polar_coords}}}
}
\value{
Calls \code{\link{polar_coords}} to return an S4 'volc3d' object
}
\description{
This function is used instead of \code{\link{polar_coords}} if you have raw
RNA-Seq count data. It takes 2 \code{DESeqDataSet} objects, extracts statistical
results and converts the results to a 'volc3d' object, which can be directly
plotted.
}
\examples{

\donttest{
  library(DESeq2)

  counts <- matrix(rnbinom(n=1500, mu=100, size=1/0.5), ncol=15)
  cond <- factor(rep(1:3, each=5), labels = c('A', 'B', 'C'))

  # object construction
  dds <- DESeqDataSetFromMatrix(counts, DataFrame(cond), ~ cond)

  # standard analysis
  dds <- DESeq(dds)

  # Likelihood ratio test
  ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)

  polar <- deseq_polar(dds, ddsLRT, "cond")
  volcano3D(polar)
  radial_ggplot(polar)
}

}
\seealso{
\code{\link{polar_coords}}, \code{\link{voom_polar}},
\code{\link[DESeq2:DESeq]{DESeq}} in the DESeq2 package
}
