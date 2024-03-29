#' Convert DESeq2 objects to a volcano3d object
#'
#' This function is used instead of \code{\link{polar_coords}} if you have raw
#' RNA-Seq count data. It takes 2 `DESeqDataSet` objects, extracts statistical
#' results and converts the results to a 'volc3d' object, which can be directly
#' plotted.
#'
#' @param object An object of class 'DESeqDataSet' with the full design formula.
#'   The function `DESeq` needs to have been run.
#' @param objectLRT An object of class 'DESeqDataSet' with the reduced design
#'   formula. The function `DESeq` needs to have been run on this object with
#'   argument `test="LRT"`.
#' @param contrast Character value specifying column within the metadata stored
#'   in the DESeq2 dataset objects is the outcome variable. This column must
#'   contain a factor with 3 levels. If not set, the function will select the
#'   last term in the design formula of `object` as per DESeq2 convention.
#' @param data Optional matrix containing gene expression data. If not supplied,
#'   the function will pull the expression data from within the DESeq2 object
#'   using the DESeq2 function `assay()`. NOTE: for consistency with gene
#'   expression datasets, genes are in rows.
#' @param pcutoff Cut-off for p-value significance
#' @param padj.method Can be any method available in `p.adjust` or `"qvalue"`.
#'   The option `"none"` is a pass-through.
#' @param filter_pairwise Logical whether adjusted p-value pairwise statistical
#'   tests are only conducted on genes which reach significant adjusted p-value
#'   cut-off on the group likelihood ratio test
#' @param ... Optional arguments passed to \code{\link{polar_coords}}
#' @return Calls \code{\link{polar_coords}} to return an S4 'volc3d' object
#' @seealso \code{\link{polar_coords}}, \code{\link{voom_polar}},
#'   \code{\link[DESeq2:DESeq]{DESeq}} in the DESeq2 package
#' @examples
#' 
#' \donttest{
#'   library(DESeq2)
#' 
#'   counts <- matrix(rnbinom(n=1500, mu=100, size=1/0.5), ncol=15)
#'   cond <- factor(rep(1:3, each=5), labels = c('A', 'B', 'C'))
#' 
#'   # object construction
#'   dds <- DESeqDataSetFromMatrix(counts, DataFrame(cond), ~ cond)
#' 
#'   # standard analysis
#'   dds <- DESeq(dds)
#' 
#'   # Likelihood ratio test
#'   ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
#' 
#'   polar <- deseq_polar(dds, ddsLRT, "cond")
#'   volcano3D(polar)
#'   radial_ggplot(polar)
#' }
#' 
#' @export

deseq_polar <- function(object, objectLRT, contrast = NULL,
                        data = NULL,
                        pcutoff = 0.05,
                        padj.method = "BH",
                        filter_pairwise = TRUE, 
                        ...) {
  
  # Check packages and input data
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Can't find package DESeq2. Try:
           BiocManager::install('DESeq2')", call. = FALSE)
  }
  if (!inherits(object, "DESeqDataSet")) stop("Not a DESeqDataSet object")
  if (!inherits(objectLRT, "DESeqDataSet")) stop("Not a DESeqDataSet object")
  termsDE <- attr(terms(DESeq2::design(object)), "term.labels")
  termsLRT <- attr(terms(DESeq2::design(objectLRT)), "term.labels")
  if (!setequal(termsDE, termsLRT)) message("Different full design formulae")
  if (is.null(contrast)) {
    contrast <- termsDE[length(termsDE)]
    message("Setting contrast to `", contrast, "`")
  }
  outcome <- object@colData[, contrast]
  if (nlevels(outcome) != 3) stop("Outcome does not have 3 levels")
  
  # Group statistics
  LRT <- DESeq2::results(objectLRT)
  LRT <- as.data.frame(LRT[, c('pvalue', 'padj')])
  if (padj.method == "qvalue") {
    LRT <- deseq_qvalue(LRT)
    index <- if (filter_pairwise) {
     LRT$qvalue < pcutoff & !is.na(LRT$qvalue)
    } else !is.na(LRT$qvalue)
    ptype <- "qvalue"
  } else {
    index <- if (filter_pairwise) {
      LRT$padj < pcutoff & !is.na(LRT$padj)
    } else !is.na(LRT$padj)
    ptype <- "padj"
  }
  groups <- levels(outcome)
  contrastlist <- list(groups[1:2], groups[c(3, 1)], groups[2:3])
  
  # Pairwise statistics
  pairres <- lapply(contrastlist, function(i) {
    res <- DESeq2::results(object, contrast = c(contrast, i))
    as.data.frame(res[, c('pvalue', 'padj')])
  })
  pvals <- cbind(LRT[, "pvalue"], pairres[[1]][, "pvalue"], 
                 pairres[[2]][, "pvalue"], pairres[[3]][, "pvalue"])
  pairadj <- lapply(pairres, function(res) {
    out <- rep_len(NA, nrow(LRT))
    out[index] <- if (padj.method == "qvalue") {
      deseq_qvalue(res[index, ])$qvalue
    } else p.adjust(res[index, 'pvalue'], method = padj.method)
    out
  })
  padj <- cbind(LRT[, ptype], pairadj[[1]], pairadj[[2]], pairadj[[3]])
  dimnames(pvals) <- dimnames(padj) <- list(rownames(LRT), 
                                            c("LRT", "AvB", "AvC", "BvC"))
  
  # Extract expression data
  if (is.null(data)) {
    vstdata <- if (nrow(object) < 1000) {
      DESeq2::varianceStabilizingTransformation(object)
    } else DESeq2::vst(object)
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Can't find package SummarizedExperiment. Try:
           BiocManager::install('SummarizedExperiment')",
           call. = FALSE)
    }
    data <- SummarizedExperiment::assay(vstdata)
  }
  out <- polar_coords(outcome, t(data), pvals, padj, pcutoff, ...)
  out@df$type <- "deseq_polar"
  out
}


deseq_qvalue <- function(df, method = "qvalue") {
  q <- qval(df$pvalue[!is.na(df$padj)], method)
  df$qvalue <- 1  # NA converted to 1
  df$qvalue[!is.na(df$padj)] <- q
  df
}

