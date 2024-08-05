#' Calculate the FGFR3 Signature Score
#'
#' @name FGFR3_signature
#' @param expr A matrix or data frame of gene expression values with genes as rows and samples as columns.
#' @return A data frame with FGFR3 signature scores for each sample.
#' @export
FGFR3_signature <- function(expr) {
  if (!is.matrix(expr) && !is.data.frame(expr)) {
    stop("The 'expr' argument must be a matrix or data frame.")
  }
  expr <- as.matrix(expr)
  if (any(expr < 0)) {
    expr <- expr + abs(min(expr(expr < 0)))
  }
  max_expr <- max(expr)
  min_expr <- min(expr)
  if (max_expr - min_expr > 30) {
    expr <- log2(expr + 1)
  }
  library(GSVA)
  load(system.file("data/FGFR3_signature_gene.rdata", package = "FGFR3signature"))
  cellMarker <- FGFR3_signature_gene
  cellMarker <- lapply(split(cellMarker,cellMarker$Term), function(x){
  dd = x$Gene
  unique(dd)
})
  gsva_data <- as.data.frame(t(gsva(expr,cellMarker, method = "ssgsea")))
  gsva_data$FGFR3_signature <- gsva_data$FGFR3_UP - gsva_data$FGFR3_DN
  FGFR3_signature <- gsva_data
  return(FGFR3_signature)
}