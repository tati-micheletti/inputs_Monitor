#' Select non-collinear predictors ranked by block-CV importance
#'
#' For each predictor, computes block-CV explained deviance via
#' `computeUnivarCV()`, then greedily drops the less-important member of
#' any pair correlated above `threshold`. Follows Wiedenroth et al.
#' (`select07_blockcv`).
#'
#' @param X data.frame or matrix of candidate predictor columns.
#' @param y Numeric vector of observed 0/1 responses.
#' @param family Character. GLM family, default "binomial".
#' @param univar Character, initial univariate model form (auto-refined
#'   internally by `computeUnivarCV()`).
#' @param threshold Numeric. Absolute correlation threshold, default 0.7.
#' @param method Character. Correlation method, default "spearman".
#' @param sequence Character vector, optional fixed variable ordering
#'   (default: importance-ranked).
#' @param weights Numeric vector of observation weights.
#' @param spBlock A blocks object (from `spatialBlockingX()` or
#'   `mimicSpatialBlocks()`) with a `$folds_ids` element.
#' @return List with `D2` (importance scores), `cor_mat`, and `pred_sel`
#'   (selected, non-collinear predictor names in importance order).
select07Blockcv <- function(X, y, family = "binomial", univar = "glm2",
                             threshold = 0.7, method = "spearman",
                             sequence = NULL, weights = NULL, spBlock) {
  ks <- spBlock$folds_ids
  imp <- apply(X, 2, computeUnivarCV, response = y, family = family,
               univar = univar, ks = ks, weights = weights)

  cm <- cor(X, method = method, use = "complete.obs")

  sortImp <- if (is.null(sequence)) {
    colnames(X)[order(imp, decreasing = TRUE)]
  } else {
    sequence
  }

  pairs <- which(abs(cm) >= threshold, arr.ind = TRUE)
  index <- which(pairs[, 1] == pairs[, 2])
  pairs <- pairs[-index, ]

  exclude <- NULL
  for (i in seq_along(sortImp)) {
    if ((sortImp[i] %in% row.names(pairs)) & ((sortImp[i] %in% exclude) == FALSE)) {
      cv <- cm[setdiff(row.names(cm), exclude), sortImp[i]]
      cv <- cv[setdiff(names(cv), sortImp[1:i])]
      exclude <- c(exclude, names(which(abs(cv) >= threshold)))
    }
  }

  predSel <- sortImp[!(sortImp %in% unique(exclude))]
  list(D2 = sort(imp, decreasing = TRUE), cor_mat = cm, pred_sel = predSel)
}
