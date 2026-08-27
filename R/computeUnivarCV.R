#' Univariate block-CV importance score for one predictor variable
#'
#' Fits a univariate GLM/GAM per fold (choosing glm1/glm2/gam based on
#' the number of unique values seen in the training fold), predicts on
#' the held-out fold, and returns the explained deviance across all
#' folds. Used to rank predictors by importance before collinearity
#' filtering. Follows Wiedenroth et al.
#'
#' @param variable Numeric vector, one predictor's values.
#' @param response Numeric vector of observed 0/1 responses.
#' @param family Character. GLM family, default "binomial".
#' @param univar Character, unused as an override -- kept for signature
#'   compatibility; the univariate model form is chosen automatically.
#' @param ks Integer vector of fold assignments (same length as `response`).
#' @param weights Numeric vector of observation weights.
#' @return Numeric, explained deviance (D2) for this variable.
computeUnivarCV <- function(variable, response, family, univar, ks, weights) {
  preds <- numeric(length(response))

  univarSelect <- function(m) {
    df <- data.frame(occ = response, env = variable)
    trainDf <- df[!ks == m, ]
    if (length(unique(trainDf[, 2])) > 4) return("gam")
    if (length(unique(trainDf[, 2])) < 5 & length(unique(trainDf[, 2])) > 2) return("glm2")
    return("glm1")
  }

  univarList <- lapply(sort(unique(ks)), univarSelect)
  if (any(univarList == "glm1")) univar <- "glm1"
  else if (any(univarList == "glm2")) univar <- "glm2"
  else univar <- "gam"

  for (n in sort(unique(ks))) {
    df <- data.frame(occ = response, env = variable)
    trainDf <- df[!ks == n, ]
    testDf <- df[ks == n, ]
    m1 <- switch(univar,
                 glm1 = glm(occ ~ env, data = trainDf, family = family, weights = weights[!ks == n]),
                 glm2 = glm(occ ~ poly(env, 2), data = trainDf, family = family, weights = weights[!ks == n]),
                 gam = mgcv::gam(occ ~ s(env, k = 4), data = trainDf, family = family, weights = weights[!ks == n]))
    preds[ks == n] <- predict(m1, newdata = testDf, type = "response")
  }
  d2 <- explDeviance(response, preds, weights = weights)
  ifelse(d2 < 0, 0, d2)
}
