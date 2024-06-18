## ----message=FALSE, warning=FALSE---------------------------------------------
library(beeca)
library(dplyr)

## -----------------------------------------------------------------------------
## Prepare the dataset for input
dat <- trial02_cdisc %>%
  ## In this case we are only interested in comparing two arms: Placebo vs Xanomeline High Dose
  dplyr::filter(TRTP %in% c("Placebo", "Xanomeline High Dose")) %>%
  ## Treatment variable must be coded as a factor
  dplyr::mutate(TRTP = factor(TRTP))

## Fit the logistic regression model adjusting for SEX, RACE and AGE
fit <- glm(AVAL ~ TRTP + SEX + RACE + AGE, family = "binomial", data = dat)

## Calculate the marginal treatment effect estimate and associated variance for a difference contrast
## using the Ge et al. method with robust HC0 sandwich variance
marginal_fit <- get_marginal_effect(fit,
  trt = "TRTP",
  method = "Ge",
  type = "HC0",
  contrast = "diff",
  reference = "Placebo"
)

## -----------------------------------------------------------------------------
## View the ARD summary
marginal_fit$marginal_results

## -----------------------------------------------------------------------------
## Extract results
marginal_results <- marginal_fit$marginal_results
diff_est <- marginal_results[marginal_results$STAT == "diff", "STATVAL"][[1]]
diff_se <- marginal_results[marginal_results$STAT == "diff_se", "STATVAL"][[1]]

## 95% confidence interval
ci_l <- diff_est - (qnorm(0.975) * diff_se)
ci_u <- diff_est + (qnorm(0.975) * diff_se)

## Two-sided p-value
z_score <- diff_est / diff_se
p_value <- 2 * (1 - pnorm(abs(z_score)))

sprintf("The risk difference is %s with 95%% CI: (%s - %s)", round(diff_est, 2), round(ci_l, 2), round(ci_u, 2))
sprintf("p-value: %s", formatC(p_value, format = "e", digits = 2))

## -----------------------------------------------------------------------------
# pre-process trial01 dataset to convert treatment arm to a factor and handle missing value
data01 <- trial01 |>
  transform(trtp = as.factor(trtp)) |>
  dplyr::filter(!is.na(aval))
fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = data01)
beeca_ge <- get_marginal_effect(object = fit, trt = "trtp", method = "Ge", 
                                contrast = "diff", reference = "0", type = "model-based")
cat("Point estimate", beeca_ge$marginal_est, "\nStandard error estimate", beeca_ge$marginal_se)

## -----------------------------------------------------------------------------
ge_var_paper <- function(glmfit, trt) {
  pder <- function(ahat, vc, x) {
    #### ahat: logistic regression parameters
    #### vc: variance-covariance matrix of ahat
    #### x: full model matrix of the logistic regression
    #### return mean of phat on x and its se
    phat <- plogis(x %*% ahat)
    pbar <- mean(phat)
    pderiv <- t(phat * (1 - phat)) %*% x / nrow(x)
    sepbar <- sqrt(pderiv %*% vc %*% t(pderiv))
    return(list(pbar = pbar, sepbar = sepbar, pderiv = pderiv))
  }

  difP <- function(glmfit) {
    #### estimate the proportion difference and its standard error
    df <- glmfit$model

    vc <- vcov(glmfit)
    df[, trt] <- 1
    mat <- model.matrix(glmfit$formula, data = df)
    pderT <- pder(coef(glmfit), vc, mat)
    df[, trt] <- 0
    mat <- model.matrix(glmfit$formula, data = df)
    pderC <- pder(coef(glmfit), vc, mat)

    difb <- pderT$pbar - pderC$pbar
    sedif <- sqrt((pderT$pderiv - pderC$pderiv) %*% vc %*% t(pderT$pderiv - pderC$pderiv))

    return(list(
      pT = pderT$pbar,
      pC = pderC$pbar,
      dif = difb,
      sedif = sedif,
      var = sedif**2
    ))
  }

  return(list(est = difP(glmfit)$dif, se = difP(glmfit)$sedif[[1]]))
}
paper_ge <- ge_var_paper(fit, "trtp")
cat("Point estimate", paper_ge$est, "\nStandard error estimate", paper_ge$se)

## -----------------------------------------------------------------------------
if (requireNamespace("margins", quietly = T)) {
  margins_ge <- margins::margins(model = fit, variables = "trtp", vcov = vcov(fit))
  cat("Point estimate", summary(margins_ge)$AME, "\nStandard error estimate", summary(margins_ge)$SE)
}

## -----------------------------------------------------------------------------
if (requireNamespace("marginaleffects", quietly = T)) {
  marginaleffects_ge <- marginaleffects::avg_comparisons(fit, variables = "trtp")
  cat("Point estimate", marginaleffects_ge$estimate, "\nStandard error estimate", marginaleffects_ge$std.error)
}

## -----------------------------------------------------------------------------
cat("Point estimate", margins_trial01$Estimate, "\nStandard error estimate", margins_trial01$StdErr)

## -----------------------------------------------------------------------------
beeca_ye <- get_marginal_effect(object = fit, trt = "trtp", method = "Ye", 
                                contrast = "diff", reference = "0")
cat("Point estimate", beeca_ye$marginal_est, "\nStandard error estimate", beeca_ye$marginal_se)

## -----------------------------------------------------------------------------
if (requireNamespace("RobinCar", versionCheck = list(name = "RobinCar", op = "==", version = "0.3.0"), quietly = T)) {
  robincar_ye <- RobinCar::robincar_glm(data.frame(fit$data), response_col = as.character(fit$formula[2]),
      treat_col = "trtp", formula = fit$formula, g_family = fit$family,
      contrast_h = "diff")$contrast$result
  cat("Point estimate", robincar_ye$estimate, "\nStandard error estimate", robincar_ye$se)
}

