#' @title Three-step distal outcome analysis
#'
#' @description
#' Regress distal outcome on attribute classifications or latent
#' profile classifications from a fitted \code{\link{GDINA}} object using a
#' three-step correction. The corrected analysis uses the estimated
#' misclassification matrix \eqn{P(W=s\mid X=t)} from \code{\link{CM}} and can
#' be fit by maximum likelihood (ML) and/or the BCH weighting approach.
#'
#' @param object An estimated \code{GDINA} object returned from
#'   \code{\link{GDINA}}.
#' @param outcome A distal outcome with one value per respondent. Supported
#'   outcomes are binary, ordinal, nominal categorical, and continuous
#'   outcomes.
#' @param formula Optional one-sided or two-sided formula specifying observed
#'   covariates to include additively in the distal outcome regression. Only
#'   the right-hand side is used.
#' @param data Optional data frame containing the covariates in
#'   \code{formula}. It must contain one row per respondent in \code{object}.
#'   Ordered factors are converted to numeric scores in level order so they
#'   enter the regression with a single linear effect.
#' @param level A character string specifying whether the corrected analysis is
#'   carried out at the attribute level (\code{"attribute"}) or the latent
#'   profile level (\code{"profile"}).
#' @param attribute Optional integer vector giving the attributes to include in
#'   the attribute-level regression when \code{level = "attribute"}. Selected
#'   attributes enter the same regression model together. If \code{NULL}, all
#'   attributes are included.
#' @param classification A character string specifying the classification rule.
#'   Supported values are \code{"MAP"}, \code{"MLE"}, and \code{"EAP"}.
#'   Alternatively, a matrix of user-supplied classifications can be provided,
#'   with one row per respondent and one column per attribute.
#' @param method A character vector specifying which corrected distal analyses
#'   to return. Supported values are \code{"ML"} and \code{"BCH"}.
#' @param outcome_type A character string specifying the outcome type.
#'   Supported values are \code{"auto"}, \code{"binary"},
#'   \code{"ordinal"}, \code{"categorical"}, \code{"nominal"}, and
#'   \code{"continuous"}. The value \code{"nominal"} is treated as an
#'   alias for \code{"categorical"}.
#' @param reference A character string specifying the reference class for
#'   profile-level models. Supported values are \code{"last"} (default) and
#'   \code{"first"}. Ignored when \code{level = "attribute"}.
#' @param profile_group Optional vector or factor used to combine latent
#'   profiles when \code{level = "profile"}. Supply one value per profile in
#'   the order returned by \code{extract(object, "attributepattern")};
#'   profiles sharing the same value are collapsed into one grouped state.
#'   A named vector is also accepted and will be matched to profile labels
#'   such as \code{"000"} or \code{"101"}. If \code{NULL}, each profile is
#'   analyzed separately.
#' @param conf.level Confidence level used for Wald intervals.
#' @param maxit Maximum number of optimization iterations.
#' @param na.action See \code{\link{na.action}} for handling missing values in
#'   observed covariates.
#'
#' @details
#' Let \eqn{W} denote the estimated latent class, \eqn{X} the unobserved true
#' class, \eqn{Y} the distal outcome, and \eqn{Z} the observed covariates. The
#' ML correction is based on the
#' mixture likelihood
#' \deqn{\sum_i \log\left\{\sum_t P(W_i=s_i\mid X_i=t)P(X_i=t)f(Y_i\mid X_i=t, Z_i)\right\}.}
#' Observed covariates enter additively in the distal outcome model together
#' with the selected latent attributes or profiles. For binary distal outcomes,
#' \eqn{f(Y_i\mid X_i=t, Z_i)} is Bernoulli with a logit link; for continuous
#' outcomes it is Gaussian with state-specific means plus additive covariate
#' effects and a common residual standard deviation; for ordinal outcomes it is
#' a cumulative logit model with proportional odds; for nominal categorical
#' outcomes it is a multinomial logit model over the distal categories.
#'
#' When \code{level = "profile"}, sparse latent profiles can be combined into
#' grouped profile states through \code{profile_group}. The grouped analysis
#' aggregates the latent-profile prior distribution, hard classifications, and
#' profile misclassification matrix before fitting the distal outcome model.
#'
#' The BCH correction replaces the full mixture likelihood with BCH weights
#' derived from the inverse misclassification matrix. Both the naive analysis
#' that ignores classification error and the requested corrected analyses are
#' returned.
#'
#' This function is currently limited to models with binary attributes.
#'
#' @return
#' A list containing the fitted distal outcome analyses. \code{results}
#' contains the requested distal outcome analysis for the selected attribute or
#' profile states. Each fit includes coefficient estimates, Wald summaries,
#' fitted class-specific outcome summaries, threshold estimates for ordinal
#' outcomes, any observed-covariate terms included in the model, and the
#' correction objects used for estimation.
#'
#' @references
#' Bakk, Z., & Kuha, J. (2021). Relating latent class membership to external
#' variables: An overview. \emph{British Journal of Mathematical and Statistical
#' Psychology, 74}(2), 340-362.
#'
#' Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent structure
#' models with categorical variables: One-step versus three-step estimators.
#' \emph{Political Analysis, 12}(1), 3-27.
#'
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
#' three-step approaches. \emph{Political Analysis, 18}(4), 450-469.
#'
#' @importFrom GDINA CM extract indlogPost personparm
#'
#' @examples
#' \dontrun{
#'
#' #########################################################
#' # example 1: binary distal outcome at the attribute level
#' #########################################################
#' library(GDINA)
#' set.seed(123)
#' N <- 3000
#'
#' Q <- sim10GDINA$simQ
#' gs <- matrix(runif(nrow(Q)*2,0,0.2), nrow(Q), 2)
#'
#' sim <- simGDINA(N, Q, gs.parm = gs, model = "GDINA")
#' fit <- GDINA(sim$dat, sim$Q, verbose = 0)
#' alpha <- extract(sim, "attribute")
#'
#' y_binary <- rbinom(N, 1, plogis(-0.8 + alpha[, 1] - 0.9 * alpha[, 2] +
#'   0.7 * alpha[, 3]))
#' distal_att <- ThreeStepDistal(fit, y_binary, level = "attribute")
#' print(distal_att)
#'
#' ###################################################################
#' # example 2: continuous distal outcome with profile-level analysis
#' ###################################################################
#'
#' y_continuous <- 0.3 + 0.9 * alpha[, 1] - 1.1 * alpha[, 2] +
#'   0.8 * alpha[, 3] + rnorm(N, sd = 1.1)
#' distal_cont <- ThreeStepDistal(fit,y_continuous,
#'   level = "profile",method = "BCH",reference = "first")
#' distal_cont
#'
#' ###################################################################
#' # example 3: Binary distal outcome with observed covariates
#' ###################################################################
#' set.seed(456)
#' covariates <- data.frame(
#'   gender = factor(sample(c("female", "male"), N, replace = TRUE)),
#'   ses = ordered(sample(c("low", "medium", "high"), N, replace = TRUE),
#'     levels = c("low", "medium", "high")),age = round(rnorm(N,40,10)))
#' y_cov <- rbinom(N,1,
#'   plogis(-1.4 - 0.8 * alpha[, 1] - 0.5 * alpha[, 2] -
#'     0.5 * (covariates$gender == "male") -
#'     0.35 * as.numeric(covariates$ses) + 0.2 * covariates$age))
#' distal_cov <- ThreeStepDistal(fit,y_cov, formula = ~ gender + ses + age,
#'   data = covariates,level = "attribute", attribute = 1:2,
#'   method = "BCH")
#' distal_cov
#'
#'
#' ###################################################################
#' # example 4: combined profile-level analysis with observed covariates
#' ###################################################################
#'
#' distal_cov_comb <- ThreeStepDistal(fit,y_cov, formula = ~ gender + ses + age,
#'   data = covariates,level = "profile",
#'   profile_group = c(0,1,1,2,2,3,3,4),
#'   method = "BCH")
#'  distal_cov_comb
#'
#'  ##### same but with labels and named vector for profile grouping
#' distal_cov_comb2 <- ThreeStepDistal(fit,y_cov, formula = ~ gender + ses + age,
#'   data = covariates,level = "profile",
#'   profile_group = c("000" = "G0","100" = "G1","010" = "G1","001" = "G2","110" = "G2",
#'                     "101" = "G3","011" = "G3","111" = "G4"),
#'   method = "BCH")
#'  distal_cov_comb2
#'
#' }
#' @export
ThreeStepDistal <- function(object,
                            outcome,
                            formula = NULL,
                            data = NULL,
                            level = c("attribute", "profile"),
                            attribute = NULL,
                            classification = "MAP",
                            method = c("BCH","ML"),
                            outcome_type = c("auto", "binary", "ordinal", "categorical", "nominal", "continuous"),
                            reference = c("last", "first"),
                            profile_group = NULL,
                            conf.level = 0.95,
                            maxit = 1000,
                            na.action = getOption("na.action")) {


  level <- match.arg(level)
  method <- unique(match.arg(method, c("ML", "BCH"), several.ok = TRUE))
  reference <- match.arg(reference)
  inputs <- distal_prepare_inputs(
    object = object,
    outcome = outcome,
    outcome_type = outcome_type,
    classification = classification,
    formula = formula,
    data = data,
    na.action = na.action
  )

  if (level == "attribute") {
    natt <- ncol(inputs$classification$hard_pattern)
    if (is.null(attribute)) {
      attribute <- seq_len(natt)
    }
    attribute <- unique(as.integer(attribute))
    if (length(attribute) == 0L || any(is.na(attribute)) || any(attribute < 1L) || any(attribute > natt)) {
      stop("attribute must contain valid attribute indices.", call. = FALSE)
    }

    state_info <- distal_attribute_state_info(
      pattern = inputs$pattern,
      hard_pattern = inputs$classification$hard_pattern,
      posterior = inputs$posterior,
      attribute = attribute
    )
    design_info <- distal_state_design(
      level = "attribute",
      n_states = nrow(state_info$state_pattern),
      labels = state_info$state_labels,
      attribute = attribute,
      state_pattern = state_info$state_pattern
    )
    fits <- distal_fit_bundle(
      outcome_info = inputs$outcome,
      state_design = design_info$design,
      covariate_design = inputs$covariates$design,
      observed_state = state_info$hard_class,
      misclassification = state_info$misclassification,
      prior = state_info$prior,
      method = method,
      conf.level = conf.level,
      maxit = maxit
    )

    results <- c(list(
      attribute = attribute,
      observed = state_info$hard_pattern,
      baseline = design_info$baseline,
      prior = stats::setNames(state_info$prior, design_info$state_labels),
      misclassification = state_info$misclassification,
      predictor_labels = distal_predictor_term_names(design_info$design, inputs$covariates$design),
      state_labels = design_info$state_labels,
      covariate_terms = colnames(inputs$covariates$design)
    ), fits)
  } else {
    state_info <- distal_profile_state_info(
      pattern = inputs$pattern,
      hard_class = inputs$classification$hard_class,
      prior = colMeans(inputs$posterior),
      misclassification = CM(
        object,
        classification = inputs$classification$classification,
        matrixtype = "profile"
      )$profile_classification,
      profile_group = profile_group
    )
    design_info <- distal_state_design(
      level = "profile",
      n_states = nrow(state_info$misclassification),
      labels = state_info$state_labels,
      reference = reference
    )
    fits <- distal_fit_bundle(
      outcome_info = inputs$outcome,
      state_design = design_info$design,
      covariate_design = inputs$covariates$design,
      observed_state = state_info$hard_class,
      misclassification = state_info$misclassification,
      prior = state_info$prior,
      method = method,
      conf.level = conf.level,
      maxit = maxit
    )

    results <- c(list(
      observed = state_info$hard_class,
      profile = inputs$classification$hard_pattern,
      baseline = design_info$baseline,
      prior = stats::setNames(state_info$prior, design_info$state_labels),
      misclassification = state_info$misclassification,
      predictor_labels = distal_predictor_term_names(design_info$design, inputs$covariates$design),
      state_labels = design_info$state_labels,
      profile_group = state_info$profile_group,
      profile_labels = state_info$profile_labels,
      reference = reference,
      covariate_terms = colnames(inputs$covariates$design)
    ), fits)
  }

  ret <- list(
    call = match.call(),
    formula = formula,
    level = level,
    classification = classification,
    method = method,
    outcome_type = inputs$outcome$type,
    outcome_levels = inputs$outcome$levels,
    na.action = inputs$covariates$na.action,
    covariate_design = inputs$covariates$design,
    results = results
  )
  class(ret) <- c("ThreeStepDistal", "threeStepDistal")
  return(ret)
}


#' @export
print.ThreeStepDistal <- function(x, ...) {
  fit_name <- distal_print_method(x)
  fit <- x$results[[fit_name]]

  cat("\nThree-step distal outcome analysis\n")
  cat("Level:", x$level, "\n")
  cat("Outcome type:", x$outcome_type, "\n")
  cat("Correction methods:", paste(x$method, collapse = ", "), "\n")
  cat("Displayed fit:", fit_name, "\n")

  if (x$level == "attribute") {
    cat("Attributes in regression:", paste0("Attribute ", x$results$attribute, collapse = ", "), "\n")
  } else {
    cat("Reference profile:", x$results$baseline, "\n")
  }

  if (!is.null(x$formula)) {
    cat("Observed covariates:", paste(attr(stats::terms(x$formula), "term.labels"), collapse = ", "), "\n")
  }

  if (!is.null(fit$table)) {
    cat("\nCoefficient table\n")
    print(three_step_format_print_table(fit$table), row.names = FALSE)
  }

  if (!is.null(fit$threshold_table)) {
    cat("\nThreshold table\n")
    print(three_step_format_print_table(fit$threshold_table), row.names = FALSE)
  }

  invisible(x)
}
