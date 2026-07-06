#' @title Three-step ML correction for covariate regression
#'
#' @description
#' Fit covariate regression models for estimated attribute classifications or
#' latent profile classifications from a fitted
#' \code{\link{GDINA}} object using the maximum-likelihood three-step
#' correction (Vermunt, 2010). The correction uses the estimated misclassification matrix
#' \eqn{P(W=s\mid X=t)} obtained from \code{\link{CM}} function.
#'
#' @param object An estimated \code{GDINA} object returned from
#'   \code{\link{GDINA}}.
#' @param formula A one-sided or two-sided formula specifying the covariates.
#'   Only the right-hand side is used.
#' @param data A data frame containing the covariates in \code{formula}. It
#'   must contain one row per respondent in \code{object}.
#' @param level A character string specifying whether the corrected regression
#'   is carried out at the attribute level (\code{"attribute"}) or the latent
#'   profile level (\code{"profile"}).
#' @param attribute Optional integer vector giving the attributes to analyze when
#'   \code{level = "attribute"}. If \code{NULL}, all attributes are analyzed.
#' @param classification A character string specifying the classification rule.
#'   Supported values are \code{"MAP"}, \code{"MLE"}, and \code{"EAP"}.
#'   Alternatively, a matrix of user-supplied classifications can be provided,
#'   with one row per respondent and one column per attribute.
#' @param reference A character string specifying the reference class for
#'   profile-level regression. Supported values are \code{"last"}
#'   (default) and \code{"first"}. Ignored when
#'   \code{level = "attribute"}.
#' @param profile_group Optional vector or factor used to combine latent
#'   profiles when \code{level = "profile"}. Supply one value per profile in
#'   the order returned by \code{extract(object, "attributepattern")};
#'   profiles sharing the same value are collapsed into one grouped state.
#'   A named vector is also accepted and will be matched to profile labels
#'   such as \code{"000"} or \code{"101"}. If \code{NULL}, each profile is
#'   analyzed separately.
#' @param conf.level Confidence level used for Wald intervals.
#' @param na.action See \code{\link{na.action}} for details on handling missing data.
#'
#' @details
#' With notations from Vermunt (2010), let \eqn{W} denote the estimated class,
#' \eqn{X} the unobserved true class,
#' and \eqn{Z} the covariates. The corrected likelihood implemented here is
#' based on
#' \deqn{\sum_i \log\left\{\sum_t P(W_i=s_i\mid X_i=t)P(X_i=t\mid Z_i)\right\},}
#' which corresponds to the three-step ML correction. The posterior
#' class weights used after fitting are based on
#' \deqn{P(X_i=t\mid W_i=s_i,Z_i) \propto P(W_i=s_i\mid X_i=t)P(X_i=t\mid Z_i),}
#' At the attribute level, a separate binary logistic
#' regression is fitted for each selected attribute. At the profile level, a
#' multinomial logistic regression is fitted using the latent profiles as the
#' outcome categories. Sparse latent profiles can be combined into grouped
#' profile states through \code{profile_group}. The grouped analysis
#' aggregates the latent-profile prior distribution, hard classifications, and
#' profile misclassification matrix before fitting the multinomial model.
#'
#' This function is currently limited to models with binary attributes.
#'
#' @return
#' A list containing the corrected and naive regression fits. For
#' \code{level = "attribute"}, \code{results} is a named list with one entry
#' per fitted attribute. For \code{level = "profile"}, \code{results} contains
#' the multinomial regression fit. Each fit includes coefficient estimates,
#' Wald summaries, fitted probabilities, posterior class weights, and the
#' misclassification matrix used for correction.
#'
#' @references
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved three-step approaches. \emph{Political analysis, 18}(4), 450-469.
#'
#' @importFrom GDINA CM extract personparm
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' N <- 1000
#' Z <- data.frame(
#'   z_cont = rnorm(N),
#'   z_bin = rbinom(N, 1, 0.45),
#'   z_cat = factor(sample(c("A", "B", "C"), N, replace = TRUE))
#' )
#'
#' p1 <- plogis(-0.3 + 0.9 * Z$z_cont + 0.8 * Z$z_bin - 0.5 * (Z$z_cat == "B") +
#'   0.4 * (Z$z_cat == "C"))
#' p2 <- plogis(0.2 - 0.8 * Z$z_cont + 0.5 * Z$z_bin + 0.6 * (Z$z_cat == "B") -
#'   0.3 * (Z$z_cat == "C"))
#' alpha <- cbind(rbinom(N, 1, p1), rbinom(N, 1, p2))
#'
#' Q <- matrix(c(
#'   1, 0,
#'   1, 0,
#'   1, 0,
#'   0, 1,
#'   0, 1,
#'   0, 1,
#'   1, 1,
#'   1, 1
#' ), byrow = TRUE, ncol = 2)
#' gs <- data.frame(guess = rep(0.15, nrow(Q)), slip = rep(0.15, nrow(Q)))
#' sim <- GDINA::simGDINA(N, Q, gs.parm = gs, model = "DINA", attribute = alpha)
#' fit <- GDINA::GDINA(sim$dat, sim$Q, model = "DINA", verbose = 0)
#'
#' # Attribute-level logistic regression
#' ts_att <- ThreeStepCov(fit, ~ z_cont + z_bin + z_cat, data = Z, level = "attribute")
#' print(ts_att)
#'
#' # Profile-level multinomial logistic regression
#' ts_profile <- ThreeStepCov(
#'   fit,
#'   ~ z_cont + z_bin + z_cat,
#'   data = Z,
#'   level = "profile",
#'   reference = "first"
#' )
#' print(ts_profile)
#'
#' # Combined profile-level regression
#' ts_profile_grouped <- ThreeStepCov(
#'   fit,
#'   ~ z_cont + z_bin + z_cat,
#'   data = Z,
#'   level = "profile",
#'   profile_group = c(0, 1, 1, 2),
#'   reference = "first"
#' )
#' print(ts_profile_grouped)
#' }
#' @export
ThreeStepCov <- function(object,
                         formula,
                         data,
                         level = c("attribute", "profile"),
                         attribute = NULL,
                         classification = "MAP",
                         reference = c("last", "first"),
                         profile_group = NULL,
                         conf.level = 0.95,
                         na.action = getOption("na.action")) {
  if (!inherits(object, "GDINA")) {
    stop("object must be a GDINA object.", call. = FALSE)
  }

  level <- match.arg(level)
  reference <- match.arg(reference)

  design_info <- three_step_design_matrix(
    formula = formula,
    data = data,
    nobs = extract(object, "nobs"),
    na.action = na.action
  )
  design <- design_info$design
  class_info <- three_step_classification(object, classification)
  if (length(design_info$keep) != nrow(class_info$hard_pattern)) {
    class_info$hard_pattern <- class_info$hard_pattern[design_info$keep, , drop = FALSE]
    class_info$hard_class <- class_info$hard_class[design_info$keep]
  }

  if (any(class_info$hard_pattern > 1)) {
    stop("ThreeStepCov is currently only available for binary attributes.", call. = FALSE)
  }

  if (level == "attribute") {
    misclassification <- CM(
      object,
      classification = class_info$classification,
      matrixtype = "attribute"
    )$att_classification
    natt <- ncol(class_info$hard_pattern)
    if (is.null(attribute)) {
      attribute <- seq_len(natt)
    }
    attribute <- as.integer(attribute)
    if (any(is.na(attribute)) || any(attribute < 1L) || any(attribute > natt)) {
      stop("attribute must contain valid attribute indices.", call. = FALSE)
    }

    results <- lapply(attribute, function(k) {
      observed <- class_info$hard_pattern[, k]
      naive_fit <- three_step_binary_fit_core(
        design = design,
        observed = observed,
        misclassification = diag(2),
        conf.level = conf.level
      )
      corrected_fit <- three_step_binary_fit_core(
        design = design,
        observed = observed,
        misclassification = misclassification[[k]],
        start = naive_fit$coefficients,
        conf.level = conf.level
      )

      list(
        attribute = k,
        observed = observed,
        misclassification = misclassification[[k]],
        naive = naive_fit,
        corrected = corrected_fit
      )
    })
    names(results) <- paste0("Attribute ", attribute)
  } else {
    state_info <- distal_profile_state_info(
      pattern = class_info$pattern,
      hard_class = class_info$hard_class,
      prior = {
        posterior <- exp(extract(object, "logposterior.i"))
        posterior <- posterior / rowSums(posterior)
        colMeans(posterior[design_info$keep, , drop = FALSE])
      },
      misclassification = CM(
        object,
        classification = class_info$classification,
        matrixtype = "profile"
      )$profile_classification,
      profile_group = profile_group
    )
    naive_fit <- three_step_multinomial_fit_core(
      design = design,
      observed = state_info$hard_class,
      misclassification = diag(nrow(state_info$misclassification)),
      class_labels = state_info$state_labels,
      reference = reference,
      conf.level = conf.level
    )
    corrected_fit <- three_step_multinomial_fit_core(
      design = design,
      observed = state_info$hard_class,
      misclassification = state_info$misclassification,
      class_labels = state_info$state_labels,
      start = c(naive_fit$coefficients),
      reference = reference,
      conf.level = conf.level
    )
    results <- list(
      observed = state_info$hard_class,
      profile = class_info$hard_pattern,
      misclassification = state_info$misclassification,
      prior = stats::setNames(state_info$prior, state_info$state_labels),
      state_labels = state_info$state_labels,
      profile_group = state_info$profile_group,
      profile_labels = state_info$profile_labels,
      reference = reference,
      naive = naive_fit,
      corrected = corrected_fit
    )
  }

  ret <- list(
    call = match.call(),
    formula = formula,
    level = level,
    classification = classification,
    reference = if (level == "profile") reference else NULL,
    profile_group = if (level == "profile") profile_group else NULL,
    na.action = design_info$na.action,
    design = design,
    results = results
  )
  class(ret) <- c("ThreeStepCov")
  return(ret)
}

#' @export
print.ThreeStepCov <- function(x, ...) {
  cat("\nThree-step covariate regression\n")
  cat("Level:", x$level, "\n")

  if (x$level == "attribute") {
    cat("Attributes:", paste(names(x$results), collapse = ", "), "\n")
  } else {
    cat("Reference profile:", x$results$corrected$reference, "\n")
    if (!is.null(x$results$profile_group)) {
      cat("Grouped profiles:", paste(unique(unname(x$results$profile_group)), collapse = ", "), "\n")
    }
  }

  cat("\nCorrected coefficients\n")
  print(three_step_format_print_table(three_step_cov_print_table(x)), row.names = FALSE)
  invisible(x)
}
