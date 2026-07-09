# Marginal Effects (Impacts) ----

#' @title Compute Marginal Effects for Spatial Binary  Models
#' 
#' @aliases impacts
#' @import spatialreg
#' @export impacts
#' 
#' @description
#' Computes marginal effects (direct, indirect, total, from-region, and cumulative) from a spatial binary response model.
#' Supports inference via the Delta method or Monte Carlo simulation.
#' 
#' @param obj An object of class \code{bingmm}, \code{binlgmm} or \code{binris}. 
#' @param object An object of class \code{impacts.bingmm} for \code{summary} methods.
#' @param x An object of class \code{impacts.bingmm} for \code{print} methods. 
#' @param data Optional dataset used for computing partial effects.
#' @param variable Variable name (string) for which marginal effects are computed.
#' @param from.unit Integer specifying the spatial unit (row) for computing impacts.
#' @param dydx Character: type of derivative approximation to use (\code{"exact"} or \code{"numeric"}).
#' @param change For numeric variables, a vector indicating the lower and upper value to compute the partial change. 
#' @param vce A string indicating what kind of variance-covariance matrix of the estimate should be computed when using \code{effect.bingmm}. For the one-step GMM estimator, the options are \code{"robust"} and \code{"ml"}. For the two-step GMM estimator, the options are \code{"robust"}, \code{"efficient"} and \code{"ml"}. The option \code{"vce = ml"} is an exploratory method that evaluates the VC of the RIS estimator using the GMM estimates.
#' @param type String indicating which method is used to compute the standard errors of the marginal effects. If \code{"mc"}, then the Monte Carlo approximation is used. If \code{"delta"}, then the Delta Method is used.
#' @param result String indicating the output format: \code{"summary"} presents the average direct, indirect or total effects; \code{"from.region"} display the total effect of the variable in \code{"variable"} form spatial unit indicated in \code{"from.unit"}; if \code{"cumulative"} the cumulative effect up to power order indicated in \code{"Q"} for the variable in \code{"variable"} form spatial unit indicated in \code{"from.unit"} is returned.
#' @param vcov Optional user-supplied variance-covariance matrix.
#' @param het Logical. If \code{TRUE} (the default), then the heteroskedasticity is taken into account when computing the marginal effects. 
#' @param empirical Logical. Argument passed to \code{mvrnorm} (default \code{FALSE}): if \code{TRUE}, the coefficients and their covariance matrix specify the empirical not population mean and covariance matrix.
#' @param R Numerical. Indicates the number of draws used in the Monte Carlo approximation if \code{type = "mc"}.
#' @param approximation Logical. If \code{TRUE} then \eqn{(I - \lambda W)^{-1}} is approximated as \eqn{I + \lambda W + \lambda^2 W^2 + \lambda^3 W^3 + ... +\lambda^q W^q}. The default is \code{FALSE}.
#' @param pw numeric. The power used for the approximation \eqn{I + \lambda W + \lambda^2 W^2 + \lambda^3 W^3 + ... +\lambda^q W^q}. The default is 5.
#' @param draws Optional matrix of parameter draws for Monte Carlo simulation.
#' @param Q Number of columns for cumulative effects (if applicable).
#' @param tol Numerical. Argument passed to \code{mvrnorm}: tolerance (relative to largest variance) for numerical lack of positive-definiteness in the coefficient covariance matrix.
#' @param verbose Logical. Display progress.
#' @param digits the number of digits.
#' @param ... further arguments. Ignored.
#' 
#' @details 
#' 
#' Let the model be:
#' 
#' \deqn{
#' y^*= X\beta + WX\gamma + \lambda W y^* + \epsilon = Z\delta + \lambda Wy^{*} + \epsilon
#' }
#' 
#' where  \eqn{y = 1} if \eqn{y^*>0} and 0 otherwise; \eqn{\epsilon \sim N(0, 1)} if \code{link = "probit"} or \eqn{\epsilon \sim L(0, \pi^2/3)} if \code{link = "logit"}. 
#' The RIS estimator assumes that \eqn{\epsilon \sim N(0, 1)}. 
#' 
#' The marginal effects respect to variable \eqn{x_r} can be computed as
#' 
#' \deqn{
#' diag(f(a))D^{-1}_{\lambda}A^{-1}_{\lambda}\left(I_n\beta_r + W\gamma_r\right) = C_r(\theta)
#' }
#' 
#' where \eqn{f()} is the pdf, which depends on the assumption of the error terms; \eqn{diag} is the operator that creates a \eqn{n \times n} diagonal matrix; \eqn{A_{\lambda}= (I -\lambda W)}; and \eqn{D_{\lambda}} is a diagonal matrix whose elements represent the square root of the diagonal elements of the variance-covariance matrix of \eqn{u = A_{\lambda}^{-1}\epsilon}. 
#' 
#' We implement these three summary measures: (1) The average total effects, \eqn{ATE_r  = n^{-1}i_n'C_{r}i_n}, (2) The average direct effects, \eqn{ADE_r  = n^{-1}tr(C_{r})}, and (3) the average indirect effects, \eqn{ATE_r - ADE_r}. 
#' 
#' The standard errors of the average total, direct and indirect effects can be estimated using either Monte Carlo (MC) approximation, which takes into account the sampling distribution of \eqn{\theta}, or Delta Method. 
#' 
#' @examples
#' \donttest{
#' # Data set
#' data(oldcol, package = "spdep")
#' 
#' # Create dependent (dummy) variable
#' COL.OLD$CRIMED <- as.numeric(COL.OLD$CRIME > 35)
#' 
#' # Two-step (Probit) GMM estimator
#' ts <- sbinaryGMM(CRIMED ~ INC * HOVAL + factor(CP) | HOVAL,
#'                 link = "probit", 
#'                 listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                 data = COL.OLD, 
#'                 type = "twostep")
#'                 
#' # Summary marginal effects using Delta Method
#' summary(impacts(ts, type = "delta"))
#' 
#' # Summary marginal effects using MC with 100 draws
#' summary(impacts(ts, type = "mc", R = 100))
#' 
#' # Marginal effects using efficient VC matrix
#' summary(impacts(ts, type = "delta", vce = "efficient"))
#' 
#' # Marginal effects using efficient VC matrix and ignoring the heteroskedasticity
#' summary(impacts(ts, type = "delta", vce = "efficient", het = FALSE))
#' 
#' # Obtain total impact when INC changes in region 20
#' imp <- impacts(ts, type = "delta", result = "from.region", from.unit = 20, variable = "INC")
#' head(imp$out, 6)
#' 
#' # Cumulative impact when INC changes in region 20 up to power = 3
#' cum <- impacts(ts, type = "delta", result = "cumulative", from.unit = 20, variable = "INC", Q = 3)
#' head(cum$out$`Q:3`, 6)
#' 
#' # Marginal effects using  RIS estimator
#' ris_sar <- sbinaryRis(CRIMED ~ INC + HOVAL, data = COL.OLD,
#'                       R = 50, 
#'                       listw = spdep::nb2listw(COL.nb, style = "W"))
#' summary(impacts(ris_sar, method = "delta"))
#' summary(impacts(ris_sar, method = "mc", R = 100))
#'}
#' @return An object of class \code{impacts.bingmm}. 
#' @seealso \code{\link[spldv]{sbinaryGMM}}, \code{\link[spldv]{sbinaryLGMM}}.
#' @author Mauricio Sarrias and Gianfranco Piras. 
#' @method impacts bingmm
#' @importFrom numDeriv jacobian
#' @importFrom MASS mvrnorm
#' @keywords marginal effects
#' @export
impacts.bingmm <- function(obj,
                           data          = NULL, 
                           variable      = NULL, 
                           from.unit     = 1, 
                           dydx          = c("exact", "numeric"),
                           change        = NULL, 
                           vce           = c("robust", "efficient", "ml"),
                           type          = c("delta", "mc"),
                           result        = c("summary", "from.region", "cumulative"),
                           vcov          = NULL,
                           het           = TRUE,
                           empirical     = FALSE,
                           approximation = FALSE,
                           R             = 50,
                           pw            = 5, 
                           draws         = NULL,
                           Q             = 3,
                           tol           = 1e-06,
                           verbose       = FALSE,
                           ...) {
  
  # Match arguments
  type   <- match.arg(type)
  vce    <- match.arg(vce)
  result <- match.arg(result)
  dydx   <- match.arg(dydx)
  
  # Get parameters
  mu       <- coef(obj)
  n_params <- length(mu)
  
  # Obtain VCOV matrix
  V <- if (is.null(vcov)) {
    if (!inherits(obj, "bingmm")){
      vcov(obj)
    } else {
      vcov(obj, vce = vce)
    }
  } else {
    if (!is.matrix(vcov) || any(dim(vcov) != c(n_params, n_params))) {
      stop("Supplied vcov has incorrect dimensions")
    }
    vcov
  }
  
  if (is.null(data)) data <- obj$data
  
  # Stop early for formula transformations for which spatial marginal effects
  # are currently not supported reliably.
  check_impacts_formula_terms(obj)
  
  # Prepare derivative arguments
  der_args <- list(
    variable      = variable, 
    data          = data, 
    obj           = obj, 
    from.unit     = from.unit, 
    result        = result, 
    dydx          = dydx, 
    change        = change, 
    het           = het, 
    approximation = approximation,
    pw            = pw, 
    Q             = Q, 
    verbose       = verbose
  )
  
  compute_dydx <- function(theta) {
    do.call(dydx.bingmm, c(list(theta = theta), der_args))
  }
  
  compute_se <- function(jacobian) {
    sqrt(diag(jacobian %*% V %*% t(jacobian)))
  }
  
  add_stats <- function(est, se, rownms) {
    z <- est / se
    p <- 2 * pnorm(-abs(z))
    data.frame(dydx = est, `Std. error` = se, `z value` = z,
               `Pr(> z)` = p, row.names = rownms)
  }
  
  # ===== Delta Method ==== #
  if (type == "delta") {
    me  <- compute_dydx(mu)
    if (result == "summary"){
      jac   <- numDeriv::jacobian(compute_dydx, mu)
      out   <- add_stats(me, compute_se(jac), names(me))
      
    } else if (result == "from.region"){
      compute_jac <- if (is.list(me)) {
        function(theta) unlist(compute_dydx(theta))
      } else {
        compute_dydx
      }
      jac <- numDeriv::jacobian(compute_jac, mu)
      se  <- compute_se(jac)
      
      if (is.list(me)) {
        lengths_me <- sapply(me, length)
        stopifnot(length(unique(lengths_me)) == 1)
        n <- lengths_me[1]
        se_split <- split(se, rep(names(me), each = n))
        
        out <- mapply(function(m, s) add_stats(m, s, rownames(data)),
                      me, se_split, SIMPLIFY = FALSE)
      } else {
        out <- add_stats(me, se, rownames(data))
      }
    } else if  (result == "cumulative"){
      compute_jac <- if (is.list(me)) {
        function(theta) unlist(compute_dydx(theta))
      } else {
        function(theta) as.vector(compute_dydx(theta))
      }
      
      jac    <- numDeriv::jacobian(compute_jac, mu)
      se_vec <- compute_se(jac)
      
      if (is.list(me)) {
        se_blocks <- Map(function(m, start) {
          dims <- dim(m)
          matrix(se_vec[start:(start + prod(dims) - 1)], nrow = dims[1])
        }, me, cumsum(c(1, head(sapply(me, function(x) prod(dim(x))), -1))))
        
        out <- Map(function(m_mat, s_mat) {
          Q <- ncol(m_mat)
          setNames(lapply(seq_len(Q), function(q)
            add_stats(m_mat[, q], s_mat[, q], rownames(data))),
            colnames(m_mat))
        }, me, se_blocks)
      } else{
        se_mat <- matrix(se_vec, nrow = nrow(me), ncol = ncol(me))
        out   <- setNames(lapply(seq_len(ncol(me)), function(q)
          add_stats(me[, q], se_mat[, q], rownames(data))),
          colnames(me))
      }
    }
  } else if (type == "mc") {
    W            <- obj$listw
    sym          <- all(W == t(W))
    omega        <- eigen(W, only.values = TRUE, symmetric = sym)
    eig_range    <- if (is.complex(omega$values)) range(Re(omega$values)) else 1 / range(omega$values)
    lambda_range <- 1 / eig_range
    
    draws        <- if (is.null(draws)) MASS::mvrnorm(n = R, mu = mu, Sigma = V, tol = tol, empirical = empirical) else draws
    lambda_pos   <- length(mu)
    valid        <- draws[, lambda_pos] > lambda_range[1] & draws[, lambda_pos] < lambda_range[2]
    valid_draws  <- draws[valid, , drop = FALSE]
    
    if (nrow(valid_draws) < R) {
      warning("Some draws discarded due to invalid lambda. R reduced to ", nrow(valid_draws), ".")
    }
    
    sres_list <- lapply(seq_len(nrow(valid_draws)), function(i) compute_dydx(valid_draws[i, ]))
    
    # === Post-processing Monte Carlo Results ===
    if (result == "summary") {
      # Convert list to matrix: each row = draw, each column = marginal effect
      mat_draws <- do.call(rbind, sres_list)
      me        <- colMeans(mat_draws)
      se        <- apply(mat_draws, 2, sd)
      z         <- me / se
      p         <- 2 * pnorm(-abs(z))
      
      out <- data.frame(
        `dydx`        = me,
        `Std. error`  = se,
        `z value`     = z,
        `Pr(> z)`     = p,
        row.names     = colnames(mat_draws)
      )
    } else if (result == "from.region") {
      # Can return list or numeric vector
      if (is.list(sres_list[[1]])) {
        names_me <- names(sres_list[[1]])
        n        <- length(sres_list[[1]][[1]])
        out      <- list()
        
        for (j in names_me) {
          # Build matrix of R x n for each factor level
          mat_j <- do.call(rbind, lapply(sres_list, function(x) x[[j]]))
          me_j  <- colMeans(mat_j)
          se_j  <- apply(mat_j, 2, sd)
          z     <- me_j / se_j
          p     <- 2 * pnorm(-abs(z))
          out[[j]] <- data.frame(
            `dydx`        = me_j,
            `Std. error`  = se_j,
            `z value`     = z,
            `Pr(> z)`     = p,
            row.names     = rownames(data)
          )
        }
      } else {
        # sres_list returns numeric vector (length n)
        mat_draws <- do.call(rbind, sres_list)  # R x n
        me        <- colMeans(mat_draws)
        se        <- apply(mat_draws, 2, sd)
        z         <- me / se
        p         <- 2 * pnorm(-abs(z))
        
        out <- data.frame(
          `dydx`        = me,
          `Std. error`  = se,
          `z value`     = z,
          `Pr(> z)`     = p,
          row.names     = rownames(data)
        )
      }
    } else if (result == "cumulative") {
      # Can return list or matrix
      if (is.list(sres_list[[1]])) {
        names_me <- names(sres_list[[1]])
        out      <- list()
        
        for (j in names_me) {
          Q        <- ncol(sres_list[[1]][[j]])
          n        <- nrow(sres_list[[1]][[j]])
          mat_list <- vector("list", Q)
          
          for (q in seq_len(Q)) {
            # R x n matrix for this factor and level
            mat_q <- do.call(rbind, lapply(sres_list, function(x) x[[j]][, q]))
            me_q  <- colMeans(mat_q)
            se_q  <- apply(mat_q, 2, sd)
            z     <- me_q / se_q
            p     <- 2 * pnorm(-abs(z))
            mat_list[[q]] <- data.frame(
              `dydx`        = me_q,
              `Std. error`  = se_q,
              `z value`     = z,
              `Pr(> z)`     = p,
              row.names     = rownames(data)
            )
          }
          names(mat_list) <- colnames(sres_list[[1]][[j]])
          out[[j]]        <- mat_list
        }
      } else {
        Q        <- ncol(sres_list[[1]])
        mat_list <- vector("list", Q)
        
        for (q in seq_len(Q)) {
          mat_q <- do.call(rbind, lapply(sres_list, function(x) x[, q]))  # R x n
          me_q  <- colMeans(mat_q)
          se_q  <- apply(mat_q, 2, sd)
          z     <- me_q / se_q
          p     <- 2 * pnorm(-abs(z))
          mat_list[[q]] <- data.frame(
            `dydx`        = me_q,
            `Std. error`  = se_q,
            `z value`     = z,
            `Pr(> z)`     = p,
            row.names     = rownames(data)
          )
        }
        names(mat_list) <- colnames(sres_list[[1]])
        out <- mat_list
      }
    }
  }
  ## Saving results
  result <- structure(
    list(
      out    = out,
      result = result
    ), 
    class = "impacts.bingmm"
  )

  return(result)
}



#' @rdname impacts.bingmm
#' @method impacts binlgmm
#' @importFrom numDeriv jacobian
#' @importFrom MASS mvrnorm
#' @export
impacts.binlgmm <- impacts.bingmm

#' @rdname impacts.bingmm
#' @method impacts binris
#' @importFrom numDeriv jacobian
#' @importFrom MASS mvrnorm
#' @export
impacts.binris <- impacts.bingmm


# Check formula terms that are not currently supported by impacts().
# These transformations can make exact and numerical spatial marginal effects
# ambiguous or unreliable because the derivative/counterfactual must be
# propagated through the full spatial multiplier.
check_impacts_formula_terms <- function(object) {
  ftxt <- paste(deparse(object$formula), collapse = " ")

  unsupported_patterns <- c(
    poly  = "(^|[^[:alnum:]_.:])(?:stats::)?poly\\s*\\(",
    scale = "(^|[^[:alnum:]_.:])(?:base::)?scale\\s*\\(",
    bs    = "(^|[^[:alnum:]_.:])(?:splines::)?bs\\s*\\(",
    ns    = "(^|[^[:alnum:]_.:])(?:splines::)?ns\\s*\\("
  )

  hits <- names(unsupported_patterns)[
    vapply(unsupported_patterns, grepl, logical(1), x = ftxt, perl = TRUE)
  ]

  if (length(hits) > 0L) {
    stop(
      paste0(
        "impacts() does not currently support marginal effects for formulas ",
        "containing the following transformed terms: ",
        paste0(hits, "()", collapse = ", "), ".\n",
        "These terms can produce unreliable marginal effects in spatial binary ",
        "models because the derivative or discrete change must be propagated ",
        "through the spatial multiplier.\n",
        "Please rewrite the formula with explicit variables before fitting the ",
        "model. For quadratic terms, use x + I(x^2) instead of ",
        "poly(x, 2, raw = TRUE). For scaled variables, create the scaled column ",
        "in the data first, for example data$x_s <- as.numeric(scale(data$x)), ",
        "fit the model using x_s, and request impacts(..., variable = 'x_s')."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' @rdname impacts.bingmm
#' @method print impacts.bingmm
#' @export
print.impacts.bingmm <- function(x, ... ){
  cat("The results are:\n")
  print(x$out)
}

#' @rdname impacts.bingmm
#' @method summary impacts.bingmm
#' @export
summary.impacts.bingmm <- function(object, ...){
  result         <- object$result
  if (result == "summary"){
    CoefTable      <- object$out
    summary        <- list(CoefTable = CoefTable)
    class(summary) <- "summary.impacts.bingmm"
    summary
  } else {
    stop("S3 method summary is only for scalar summary statistics")
  }
}

#' @rdname impacts.bingmm
#' @method print summary.impacts.bingmm
#' @export
print.summary.impacts.bingmm <- function(x, digits = max(3, getOption("digits") - 3), ...){
  k <- nrow(x$CoefTable) / 3
  cat("------------------------------------------------------", fill = TRUE)
  cat("(a) Total effects :\n")
  cat("------------------------------------------------------",fill = TRUE)
  
  printCoefmat(x$CoefTable[1:k, , drop = FALSE], digits = digits, signif.legend = FALSE)
  
  cat("\n------------------------------------------------------", fill = TRUE)
  cat("(b) Direct effects :\n")
  cat("------------------------------------------------------",fill = TRUE)
  
  printCoefmat(x$CoefTable[(k + 1):(k*2), , drop = FALSE], digits = digits, signif.legend = FALSE)
  
  cat("\n------------------------------------------------------", fill = TRUE)
  cat("(c) Indirect effects :\n")
  cat("------------------------------------------------------",fill = TRUE)
  
  printCoefmat(x$CoefTable[(2*k + 1):(k*3), , drop = FALSE], digits = digits)
}


dydx.bingmm <- function(theta, 
                        variable      = NULL, 
                        data          = NULL, 
                        obj           = NULL, 
                        from.unit     = 1, 
                        result        = c("summary", "from.region", "cumulative"), 
                        dydx          = c("exact", "numeric"), 
                        change        = NULL,
                        het           = TRUE, 
                        approximation = FALSE, 
                        pw            = 5,
                        Q             = 3,
                        verbose       = TRUE,
                        ...){
  # This function generates the partial effects for binary spatial models
  # The first argument must be the vector of coefficients to construct also the Jacobian
  
  if (is.null(obj)) stop("object NULL")
  check_impacts_formula_terms(obj)
  
  # Match argument and check conditions
  result <- match.arg(result)
  dydx   <- match.arg(dydx)
  
  # Preliminary checks
  if (result %in% c("from.region", "cumulative") && is.null(variable)) {
    stop(sprintf("'%s' estimates require a specified variable.", result))
  }
  if (!is.null(change)) {
    if (result == "summary") {
      stop("'change' parameter is only allowed when dydx = 'numeric' and result != 'summary'")
    }
    # Additional validation for change values could go here
    # }
  }
  
  n <- nrow(obj$X)
  if (result %in% c("from.region", "cumulative")) {
    if (!is.numeric(from.unit) || length(from.unit) != 1L ||
        is.na(from.unit) || !is.finite(from.unit) ||
        from.unit < 1L || from.unit > n || from.unit != floor(from.unit)) {
      stop("'from.unit' must be a single valid spatial unit index between 1 and n.")
    }
    from.unit <- as.integer(from.unit)
  }
  
  # Data
  if (is.null(data)) data <- obj$data
  
  # Compute Sinv only if needed (not for cumulative)
  Sinv <- NULL
  if (result != "cumulative") {
    lambda <- theta["lambda"]
    A      <- Matrix::Diagonal(n) - lambda * obj$listw
    Sinv   <- if (approximation) app_W(obj$listw, lambda, pw) else Matrix::solve(A)
  }
  
  
  # Get classes of variables: numeric, logical and factor
  tvars               <- get.class.var(obj)
  var_type_map        <- rep(names(tvars), times = sapply(tvars, length))
  names(var_type_map) <- unlist(tvars)
  var_names           <- names(var_type_map)
  
  if (result %in% c("from.region", "cumulative")) {
    if (!is.character(variable) || length(variable) != 1L ||
        is.na(variable) || !(variable %in% var_names)) {
      stop("'variable' must be the name of a variable used in the fitted model.")
    }
  }
  
  # --- 1) Summary effect (numeric) ----
  if (result == "summary" && dydx == "numeric"){
    TE_accum <- numeric()
    DE_accum <- numeric()
    names_accum <- character()
    
    # Get variable types
    #var_names <- names(var_type_map)
    for (v in var_names){
      type.var <- var_type_map[[v]]
      
      if (type.var == "fnames"){
        levs <- nlevels(factor(data[[v]])) - 1L
        for (j in seq_len(levs)){
          key <- paste0(v, "_lev", j)
          TE_accum[key] <- 0
          DE_accum[key] <- 0
        }
      } else {
        TE_accum[v] <- 0
        DE_accum[v] <- 0
      }
    }
    
    #Loop over observations
    for (i in seq_len(n)){
      if (verbose && i %% 100 == 0) cat("Processing unit", i, "of", n, "\n")
      
      for (v in var_names){
        type.var <- var_type_map[[v]]
        
        if (type.var == "nnames"){
          der <- dydx.num.sp(theta = theta, object = obj, from.unit = i, 
                             variable = v, Sinv  = Sinv, het = het, data = data)
          TE_accum[v] <- TE_accum[v] + sum(der)
          DE_accum[v] <- DE_accum[v] + der[i]
          
        } else if (type.var == "lnames"){
          der <- dydx.logical.sp(theta = theta, object = obj, from.unit = i, 
                                 variable = v, Sinv  = Sinv, het = het, data = data)
          TE_accum[v] <- TE_accum[v] + sum(der)
          DE_accum[v] <- DE_accum[v] + der[i]
          
        } else if (type.var == "fnames"){
          der_list <- dydx.factor.sp(theta = theta, object = obj, from.unit = i, 
                                     variable = v, Sinv  = Sinv, het = het, data = data)
          for(j in seq_along(der_list)){
            key <- paste0(v, "_lev", j)
            TE_accum[key] <- TE_accum[key] + sum(der_list[[j]])
            DE_accum[key] <- DE_accum[key] + der_list[[j]][i]
          }
        }
      }
    }
    TE <- TE_accum / n
    DE <- DE_accum / n
    IE <- TE - DE
    names(TE) <- paste0("TE:", names(TE))
    names(DE) <- paste0("DE:", names(DE))
    names(IE) <- paste0("IE:", names(IE))
    out <- c(TE, DE, IE)
    return(out)
  }
  # --- 2) Summary effect (exact) ----
  if (result == "summary" && dydx == "exact"){
    out <- dydx.exact.sp(theta         = theta, 
                         object        = obj,
                         data          = data,
                         Sinv          = Sinv,
                         het           = het, 
                         approximation = approximation, 
                         pw = pw)
  }  
  # --- 3) From region effect ----
  if (result == "from.region"){
    # Return a vector of marginal impact from region = `from.unit`
    type.var <- var_type_map[[variable]]
    switch(type.var, 
           nnames = {
             if (dydx == "exact" && is.null(change)) {
               out <- dydx.exact.num.from.sp(theta = theta, object = obj,
                                             from.unit = from.unit,
                                             variable = variable, Sinv = Sinv,
                                             het = het, data = data)
             } else {
               out <- dydx.num.sp(theta = theta, object = obj, from.unit = from.unit,
                                  variable = variable, Sinv = Sinv, het = het,
                                  data = data, change = change)
             }
           }, 
           lnames = {
             out <- dydx.logical.sp(theta = theta, object = obj, from.unit = from.unit, 
                                    variable = variable, Sinv  = Sinv, het = het, data = data)
           }, 
           fnames = {
             #levs <- length(levels(factor(data[[variable]]))) -  1L
             #out  <- vector("list", levs)
             out <- dydx.factor.sp(theta = theta, object = obj, from.unit = from.unit, 
                                   variable = variable, Sinv = Sinv, het = het, data = data)
           },
           stop("Unhandled variable type: ", type.var)
    )
    #names(out) <- rownames(data)
  }
  # --- 4) From region cumulative effect ----
  if (result == "cumulative"){
    type.var <- var_type_map[[variable]]
    levs     <- if (type.var == "fnames") length(levels(factor(data[[variable]]))) - 1L else NULL
    
    # Build W once as matrix
    #W_mat <- as.matrix(obj$listw)
    
    # Precompute power of W
    W_powers <- vector("list", Q)
    W <- obj$listw
    if (Q > 0) {
      W_powers[[1]] <- W
      if (Q > 1){
        for (q in 2:Q) W_powers[[q]] <- W_powers[[q-1]] %*% W
      }
    }
    
    # Result structure
    out <- switch(type.var,
                  "fnames" = {
                    mat_list <- lapply(seq_len(levs), function(j) {
                      matrix(NA_real_, n, Q + 1L, dimnames = list(rownames(data), paste0("Q:", 0:Q)))
                    })
                    factor_levels <- levels(factor(data[[variable]]))[-1L]
                    names(mat_list) <- factor_levels
                    mat_list
                  },
                  matrix(NA_real_, n, Q + 1L, dimnames = list(rownames(data), paste0("Q:", 0:Q)))
    )
    
    lambda <- theta["lambda"]
    
    # Accumulate powers without computing inverse
    Sinv_q_prev <- Matrix::Diagonal(n)
    
    for (q in 0:Q) {
      if (q == 0) {
        Sinv_q <- Matrix::Diagonal(n)
      } else {
        # Use precomputed sparse powers with scalar multiplication
        Sinv_q <- Sinv_q_prev + (lambda^q) * W_powers[[q]]
        Sinv_q_prev <- Sinv_q
      }
      
      if (type.var == "nnames") {
        if (dydx == "exact" && is.null(change)) {
          out[, q + 1] <- dydx.exact.num.from.sp(theta = theta, object = obj,
                                                 from.unit = from.unit,
                                                 variable = variable, Sinv = Sinv_q,
                                                 het = het, data = data)
        } else {
          out[, q + 1] <- dydx.num.sp(theta = theta, object = obj, from.unit = from.unit,
                                      variable = variable, Sinv = Sinv_q, het = het,
                                      data = data, change = change)
        }
      } else if (type.var == "lnames") {
        out[, q + 1] <- dydx.logical.sp(theta = theta, object = obj, from.unit = from.unit, 
                                        variable = variable, Sinv = Sinv_q, het = het, data = data)
      } else if (type.var == "fnames") {
        der_list <- dydx.factor.sp(theta = theta, object = obj, from.unit = from.unit, 
                                   variable = variable, Sinv = Sinv_q, het = het, data = data)
        for (j in seq_len(levs)){
          out[[j]][, q + 1] <- der_list[[j]]
        }
      }
    }
  }
  return(out)
}


# Helper functions  ----
get.class.var <- function(obj, ...){
  # Identify classes of variables in model (by mf)
  classes <- attributes(terms(obj$mf))[["dataClasses"]][-1L]
  
  # Drop if "(weights)" variables, just in case 
  classes <- classes[!names(classes) %in% "(weights)"]
  
  # handle character variables as factors
  classes[classes == "character"] <- "factor"
  
  # Clean names in terms (From margins)
  clean_terms <- function(terms) {
    # the use of paste("`", x, "`") is a hack to deal with variables that have spaces in their names
    unlist(lapply(terms, function(x) all.vars(formula(paste0("~", ifelse(grepl("[[:alnum:]_.] [[:alnum:]_.]", x), paste0("`", x, "`"), x))))))
  }
  names(classes) <- clean_terms(names(classes))
  
  # List identifying numeric, logical and factors
  vars <- list(
    nnames = unique(names(classes)[!classes %in% c("factor", "ordered", "logical")]),
    lnames = unique(names(classes)[classes == "logical"]),
    fnames = unique(names(classes)[classes %in% c("factor", "ordered")])
  )
  
  # check whether the list is completely NULL
  if (is.null(unlist(vars))) {
    stop("No variables found in model.")
  }
  return(vars)
}

resolve_variable <- function(other_var, data, model_data) {
  # Case 1: Regular variable in model_data
  if (other_var %in% colnames(model_data)) {
    return(model_data[, other_var])
  }
  
  # Case 2: Directly in raw data
  if (other_var %in% names(data)) {
    return(data[[other_var]])
  }
  
  # Case 3: factor(var)level
  factor_match <- regexec("^factor\\(([^\\)]+)\\)(.+)$", other_var)
  factor_tokens <- regmatches(other_var, factor_match)[[1]]
  if (length(factor_tokens) == 3) {
    base <- factor_tokens[2]
    level <- factor_tokens[3]
    if (base %in% names(data)) {
      return(as.numeric(as.factor(data[[base]]) == level))
    }
  }
  
  # Case 4: I(expression)
  I_match <- regexec("^I\\((.+)\\)$", other_var)
  I_tokens <- regmatches(other_var, I_match)[[1]]
  if (length(I_tokens) == 2) {
    expr_text <- I_tokens[2]
    expr <- try(parse(text = expr_text)[[1]], silent = TRUE)
    if (!inherits(expr, "try-error")) {
      env <- as.environment(data)
      return(eval(expr, envir = env))
    }
  }
  
  # Case 5: poly(var, degree)index
  poly_match <- regexec("^poly\\(([^,]+),\\s*(\\d+)\\)(\\d+)$", other_var)
  poly_tokens <- regmatches(other_var, poly_match)[[1]]
  if (length(poly_tokens) == 4) {
    base <- poly_tokens[2]
    degree <- as.numeric(poly_tokens[3])
    index <- as.numeric(poly_tokens[4])
    if (base %in% names(data)) {
      poly_basis <- poly(data[[base]], degree = degree, raw = TRUE)
      return(poly_basis[, index])
    }
  }
  
  warning(sprintf("Could not resolve variable '%s'", other_var))
  return(rep(NA_real_, nrow(data)))
}

reconstruct.factor.interactions <- function(f.lev, terms, data) {
  Xvals <- list()
  
  for (tn in terms) {
    parts <- strsplit(tn, ":")[[1]]
    
    if (length(parts) == 1) {
      # Simple term (factor or lagged factor)
      if (startsWith(tn, "lag_")) {
        # E.g., lag_factor(z)1
        fvar <- gsub("^lag_factor\\((.+)\\)(.+)$", "\\1", tn)
        lev  <- gsub("^lag_factor\\(.+\\)(.+)$", "\\1", tn)
        tmp  <- 1
        Xvals[[tn]] <- tmp
      } else {
        # E.g., factor(z)1
        fvar <- gsub("^factor\\((.+)\\)(.+)$", "\\1", tn)
        lev  <- gsub("^factor\\(.+\\)(.+)$", "\\1", tn)
        tmp  <- 1
        Xvals[[tn]] <- tmp
      }
      
    } else if (length(parts) == 2) {
      # Interaction terms (x:factor or lag_factor:x, etc.)
      first <- parts[1]
      second <- parts[2]
      
      if (startsWith(first, "lag_")) {
        # Case: lag_factor(z)1:x or similar
        fvar <- gsub("^lag_factor\\((.+)\\)(.+)$", "\\1", first)
        lev  <- gsub("^lag_factor\\(.+\\)(.+)$", "\\1", first)
        tmp1 <- 1
        tmp2 <- data[[second]]
        Xvals[[tn]] <- tmp1 * tmp2
        
      } else if (startsWith(second, "lag_")) {
        # Case: x:lag_factor(z)1
        fvar <- gsub("^lag_factor\\((.+)\\)(.+)$", "\\1", second)
        lev  <- gsub("^lag_factor\\(.+\\)(.+)$", "\\1", second)
        tmp1 <- data[[first]]
        tmp2 <- 1
        Xvals[[tn]] <- tmp1 * tmp2
        
      } else if (startsWith(first, "factor(")) {
        # Case: factor(z)1:x
        fvar <- gsub("^factor\\((.+)\\)(.+)$", "\\1", first)
        lev  <- gsub("^factor\\(.+\\)(.+)$", "\\1", first)
        tmp1 <- 1
        tmp2 <- data[[second]]
        Xvals[[tn]] <- tmp1 * tmp2
        
      } else if (startsWith(second, "factor(")) {
        # Case: x:factor(z)1
        fvar <- gsub("^factor\\((.+)\\)(.+)$", "\\1", second)
        lev  <- gsub("^factor\\(.+\\)(.+)$", "\\1", second)
        tmp1 <- data[[first]]
        tmp2 <- 1
        Xvals[[tn]] <- tmp1 * tmp2
        
      } else {
        stop(paste("Unsupported interaction format in:", tn))
      }
      
    } else {
      stop(paste("Too many interaction components in term:", tn))
    }
  }
  
  # Return a data.frame with the terms
  as.data.frame(Xvals)
}

# Functions for numerical derivatives ----
dydx.num.sp <- function(theta = NULL,
                        data, 
                        object, 
                        variable,
                        from.unit = 1, 
                        change = NULL,
                        Sinv = NULL, 
                        het  = TRUE,
                        approximation = FALSE,
                        pw  = 5,
                        eps = 1e-7, ...){
  # Check if change is numeric with two values
  if (is.numeric(change)){
    stopifnot(length(change) == 2)
    lwr <- change[1]
    upr <- change[2]
    change <- "discrete"
  } else {
    change <- "numerical"
  }
  
  # Data
  if (missing(data) || is.null(data)) data <- object$data
  
  
  # Set value of h based on eps (based on margins)
  # set.h <- function(x){
  #   x + (max(abs(x), 1, na.rm = TRUE) * sqrt(eps)) - x
  # }
  
  # Generate S matrix if `NULL` to speed up
  if (is.null(Sinv)){
    if(is.null(theta)) lambda <- coef(object)["lambda"] else lambda <- theta["lambda"]
    n         <- nrow(object$listw)
    A         <- Diagonal(n) - lambda * object$listw
    Sinv      <- if (approximation) app_W(object$listw, lambda, pw)  else Matrix::solve(A)
  }
  
  # Set value of h based on eps
  # FIXME: Improve how h is computed (see margins)
  #set.h <- function(x){
  #  (abs(mean(x)) + eps) * eps
  #}
  
  x_orig <- data[[variable]]
  
  if (change == "numerical"){
    h <- (abs(mean(x_orig)) + eps) * eps
    d0 <- x_orig
    d1 <- x_orig
    d0[from.unit] <- d0[from.unit] - h
    d1[from.unit] <- d1[from.unit] + h
    data[[variable]] <- d0
    y.hat.0 <- predict(object = object, newdata = data, Sinv = Sinv, het = het, theta = theta)
    data[[variable]] <- d1
    y.hat.1 <- predict(object = object, newdata = data, Sinv = Sinv, het = het, theta = theta)
    partial.ef <- (y.hat.1 - y.hat.0) / (2 * h)
  } else if (change == "discrete"){
    d0 <- x_orig
    d1 <- x_orig
    d0[from.unit] <- lwr
    d1[from.unit] <- upr
    data[[variable]] <- d0
    y.hat.0 <- predict(object = object, newdata = data, Sinv = Sinv, het = het, theta = theta)
    data[[variable]] <- d1
    y.hat.1 <- predict(object = object, newdata = data, Sinv = Sinv, het = het, theta = theta)
    partial.ef <- y.hat.1 - y.hat.0
  } else if (change == "exact"){
    stop("not yet implemented")
  }
  return(partial.ef)
}

dydx.factor.sp <- function(theta = NULL, 
                           object, 
                           variable,
                           data, 
                           from.unit = 1, 
                           Sinv = NULL, 
                           het  = TRUE,
                           approximation = FALSE,
                           pw  = 5,
                           ...){
  # Extract data and basic setup
  if (missing(data) || is.null(data)) data <- object$data
  x         <- data[[variable]]
  levs_full <- levels(as.factor(x))
  base      <- levs_full[1L]
  levs      <- levs_full[-1L]
  n         <- nrow(data)
  L         <- length(levs)
  
  # # Extract levels and set base
  # levs <- levels(as.factor(data[[variable]]))
  # base <- levs[1L]
  # levs <- levs[-1L] # Exclude base
  # 
  
  # Generate S matrix if `NULL` to speed up
  if (is.null(Sinv)){
    if(is.null(theta)) lambda <- coef(object)["lambda"] else lambda <- theta["lambda"]
    A         <- Matrix::Diagonal(n) - lambda * object$listw
    Sinv      <- if (approximation) app_W(object$listw, lambda, pw)  else Matrix::solve(A)
  }
  
  # Compute baseline prediction
  data[[variable]] <- base
  y.hat.0 <- predict(object = object, newdata = data, Sinv = Sinv, het = het, theta = theta)
  
  #Preallocate output
  partial.ef <- vector("list", L)
  
  #Loop over factor levels (excluding base)
  for (l in seq_len(L)){
    d1 <- data
    d1[from.unit, variable] <- levs[l]
    y.hat.1 <- predict(object = object, newdata = d1, Sinv = Sinv, het = het, theta = theta)
    partial.ef[[l]] <- y.hat.1 - y.hat.0
  }
  names(partial.ef) <- paste0("factor(", variable, ")", levs)
  return(partial.ef)
}

dydx.ordered.sp <- dydx.factor.sp

dydx.logical.sp <- function(theta = NULL,
                            object, 
                            variable,
                            data,
                            from.unit = 1, 
                            Sinv = NULL, 
                            het  = TRUE,
                            approximation = FALSE,
                            pw  = 5,
                            ...){
  # Data
  if (missing(data) || is.null(data)) data <- object$data
  
  
  # Set value of h based on eps (based on margins)
  # set.h <- function(x){
  #   x + (max(abs(x), 1, na.rm = TRUE) * sqrt(eps)) - x
  # }
  
  # Generate S matrix if `NULL` to speed up
  if (is.null(Sinv)){
    if(is.null(theta)) lambda <- coef(object)["lambda"] else lambda <- theta["lambda"]
    n         <- nrow(object$listw)
    A         <- Diagonal(n) - lambda * object$listw
    Sinv      <- if (approximation) app_W(object$listw, lambda, pw)  else Matrix::solve(A)
  }
  
  # Data with changes
  d0 <- data
  d0[[variable]] <- FALSE
  y.hat.0 <- predict(object = object, newdata = d0, Sinv = Sinv, het = het, theta = theta)
  
  d1 <- d0
  d1[from.unit, variable] <- TRUE
  y.hat.1 <- predict(object = object, newdata = d1, Sinv = Sinv, het = het, theta = theta)
  partial.ef <- y.hat.1 - y.hat.0
  
  return(partial.ef)
}

dydx.exact.num.from.sp <- function(theta = NULL,
                                    object,
                                    variable,
                                    from.unit = 1,
                                    data = NULL,
                                    Sinv = NULL,
                                    het = TRUE,
                                    approximation = FALSE,
                                    pw = 5,
                                    ...){
  # Exact continuous from-region impact. This returns the column of the
  # analytical probability-impact matrix associated with `from.unit`.
  # Discrete changes, logical variables, and factors continue to use the
  # predictive contrast path.
  if (missing(data) || is.null(data)) data <- object$data
  if (is.null(theta)) theta <- coef(object)

  lambda <- theta["lambda"]
  betas  <- theta[names(theta) != "lambda"]
  X      <- object$X
  n      <- nrow(X)
  W      <- object$listw

  dfun <- switch(object$link, "probit" = dnorm, "logit" = dlogis)

  if (is.null(Sinv)) {
    A    <- Matrix::Diagonal(n) - lambda * W
    Sinv <- if (approximation) app_W(W, lambda, pw) else Matrix::solve(A)
  }

  Di <- NULL
  if (het) {
    Di <- Matrix::Diagonal(x = 1 / sqrt(Matrix::rowSums(Sinv^2)))
  }

  XB <- as.vector(X %*% betas)
  a  <- if (het) as.vector(Di %*% (Sinv %*% XB)) else as.vector(Sinv %*% XB)
  fa <- dfun(a)
  dM <- Matrix::Diagonal(x = fa) %*% (if (het) Di %*% Sinv else Sinv)

  chain.coef <- dydx.chain.num.sp(var = variable, theta = theta,
                                  model_data = X, data = data, W = W)
  out <- as.numeric(dM %*% chain.coef[, from.unit, drop = FALSE])
  names(out) <- rownames(data)
  return(out)
}

# Functions for exact derivatives ----
dydx.exact.sp <- function(theta, 
                          object,
                          data          = NULL,
                          Sinv          = NULL,
                          het           = TRUE, 
                          approximation = FALSE, 
                          pw = 5, 
                          ...){
  # Extract coefficients 
  lambda    <- theta["lambda"]
  betas     <- theta[names(theta) != "lambda"]
  
  # Retrieve data
  if (is.null(data)) data <- object$data
  X        <- object$X
  N        <- nrow(X)
  W        <- object$listw
  
  # Set density and distribution functions
  dfun <- switch(object$link, "probit" = dnorm, "logit"  = dlogis)
  pfun <- switch(object$link, "probit" = pnorm, "logit"  = plogis)
  
  
  # Compute Sinv if not supplied
  if (is.null(Sinv)){
    A         <- Matrix::Diagonal(N) - lambda * W
    Sinv      <- if (approximation) app_W(W, lambda, pw)  else Matrix::solve(A)
  }
  
  # Variable classes
  tvars               <- get.class.var(object)
  var_type_map        <- rep(names(tvars), times = sapply(tvars, length))
  names(var_type_map) <- unlist(tvars)
  var_names           <- names(var_type_map)
  has_numeric         <- any(var_type_map == "nnames")
  
  # Compute heteroskedasticity adjustment matrix Di (if needed)
  Di <- NULL
  if (het) {
    Di <- Matrix::Diagonal(x = 1 / sqrt(Matrix::rowSums(Sinv^2)))
  }
  
  # Precompute dM for numeric variables only
  if (has_numeric) {
    XB <- as.vector(X %*% betas)
    a  <- if (het) as.vector(Di %*% (Sinv %*% XB)) else as.vector(Sinv %*% XB)
    fa <- dfun(a)
    dM <- Matrix::Diagonal(x = fa) %*% (if (het) Di %*% Sinv else Sinv)
  } else {
    dM <- NULL
  }
  
  # Initialize lists for effects
  TE <- DE <- list()
  
  for (v in var_names){
    type.var <- var_type_map[[v]]
    
    if (type.var == "nnames"){
      chain.coef <- dydx.chain.num.sp(var = v, theta = theta, model_data = X, data = data, W = W)
      ME <- dM %*% chain.coef
      TE[[v]] <- sum(ME) / N
      DE[[v]] <- sum(diag(ME)) / N
      
    } else if (type.var == "fnames"){
      res.factor <- dydx.chain.fac.sp(v = v, X = X, betas = betas,
                                      Sinv = Sinv, Di = Di, data = data,
                                      N = N, pfun = pfun, het = het, W = W)
      TE <- c(TE, res.factor$TE)
      DE <- c(DE, res.factor$DE)
    }
  }
  # Compute IE
  IE <- mapply(`-`, TE, DE, SIMPLIFY = FALSE)
  
  # Rename elements
  names(TE) <- paste0("TE:", names(TE))
  names(DE) <- paste0("DE:", names(DE))
  names(IE) <- paste0("IE:", names(IE))
  
  # Output
  out <- c(unlist(TE), unlist(DE), unlist(IE))
  return(out)
}

dydx.chain.num.sp <- function(var, typevar, theta, model_data, data, W, ...) {
  n           <- nrow(model_data)
  theta_names <- names(theta)
  da_dx       <- Matrix::Matrix(0, n, n, sparse = TRUE)
  
  x_vals <- data[[var]]
  
  # 1. Linear term: var ----
  if (var %in% theta_names) {
    da_dx <- da_dx + Matrix::Diagonal(n, theta[[var]])
  }
  
  # 2. Linear term: lag_var (WX) ----
  lag_var <- paste0("lag_", var)
  if (lag_var %in% theta_names) {
    da_dx <- da_dx + theta[[lag_var]] * W
  }
  
  # --- 3. Quadratic terms like I(x^2), I(x^3), etc. ----
  idx_I <- grep(paste0("^I\\(", var, "\\^(\\d+)\\)$"), theta_names)
  if (length(idx_I) > 0) {
    powers <- as.numeric(sub(paste0("^I\\(", var, "\\^(\\d+)\\)$"), "\\1", theta_names[idx_I]))
    coefs  <- theta[idx_I]
    deriv_vals <- rowSums(mapply(function(c, p) c * p * x_vals^(p - 1), coefs, powers))
    da_dx <- da_dx + Matrix::Diagonal(n, deriv_vals)
  }
  
  # --- 4. Quadratic terms in WX ----
  idx_lagI <- grep(paste0("^lag_I\\(", var, "(\\^\\d+)?\\)$"), theta_names)
  if (length(idx_lagI) > 0) {
    raw_powers <- sub(paste0("^lag_I\\(", var, "(\\^\\d+)?\\)$"), "\\1", theta_names[idx_lagI])
    powers <- ifelse(raw_powers == "", 1, as.numeric(gsub("\\^", "", raw_powers)))
    coefs  <- theta[idx_lagI]
    for (i in seq_along(powers)) {
      da_dx <- da_dx + coefs[i] * powers[i] * W %*% Matrix::Diagonal(x = x_vals^(powers[i] - 1))
    }
  }
  
  # --- 5. poly(x, degree) terms ----
  idx_poly <- grep(paste0("poly\\(", var, ",\\s*\\d+\\)(\\d+)"), theta_names)
  if (length(idx_poly) > 0) {
    warning("poly() terms detected. Marginal effects are approximate unless original basis is available.")
    for (name in theta_names[idx_poly]) {
      if (name %in% names(model_data)) {
        da_dx <- da_dx + Matrix::Diagonal(n, theta[[name]] * model_data[[name]])
      }
    }
  }
  
  # --- 6. Interactions involving var or lag_var ---
  interaction_patterns <- c(
    paste0("(^", var, ":)|(:", var, "$)|(:", var, ":)"),
    paste0("(^lag_", var, ":)|(:lag_", var, "$)|(:lag_", var, ":)"),
    paste0("lag_.*:", var), paste0(var, ":lag_.*")
  )
  idx_int <- unique(unlist(lapply(interaction_patterns, function(p) grep(p, theta_names))))
  
  if (length(idx_int) > 0) {
    for (iname in theta_names[idx_int]) {
      vars_in_term <- unlist(strsplit(iname, ":"))
      is_lagged <- grepl("^lag_", vars_in_term)
      clean_vars <- gsub("^lag_", "", vars_in_term)
      
      if (any(clean_vars == var & is_lagged)) {
        # lag_var:x or lag_var:other
        other_var <- vars_in_term[clean_vars != var]
        other_var_clean <- gsub("^lag_", "", other_var)
        x_other <- resolve_variable(other_var_clean, data, model_data)
        
        if (!any(is.na(x_other))) {
          da_dx <- da_dx + theta[[iname]] * W %*% Matrix::Diagonal(x = x_other)
        } else {
          warning(sprintf("Could not resolve '%s' in interaction term '%s'", other_var, iname))
        }
        
      } else if (any(clean_vars == var & !is_lagged)) {
        # var:x or var:lag_other
        other_var <- vars_in_term[clean_vars != var]
        other_var_clean <- gsub("^lag_", "", other_var)
        x_other <- resolve_variable(other_var_clean, data, model_data)
        
        if (!any(is.na(x_other))) {
          if (grepl("^lag_", other_var)) {
            da_dx <- da_dx + theta[[iname]] * W %*% Matrix::Diagonal(x = x_other)
          } else {
            da_dx <- da_dx + Matrix::Diagonal(n, theta[[iname]] * x_other)
          }
        } else {
          warning(sprintf("Could not resolve '%s' in interaction term '%s'", other_var, iname))
        }
      }
    }
  }
  return(da_dx)
}

dydx.chain.fac.sp <- function(v, X, betas, Sinv, Di, data, N, pfun, het, W, ...) {
  TE <- DE <- list()
  
  # Extract original values and levels
  x         <- data[[v]]
  levs_full <- levels(as.factor(x))
  base      <- levs_full[1L]
  levs      <- levs_full[-1L]
  L         <- length(levs)
  
  # Locate all coefficients related to the factor variable
  f.base        <- paste0("factor(", v, ")")
  f.idxs        <- grep(f.base, names(betas), fixed = TRUE)
  if (length(f.idxs) == 0) {
    stop(paste0("No coefficients found for ", f.base, ". Check your model matrix."))
  }
  
  # Build base prediction (all dummies  = 0)
  X.0           <- X
  X.0[, f.idxs] <- 0
  b.f.0         <- Matrix::Diagonal(x = drop(X.0[, f.idxs, drop = FALSE] %*% betas[f.idxs]))
  b.rest        <- betas[-f.idxs]
  X.b.rest      <- drop(X[, -f.idxs, drop = FALSE] %*% b.rest)
  R.X.b.rest    <- matrix(rep(X.b.rest, N), nrow = N)
  
  T.0 <- Sinv %*% (R.X.b.rest + b.f.0)
  P.0 <- if (het) pfun(as.matrix(Di %*% T.0)) else pfun(as.matrix(T.0))
  
  for (l in 1:L) {
    lev_name      <- levs[l]
    # Key base names
    f.lev            <- paste0("factor(", v, ")", lev_name)
    lag.f.lev        <- paste0("lag_factor(", v, ")", lev_name)
    
    # Get all related variable names from the beta vector
    all.relevant <- names(betas)[grepl(f.lev, names(betas), fixed = TRUE) |
                                   grepl(lag.f.lev, names(betas), fixed = TRUE)]
    
    if (length(all.relevant) == 0) {
      warning(paste0("No coefficients found for level ", lev_name, ". Skipping."))
      next
    }
    
    # Reconstruct X values for interactions, etc.
    X.l <- X.0
    vals.int <- tryCatch({
      reconstruct.factor.interactions(f.lev, all.relevant, data)
    }, error = function(e) {
      stop(paste0("Error reconstructing interactions for ", f.lev, ": ", e$message))
    })
    
    X.l[, all.relevant] <- as.matrix(vals.int)
    
    # Separate terms into lagged vs non-lagged
    lag.terms    <- grep("^lag_", all.relevant, value = TRUE)
    nonlag.terms <- setdiff(all.relevant, lag.terms)
    
    # Compute contributions
    b.f.l <- Matrix::Matrix(0, N, N, sparse = TRUE)
    
    if (length(nonlag.terms) > 0) {
      b.f.l <- b.f.l + Matrix::Diagonal(x = drop(X.l[, nonlag.terms, drop = FALSE] %*% betas[nonlag.terms, drop = TRUE]))
    }
    
    if (length(lag.terms) > 0) {
      for (kk in lag.terms) {
        xval <- drop(X.l[, kk])
        if (!all(xval == 0)) {  # skip useless computation
          b.f.l <- b.f.l + betas[kk] * W %*% Matrix::Diagonal(x = xval)
        }
      }
    }
    
    # Compute new predictions 
    T.l           <- Sinv %*% (R.X.b.rest + b.f.l)
    P.l           <- if (het) pfun(as.matrix(Di %*% T.l)) else pfun(as.matrix(T.l))
    ME.l          <- P.l - P.0
    key           <- paste0(v, "_lev", levs[[l]])
    TE[[key]]     <- sum(ME.l) / N
    DE[[key]]     <- sum(Matrix::diag(ME.l)) / N
  }
  
  list(TE = TE, DE = DE)
}



