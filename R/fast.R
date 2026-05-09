#
# New file: fast / Rcpp-backed variants of core flexmix routines.
# Original functions in flexmix.R are unchanged so behaviour can be
# compared and verified.
#

## Fast replacement for the internal log_row_sums()
log_row_sums_fast <- function(m) {
  cpp_log_row_sums(m)
}

## Fast replacement for groupPosteriors()
groupPosteriors_fast <- function(x, group) {
  if (length(group) > 0) {
    g <- as.integer(group)
    cpp_group_posteriors(x, g, length(unique(g)))
  } else x
}

## Fast classification step (weighted/hard/CEM only; SEM/random fall back)
.FLXgetOK_fast <- function(p, control, weights) {
  cls <- control@classify
  if (is.null(weights)) {
    if (cls == "weighted") return(p)
    if (cls %in% c("CEM", "hard"))
      return(cpp_hard_assign(p, numeric(0)))
  } else {
    if (cls == "weighted") return(p * weights)
    if (cls %in% c("CEM", "hard"))
      return(cpp_hard_assign(p, as.numeric(weights)))
  }
  ## SEM / random: defer to original
  flexmix:::.FLXgetOK(p, control, weights)
}

###**********************************************************
## Fast FLXfit: same algorithm as setMethod("FLXfit","list",...) in flexmix.R
## but using cpp_estep / cpp_hard_assign / cpp_group_posteriors.
##
## Returns a standard "flexmix" S4 object so all downstream methods work.
###**********************************************************
FLXfit_fast <- function(model, concomitant, control,
                        postunscaled = NULL, groups, weights) {
  k <- ncol(postunscaled)
  N <- nrow(postunscaled)
  control <- allweighted(model, control, weights)
  if (control@verbose > 0)
    cat("Classification:", control@classify, " (fast)\n")
  
  group       <- groups$group
  groupfirst  <- groups$groupfirst
  gf_idx0     <- which(groupfirst) - 1L  # 0-based for C++
  
  if (length(group) > 0)
    postunscaled <- groupPosteriors_fast(postunscaled, group)
  
  logpostunscaled <- log(postunscaled)
  postscaled <- exp(logpostunscaled - log_row_sums_fast(logpostunscaled))
  
  llh <- -Inf
  converged <- FALSE
  components <- rep(list(rep(list(new("FLXcomponent")), k)), length(model))
  
  for (iter in seq_len(control@iter.max)) {
    
    ## --- classification / weighting step ---
    postscaled <- .FLXgetOK_fast(postscaled, control, weights)
    
    ## --- concomitant fit -> prior ---
    prior <- if (is.null(weights))
      ungroupPriors(concomitant@fit(concomitant@x,
                                    postscaled[groupfirst, , drop = FALSE]),
                    group, groupfirst)
    else
      ungroupPriors(concomitant@fit(concomitant@x,
                                    (postscaled / weights)[groupfirst & weights > 0, , drop = FALSE],
                                    weights[groupfirst & weights > 0]),
                    group, groupfirst)
    
    ## --- min.prior pruning (unchanged from original) ---
    nok <- if (nrow(prior) == 1) which(prior < control@minprior) else {
      if (is.null(weights))
        which(colMeans(prior[groupfirst, , drop = FALSE]) < control@minprior)
      else
        which(colSums(prior[groupfirst, ] * weights[groupfirst]) /
                sum(weights[groupfirst]) < control@minprior)
    }
    if (length(nok)) {
      if (control@verbose > 0)
        cat("*** Removing", length(nok), "component(s) ***\n")
      prior <- prior[, -nok, drop = FALSE]
      prior <- prior / rowSums(prior)
      postscaled <- postscaled[, -nok, drop = FALSE]
      zeros <- rowSums(postscaled) == 0
      postscaled[zeros, ] <- if (nrow(prior) > 1) prior[zeros, ]
      else prior[rep(1, sum(zeros)), ]
      postscaled <- postscaled / rowSums(postscaled)
      if (!is.null(weights)) postscaled <- postscaled * weights
      k <- ncol(prior)
      if (k == 0) stop("all components removed")
      model      <- lapply(model, FLXremoveComponent, nok)
      components <- lapply(components, "[", -nok)
    }
    
    ## --- M-step (unchanged) ---
    components <- lapply(seq_along(model),
                         function(i) FLXmstep(model[[i]], postscaled, components[[i]]))
    
    postunscaled <- matrix(0, nrow = N, ncol = k)
    for (n in seq_along(model))
      postunscaled <- postunscaled +
      FLXdeterminePostunscaled(model[[n]], components[[n]])
    
    if (length(group) > 0)
      postunscaled <- groupPosteriors_fast(postunscaled, group)
    
    ## --- E-step (fused in C++) ---
    es <- cpp_estep(postunscaled,
                    prior,
                    gf_idx0,
                    if (is.null(weights)) numeric(0) else as.numeric(weights))
    postscaled  <- es$postscaled
    postunscaled <- es$postunscaled
    llh.old <- llh
    llh     <- es$llh
    
    if (is.na(llh) || is.infinite(llh))
      stop(paste(formatC(iter, width = 4), "Log-likelihood:", llh))
    
    if (abs(llh - llh.old) / (abs(llh) + 0.1) < control@tolerance) {
      if (control@verbose > 0) { printIter(iter, llh); cat("converged\n") }
      converged <- TRUE
      break
    }
    if (control@verbose && (iter %% control@verbose == 0))
      printIter(iter, llh)
  }
  
  ## --- assemble return object (mirrors original) ---
  components <- lapply(seq_len(k),
                       function(i) lapply(components, function(x) x[[i]]))
  names(components) <- paste("Comp", seq_len(k), sep = ".")
  cluster <- max.col(postscaled)
  size <- if (is.null(weights)) tabulate(cluster, nbins = k)
  else tabulate(rep(cluster, weights), nbins = k)
  names(size) <- seq_len(k)
  concomitant <- FLXfillConcomitant(concomitant,
                                    postscaled[groupfirst, , drop = FALSE],
                                    weights[groupfirst])
  df <- concomitant@df(concomitant@x, k) +
    sum(sapply(components, sapply, slot, "df"))
  control@nrep <- 1
  prior <- if (is.null(weights))
    colMeans(postscaled[groupfirst, , drop = FALSE])
  else
    colSums(postscaled[groupfirst, , drop = FALSE] * weights[groupfirst]) /
    sum(weights[groupfirst])
  
  new("flexmix", model = model, prior = prior,
      posterior = list(scaled = postscaled, unscaled = postunscaled),
      weights = weights,
      iter = iter, cluster = cluster, size = size,
      logLik = llh, components = components,
      concomitant = concomitant,
      control = control, df = df, group = group,
      k = as(k, "integer"),
      converged = converged)
}

###**********************************************************
## User-facing entry point: flexmix_fast()
## Mirrors flexmix() but routes through FLXfit_fast().
## NOTE: only handles classify in {auto, weighted, hard, CEM}.
##       SEM / random fall back to the standard implementation.
###**********************************************************
flexmix_fast <- function(formula, data = list(), k = NULL, cluster = NULL,
                         model = NULL, concomitant = NULL,
                         control = NULL, weights = NULL) {
  mycall  <- match.call()
  control <- as(control, "FLXcontrol")
  if (control@classify %in% c("SEM", "random")) {
    ## fall back: same result, just slower
    z <- flexmix(formula = formula, data = data, k = k, cluster = cluster,
                 model = model, concomitant = concomitant,
                 control = control, weights = weights)
    z@call <- mycall
    return(z)
  }
  
  if (is.null(model))            model <- list(FLXMRglm())
  else if (is(model, "FLXM"))    model <- list(model)
  if (!is(concomitant, "FLXP"))  concomitant <- FLXPconstant()
  
  groups <- flexmix:::.FLXgetGrouping(formula, data)
  model  <- lapply(model, FLXcheckComponent, k, cluster)
  k <- unique(unlist(sapply(model, FLXgetK, k)))
  if (length(k) > 1) stop("number of clusters not specified correctly")
  
  model <- lapply(model, FLXgetModelmatrix, data, formula)
  
  groups$groupfirst <-
    if (length(groups$group)) flexmix:::groupFirst(groups$group)
  else rep(TRUE, FLXgetObs(model[[1]]))
  
  if (is(weights, "formula"))
    weights <- model.frame(weights, data = data, na.action = NULL)[, 1]
  if (!is.null(weights) && !identical(weights, as.integer(weights)))
    stop("only integer weights allowed")
  
  postunscaled <- flexmix:::initPosteriors(k, cluster,
                                           FLXgetObs(model[[1]]), groups)
  if (ncol(postunscaled) == 1L) concomitant <- FLXPconstant()
  concomitant <- FLXgetModelmatrix(concomitant, data = data, groups = groups)
  
  z <- FLXfit_fast(model = model, concomitant = concomitant,
                   control = control, postunscaled = postunscaled,
                   groups = groups, weights = weights)
  z@formula <- formula
  z@call    <- mycall
  z@k0      <- as.integer(k)
  z
}

###**********************************************************
## stepFlexmix_fast: drop-in style replacement that just calls flexmix_fast()
## repeatedly, keeping best per k.  Single R session, no foreach.
###**********************************************************
stepFlexmix_fast <- function(..., k, nrep = 3, verbose = TRUE,
                             drop = TRUE, unique = FALSE) {
  k <- as.integer(k)
  results <- vector("list", length(k))
  names(results) <- as.character(k)
  
  for (i in seq_along(k)) {
    best <- NULL
    for (r in seq_len(nrep)) {
      if (verbose) cat(sprintf("k = %d, rep = %d ", k[i], r))
      fit <- try(flexmix_fast(..., k = k[i]), silent = TRUE)
      if (inherits(fit, "try-error")) { if (verbose) cat("FAIL\n"); next }
      if (verbose) cat(sprintf(": logLik = %.4f\n", logLik(fit)))
      if (is.null(best) || logLik(fit) > logLik(best)) best <- fit
    }
    results[[i]] <- best
  }
  
  if (drop && length(k) == 1) return(results[[1]])
  ## Wrap as a stepFlexmix object if available, else return list
  if (isClass("stepFlexmix")) {
    new("stepFlexmix",
        models = results,
        k = k,
        nrep = as.integer(nrep),
        call = match.call())
  } else results
}