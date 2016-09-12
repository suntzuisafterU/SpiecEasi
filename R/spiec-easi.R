#' SPIEC-EASI
#'
#' Run the whole pipeline
#'
#' @param data OTU count matrix with samples on rows and features/OTUs in columns or a phyloseq object with an OTU count table
#' @param method estimation method to use as a character string. Currently either 'glasso' or 'mb' (meinshausen-buhlmann)
#' @param sel.criterion character string specifying criterion/method for model selection accepts 'stars' [default] or 'gstars'
#' @param pulsar.select flag to use pulsar for model selection
#' @param pulsar.params named list of arguments to pass to pulsar. See details.
#' @param icov.select defunct. can be set for backwards compatability, but ignored if pulsar.select and/or pulsar.params are given
#' @param icov.select.params defunct. can be set for backwards compatability, but ignored if pulsar.select and/or pulsar.params are given
#' @param ... further arguments
#' @export
spiec.easi <- function(data, ...) {
  UseMethod('spiec.easi', data)
}

#' @rdname spiec.easi
#' @import phyloseq
spiec.easi.phyloseq <- function(data, ...) {
  OTU <- otu_table(data)@.Data
  if (otu_table(data)@taxa_are_rows) OTU <- t(OTU)
  spiec.easi.default(OTU, ...)
}


#' @rdname spiec.easi
#' @importFrom pulsar pulsar getLamPath getMaxCov refit get.opt.index opt.index
#' @export
spiec.easi.default <- function(data, method='glasso', sel.criterion='stars', verbose=TRUE,
                               pulsar.select=TRUE, pulsar.params=list(),
                               icov.select=pulsar.select, icov.select.params=pulsar.params, ...) {

  args <- list(...)
  if (verbose) message("Applying data transformations...")

  switch(method,
        glasso = { estFun <- "sparseiCov" ;  args$method <- method
                   X <- t(clr(data+1, 1)) ; maxlam <- getMaxCov(cor(X)) },

        mb     = { estFun <- "sparseiCov" ;  args$method <- method
                   X <- t(clr(data+1, 1)) ; maxlam <- getMaxCov(cor(X)) },

        slr    = { estFun <- "sparseLowRankiCov" ;
                   X <- t(clr(data+1, 1)) ; maxlam <- getMaxCov(cor(X)) },

        coat   = { estFun <- "coat" ; X <- t(clr(data+1, 1)) ; maxlam <- getMaxCov(X) },

        ising  = { estFun <- "neighborhood.net" ; args$method <- method ; 
                   X <- sign(data) ; maxlam <- max(abs(t(scale(X)) %*% X)) / nrow(X) },

        poisson= { estFun <- "neighborhood.net" ; args$method <- method
                   X <- data ; maxlam <- max(abs(t(scale(X)) %*% X)) / nrow(X) },

        loglin = { estFun <- "neighborhood.net" ; args$method <- method
                   X <- data ; maxlam <- max(abs(t(scale(X)) %*% X)) / nrow(X) }
    )

  if (is.null(args[[ "lambda" ]])) {
    if (is.null(args[[ "lambda.min.ratio" ]])) args$lambda.min.ratio <- 1e-3
    if (is.null(args[[ "nlambda" ]])) args$nlambda <- 20
    args$lambda <- getLamPath(maxlam, maxlam*args$lambda.min.ratio, args$nlambda, log=TRUE)
    args$lambda.min.ratio <- NULL ; args$nlambda <- NULL
  }

  ocall <- match.call(expand.dots=FALSE)
  ## if pulsar options are not specified, check for deprecated icov.select options are
  if (is.null(ocall[[ "pulsar.select" ]]) && is.null(ocall[[ "pulsar.params" ]])) {
      pulsar.select <- icov.select
      pulsar.params <- icov.select.params
  }

  pulsar.params$criterion <- 
    switch(sel.criterion,
           stars = "stars",
          bstars = "stars",
          gstars = c("stars", "gcd"),
          stop("Unknown selection criterion"))

  ## process pulsar.params defaults
  if (sel.criterion %in% c("bstars", "gstars"))
    pulsar.params$lb.stars <- pulsar.params$ub.stars <- TRUE
  if (is.null(pulsar.params[[ "thresh" ]])) pulsar.params$thresh <- 0.05

  call <- quote(pulsar(data=X, fun=match.fun(estFun), fargs=args))
  call <- do.call('update', c(pulsar.params, list(object=list(call=call), evaluate=FALSE)))

  if (pulsar.select) {
    if (verbose) message("Selecting model with pulsar using ", sel.criterion, "...")
    est <- eval(call, environment())
    if (sel.criterion == "gstars") 
        opt.index <- pulsar::opt.index(est) <- get.opt.index(est, 'gcd')
    else
        opt.index <- opt.index(est, 'stars')
  } else
    est <- structure(list(call=call, envir=environment()), class='pulsar')

  if (verbose) message("Fitting final estimate with ", method, "...")
  suppressWarnings(
  fit <- refit(est)
  )
  if (pulsar.select) {
    fit$select <- list(
     merge     = est$stars$merge,
     summary   = est$stars$summary,
     opt.index = opt.index,
     lb.index  = est$stars$lb.index,
     ub.index  = est$stars$ub.index
    )
    fit$refit <- if (sel.criterion=="gstars") fit$refit$gcd else fit$refit$stars
    attr(fit$refit, 'names') <- method
  }
  fit$lambda <- args$lambda
  fit$fun    <- call(estFun)[[1]]
  if (verbose) message('done')

  return(fit)
}
