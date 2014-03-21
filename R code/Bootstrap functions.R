################################################################################
##
## bootMer functions for linear calibration with grouped data
##
################################################################################

## Parallel version of ammended bootMer2 function
bootMer_parallel <- function(x, FUN, nsim = 1, seed = NULL, use.u = FALSE,
                             type = c("parametric", "semiparametric"),
                             verbose = FALSE,
                             .progress = "none", PBargs = list(),
                             parallel = c("no", "multicore", "snow"),
                             ncpus = getOption("boot.ncpus", 1L), cl = NULL)
{
  stopifnot((nsim <- as.integer(nsim[1])) > 0)
  if (.progress!="none") { ## progress bar
    pbfun <- get(paste0(.progress,"ProgressBar"))
    setpbfun <- get(paste0("set",.simpleCap(.progress),"ProgressBar"))
    pb <- do.call(pbfun,PBargs)
  }
  if (missing(parallel)) parallel <- getOption("boot.parallel", "no")
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") have_snow <- TRUE
    if (!have_mc && !have_snow) ncpus <- 1L
  }
  do_parallel <- (ncpus > 1L && (have_mc || have_snow))
  if (do_parallel & .progress!="none")
    message("progress bar disabled for parallel operations")
  
  FUN <- match.fun(FUN)
  type <- match.arg(type)
  if(!is.null(seed)) set.seed(seed)
  else if(!exists(".Random.seed", envir = .GlobalEnv))
    runif(1) # initialize the RNG if necessary
  
  mc <- match.call()
  t0 <- FUN(x)
  if (!is.numeric(t0))
    stop("bootMer_parallel currently only handles functions that return numeric vectors")
  
  mle <- list(beta = getME(x,"beta"), theta = getME(x,"theta"))
  if (isLMM(x)) mle <- c(mle,list(sigma = sigma(x)))
  ## FIXME: what about GLMMs with scale parameters??
  ## FIXME: remove prefix when incorporated in package
  
  if (type=="parametric") {
    ss <- simulate(x, nsim=nsim, use.u=use.u, na.action=na.exclude)
  } else {
    if (use.u) {
      if (isGLMM(x)) warning("semiparametric bootstrapping is questionable for GLMMs")
      ss <- replicate(nsim,fitted(x)+sample(residuals(x,"response")),
                      simplify=FALSE)
    } else {
      stop("semiparametric bootstrapping with use.u=FALSE not yet implemented")
    }
  }
  
  # define ffun as a closure containing the referenced variables 
  # in its scope to avoid explicit clusterExport statement
  # in the PSOCKcluster case 
  ffun <- local({
    FUN 
    refit 
    x 
    ss 
    verbose 
    do_parallel
    length.t0 <- length(t0)
    function(i) {
      foo <- try(FUN(refit(x,ss[[i]])),silent=TRUE)
      if(verbose) { cat(sprintf("%5d :",i)); str(foo) }
      if (!do_parallel && .progress!="none") { setpbfun(pb,i/nsim) }
      if (inherits(foo, "try-error")) rep(NA, length.t0) else foo
    }})
  
  simvec <- seq_len(nsim)
  res <- if (do_parallel) {
    if (have_mc) {
      parallel::mclapply(simvec, ffun, mc.cores = ncpus)
    } else if (have_snow) {
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
        ## explicit export of the lme4 namespace since most FUNs will probably
        ## use some of them
        parallel::clusterExport(cl, varlist=getNamespaceExports("lme4"))
        if(RNGkind()[1L] == "L'Ecuyer-CMRG")
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parLapply(cl, simvec, ffun)
        parallel::stopCluster(cl)
        res
      } else parallel::parLapply(cl, simvec, ffun)
    }
  } else lapply(simvec, ffun)
  
  t.star <- do.call(cbind,res)
  rownames(t.star) <- names(t0)
  if ((numFail <- sum(apply(is.na(t.star),2,all)))>0) {
    warning("some bootstrap runs failed (",numFail,"/",nsim,")")
  }
  ## boot() ends with the equivalent of
  ## structure(list(t0 = t0, t = t.star, R = R, data = data, seed = seed,
  ##  	      statistic = statistic, sim = sim, call = call,
  ##		      ran.gen = ran.gen, mle = mle),
  ##		 class = "boot")
  s <- structure(list(t0 = t0, t = t(t.star), R = nsim, data = x@frame,
                      seed = .Random.seed,
                      statistic = FUN, sim = "parametric", call = mc,
                      ## these two are dummies
                      ran.gen = "simulate(<lmerMod>, 1, *)", mle = mle),
                 class = "boot")
  attr(s,"bootFail") <- numFail
  s
} ## {bootMer_parallel}


## Ammended bootMer function 
bootMer2 <- function(x, FUN, FUN0, nsim = 1, seed = NULL, use.u = FALSE,
                     type=c("parametric","semiparametric"),
                     verbose = FALSE,
                     .progress="none", PBargs=list(),
                     parallel = c("no", "multicore", "snow"),
                     ncpus = getOption("boot.ncpus", 1L), cl = NULL)
{
  stopifnot((nsim <- as.integer(nsim[1])) > 0)
  if (.progress!="none") { ## progress bar
    pbfun <- get(paste0(.progress,"ProgressBar"))
    setpbfun <- get(paste0("set",.simpleCap(.progress),"ProgressBar"))
    pb <- do.call(pbfun,PBargs)
  }
  if (missing(parallel)) parallel <- getOption("boot.parallel", "no")
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") have_snow <- TRUE
    if (!have_mc && !have_snow) ncpus <- 1L
  }
  do_parallel <- (ncpus > 1L && (have_mc || have_snow))
  if (do_parallel & .progress!="none")
    message("progress bar disabled for parallel operations")
  
  FUN <- match.fun(FUN)
  FUN0 <- match.fun(FUN0)
  type <- match.arg(type)
  if(!is.null(seed)) set.seed(seed)
  else if(!exists(".Random.seed", envir = .GlobalEnv))
    runif(1) # initialize the RNG if necessary
  
  mc <- match.call()
  t0 <- FUN0(x)
  if (!is.numeric(t0))
    stop("bootMer currently only handles functions that return numeric vectors")
  
  mle <- list(beta = getME(x,"beta"), theta = getME(x,"theta"))
  if (isLMM(x)) mle <- c(mle,list(sigma = sigma(x)))
  ## FIXME: what about GLMMs with scale parameters??
  ## FIXME: remove prefix when incorporated in package
  
  if (type=="parametric") {
    ss <- simulate(x, nsim=nsim, use.u=use.u)
  } else {
    if (use.u) {
      if (isGLMM(x)) warning("semiparametric bootstrapping is questionable for GLMMs")
      ss <- replicate(nsim,fitted(x)+sample(residuals(x,"response")),
                      simplify=FALSE)
    } else {
      stop("semiparametric bootstrapping with use.u=FALSE not yet implemented")
    }
  }
  
  # define ffun as a closure containing the referenced variables 
  # in its scope to avoid explicit clusterExport statement
  # in the PSOCKcluster case 
  ffun <- local({
    FUN 
    refit 
    x 
    ss 
    verbose 
    do_parallel
    length.t0 <- length(t0)
    function(i) {
      foo <- try(FUN(refit(x,ss[[i]])),silent=TRUE)
      if(verbose) { cat(sprintf("%5d :",i)); str(foo) }
      if (!do_parallel && .progress!="none") { setpbfun(pb,i/nsim) }
      if (inherits(foo, "try-error")) rep(NA, length.t0) else foo
    }})
  
  simvec <- seq_len(nsim)
  res <- if (do_parallel) {
    if (have_mc) {
      parallel::mclapply(simvec, ffun, mc.cores = ncpus)
    } else if (have_snow) {
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
        ## explicit export of the lme4 namespace since most FUNs will probably
        ## use some of them
        parallel::clusterExport(cl, varlist=getNamespaceExports("lme4"))
        if(RNGkind()[1L] == "L'Ecuyer-CMRG")
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parLapply(cl, simvec, ffun)
        parallel::stopCluster(cl)
        res
      } else parallel::parLapply(cl, simvec, ffun)
    }
  } else lapply(simvec, ffun)
  
  t.star <- do.call(cbind,res)
  rownames(t.star) <- names(t0)
  if ((numFail <- sum(apply(is.na(t.star),2,all)))>0) {
    warning("some bootstrap runs failed (",numFail,"/",nsim,")")
  }
  ## boot() ends with the equivalent of
  ## structure(list(t0 = t0, t = t.star, R = R, data = data, seed = seed,
  ##		      statistic = statistic, sim = sim, call = call,
  ##		      ran.gen = ran.gen, mle = mle),
  ##		 class = "boot")
  s <- structure(list(t0 = t0, t = t(t.star), R = nsim, data = x@frame,
                      seed = .Random.seed,
                      statistic = FUN, sim = "parametric", call = mc,
                      ## these two are dummies
                      ran.gen = "simulate(<lmerMod>, 1, *)", mle = mle),
                 class = "boot")
  attr(s,"bootFail") <- numFail
  s
} ## {bootMer}

## Parallel version of ammended bootMer2 function
bootMer2_parallel <- function(x, FUN, FUN0, nsim = 1, seed = NULL, use.u = FALSE,
                              type=c("parametric","semiparametric"),
                              verbose = FALSE,
                              .progress="none", PBargs=list(),
                              parallel = c("no", "multicore", "snow"),
                              ncpus = getOption("boot.ncpus", 1L), cl = NULL)
{
  stopifnot((nsim <- as.integer(nsim[1])) > 0)
  if (.progress!="none") { ## progress bar
    pbfun <- get(paste0(.progress,"ProgressBar"))
    setpbfun <- get(paste0("set",.simpleCap(.progress),"ProgressBar"))
    pb <- do.call(pbfun,PBargs)
  }
  if (missing(parallel)) parallel <- getOption("boot.parallel", "no")
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") have_snow <- TRUE
    if (!have_mc && !have_snow) ncpus <- 1L
  }
  do_parallel <- (ncpus > 1L && (have_mc || have_snow))
  if (do_parallel & .progress!="none")
    message("progress bar disabled for parallel operations")
  
  FUN <- match.fun(FUN)
  type <- match.arg(type)
  if(!is.null(seed)) set.seed(seed)
  else if(!exists(".Random.seed", envir = .GlobalEnv))
    runif(1) # initialize the RNG if necessary
  
  mc <- match.call()
  t0 <- FUN0(x)
  if (!is.numeric(t0))
    stop("bootMer2_parallel currently only handles functions that return numeric vectors")
  
  mle <- list(beta = getME(x,"beta"), theta = getME(x,"theta"))
  if (isLMM(x)) mle <- c(mle,list(sigma = sigma(x)))
  ## FIXME: what about GLMMs with scale parameters??
  ## FIXME: remove prefix when incorporated in package
  
  if (type=="parametric") {
    ss <- simulate(x, nsim=nsim, use.u=use.u, na.action=na.exclude)
  } else {
    if (use.u) {
      if (isGLMM(x)) warning("semiparametric bootstrapping is questionable for GLMMs")
      ss <- replicate(nsim,fitted(x)+sample(residuals(x,"response")),
                      simplify=FALSE)
    } else {
      stop("semiparametric bootstrapping with use.u=FALSE not yet implemented")
    }
  }
  
  # define ffun as a closure containing the referenced variables 
  # in its scope to avoid explicit clusterExport statement
  # in the PSOCKcluster case 
  ffun <- local({
    FUN 
    refit 
    x 
    ss 
    verbose 
    do_parallel
    length.t0 <- length(t0)
    function(i) {
      foo <- try(FUN(refit(x,ss[[i]])),silent=TRUE)
      if(verbose) { cat(sprintf("%5d :",i)); str(foo) }
      if (!do_parallel && .progress!="none") { setpbfun(pb,i/nsim) }
      if (inherits(foo, "try-error")) rep(NA, length.t0) else foo
    }})
  
  simvec <- seq_len(nsim)
  res <- if (do_parallel) {
    if (have_mc) {
      parallel::mclapply(simvec, ffun, mc.cores = ncpus)
    } else if (have_snow) {
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
        ## explicit export of the lme4 namespace since most FUNs will probably
        ## use some of them
        parallel::clusterExport(cl, varlist=getNamespaceExports("lme4"))
        if(RNGkind()[1L] == "L'Ecuyer-CMRG")
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parLapply(cl, simvec, ffun)
        parallel::stopCluster(cl)
        res
      } else parallel::parLapply(cl, simvec, ffun)
    }
  } else lapply(simvec, ffun)
  
  t.star <- do.call(cbind,res)
  rownames(t.star) <- names(t0)
  if ((numFail <- sum(apply(is.na(t.star),2,all)))>0) {
    warning("some bootstrap runs failed (",numFail,"/",nsim,")")
  }
  ## boot() ends with the equivalent of
  ## structure(list(t0 = t0, t = t.star, R = R, data = data, seed = seed,
  ##		      statistic = statistic, sim = sim, call = call,
  ##		      ran.gen = ran.gen, mle = mle),
  ##		 class = "boot")
  s <- structure(list(t0 = t0, t = t(t.star), R = nsim, data = x@frame,
                      seed = .Random.seed,
                      statistic = FUN, sim = "parametric", call = mc,
                      ## these two are dummies
                      ran.gen = "simulate(<lmerMod>, 1, *)", mle = mle),
                 class = "boot")
  attr(s,"bootFail") <- numFail
  s
} ## {bootMer2_parallel}

##' @S3method as.data.frame boot
as.data.frame.boot <- function(x,...) {
  as.data.frame(x$t)
}

