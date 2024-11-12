#' SVAR_Estimation_helpers
#' @rdname SVAR_Estimation_helpers
#' @name Initial_GARCH_guess
#' @title Helper functions for SVAR Estimation.
#' @description A wrapper around a function call to rugarch's ugarch fit, which extracts the fit parameters. To be used as a starting point for the numerical maximum likelihood estimation.
#' @param Residuals The residuals of the time series system.
#' @param variance_targeting Governs output formatting according to whether variance targeting is used.
#' @seealso 
#' [rugarch::ugarchspec] for specifying a GARCH model to fit.
#' [rugarch::ugarchfit] for fitting a GARCH model. 
#' @importFrom rugarch ugarchspec ugarchfit
#' @importFrom stats coef
#' @return Parameters estimated by rugarch
Initial_GARCH_guess <- function(
    Residuals,
    variance_targeting = T){

  # Specifying the GARCH(1,1) univariate
  #-------------------------------------
  spec        <- ugarchspec(
    variance.model = list(model = "fGARCH", submodel = "GARCH"),
    mean.model = list(armaOrder = c(0, 0)))

  nshocks <- ncol(Residuals)
  for (i in 1:nshocks){
    assign(paste0("estim_0_res",i),
           ugarchfit(spec, Residuals[,i], solver = "hybrid"))
  }

  params_0 <- matrix(NA, nrow = nshocks, ncol = 3)
  for(i in 1:nshocks){
    params_0[i,1] <- get(
      paste0("estim_0_res",i))@fit$coef['alpha1'] # Get alphas

    params_0[i,2] <-  get(
      paste0("estim_0_res",i))@fit$coef['beta1'] # Get betas

    params_0[i,3] <-  get(
      paste0("estim_0_res",i))@fit$coef['omega'] # Get betas
  }

  if(variance_targeting == T){
    paramlist <- list(
      "alpha" = params_0[,1],
      "beta"  = params_0[,2],
      "omega" = params_0[,3],
      "theta" =  rep(0, nshocks * (nshocks - 1) / 2) # Store alpha, beta, and zero as initial values for theta
    )
  } else {
    # Computing the initial rotation matrix
    filled_params <-
      sqrt(coef(get("estim_0_res1"))["omega"]/(1 - (coef(get("estim_0_res1"))["alpha1"] + coef(get("estim_0_res1"))["beta1"])))
    for(i in 2:nshocks){
      filled_params <- c(
        filled_params,
        rep(0, nshocks),
        sqrt(coef(get(paste0("estim_0_res",i)))["omega"]/(
          1 - (coef(get(paste0("estim_0_res",i)))["alpha1"] + coef(get(paste0("estim_0_res",i)))["beta1"])))
      )
    }

    paramlist <- list(
      "alpha" = params_0[,1],
      "beta" =  params_0[,2],
      "omega" = params_0[,3],
      "theta" = matrix(filled_params, nrow=nshocks, byrow = T) # Store alpha, beta, and diagonal entries for rotation matrix
    )
  }

  return(paramlist)



}


#' Set_GARCH_estim
#' @rdname SVAR_Estimation_helpers
#' @name Set_GARCH_estim
#' @description Various checks and organizers for the SVAR estimation function. Ensures that parameters are entered correctly.
#' @inheritParams Estimate_SVAR
#' @importFrom stats qlogis
#' @return Parameters refined for MLE
Set_GARCH_estim <- function(
    Residuals,
    variable_names=NULL,
    initialpars=NULL,
    fixedpars=NULL,
    iter = 1,
    fixvariance = TRUE,
    variance_targeting = T,
    trace_param = 0
) {


  # Get number of columns
  nvar <- ncol(Residuals)

  # Checks for initial parameter guesses
  if (!is.null(initialpars)) {
    names(initialpars) <- tolower(names(initialpars))

    if(!is.list(initialpars)) {
      stop("Parameters guess must be a list")
    }
  }

  if (is.null(initialpars)) {

    initialpars <- Initial_GARCH_guess(Residuals, variance_targeting = variance_targeting)

  } else if (!(length(setdiff(c("alpha", "beta", "omega"), names(initialpars))) == 0 & length(initialpars$alpha)==nvar & length(initialpars$beta)==nvar & length(initialpars$omega)==nvar)) {

    stop("Parameter guess must be empty or a list with 'alpha', 'beta' and 'omega' for each series. Can optionally include 'theta'")

  } else if (!(length(initialpars) %in% c(3, 4))) {

    stop("Parameter guess must be empty or a list with 'alpha', 'beta' and 'omega' for each series. Can optionally include 'theta'")

  } else if(length(initialpars)==4) {

    if(!("theta" %in% names(initialpars))) {
      stop("Parameter guess must be empty or a list with 'alpha', 'beta' and 'omega' for each series. Can optionally include 'theta'")
    }

  }else if (variance_targeting==F & !is.null(initialpars$theta) ) {

    if (!(is.matrix(initialpars$theta) & ncol(initialpars$theta)==nvar & nrow(initialpars$theta==nvar))) {
      stop("If variance targeting is not specified, non-empty guess for theta must be fully-specified rotation matrix")
    }

  } else if (variance_targeting==T  & !is.null(initialpars$theta) & !(length(initialpars$theta)==nvar*(nvar-1)/2)) {

    stop("If variance targeting is specified, non-empty guess for theta must be n * (n-1)/2 inputs for Givens rotations matrices")
  }

  # Set theta to zero if previously not set (corresponds to theta of pi for the Givens Rotation Matrices)
  if(is.null(initialpars$theta)) {

    initialpars$theta <- rep(0, nvar*(nvar-1)/2)

  }

  # Checks for fixed parameters
  if (!is.null(fixedpars)) {

    if(!is.list(initialpars)) {
      stop("Fixed parameters must be a list")
    }

    names(fixedpars) <- tolower(names(fixedpars))

    if (length(setdiff(c("alpha", "beta", "omega"), names(fixedpars))) > 1) {
      stop("Fixed parameters must be empty or either some combination alpha, beta, and omega of each series.")
    }

    if (!is.null(fixedpars$alpha) & length(fixedpars$alpha)!=nvar) {
      stop("Fixed parameters must be empty or either some combination alpha, beta, and omega of each series.")
    }

    if (!is.null(fixedpars$beta) & length(fixedpars$beta)!=nvar) {
      stop("Fixed parameters must be empty or either some combination alpha, beta, and omega of each series.")
    }

    if (!is.null(fixedpars$omega) & length(fixedpars$omega)!=nvar) {
      stop("Fixed parameters must be empty or either some combination alpha, beta, and omega of each series.")
    }

    if(length(setdiff(names(fixedpars), c("alpha", "beta", "omega")))>0) {
      stop("Valid entries for fixed parameters are alpha, beta, and omega. All named vectors in a list.")
    }

  }

  # Checking that fixed parameters agree with fixvariance argument
  if(fixvariance==TRUE & !is.null(fixedpars)) {

    omega_exp <- rep(1, nvar) - fixedpars$alpha - fixedpars$beta

    if( all( abs(omega_exp - fixedpars$omega) > 1e-12) ) {

      stop("FixVariance is specified, but fixed parameters do not satisfy that unconditional variance is equal to 1")

    }
  }

  # Propagating fixed parameters to initial guess
  if (!is.null(fixedpars$alpha)) {
    initialpars$alpha=fixedpars$alpha
  }

  if (!is.null(fixedpars$beta)) {
    initialpars$beta=fixedpars$beta
  }

  if (!is.null(fixedpars$omega)) {
    initialpars$omega=fixedpars$omega
  }

  #### Transforming parameters to make the optimization smoother
  alpha_beta <- initialpars$alpha + initialpars$beta
  alphashare <- initialpars$alpha/alpha_beta

  # Convert input parameters to their quantile outputs for the standard logistic function.This will ensure that the optimizer does not venture out of stationary territory
  # These will preserve the values of the initial guesses as they are plugged into the CDF of the logistic function, normalizing the unconditional variance to 1.

  alpha_beta <- qlogis(alpha_beta)
  alphashare <- qlogis(alphashare)

  if(!is.null(fixedpars)) {


    fixedpars$alphabeta = qlogis(fixedpars$alpha + fixedpars$beta)
    fixedpars$alphashare = qlogis(fixedpars$alpha/(fixedpars$alpha + fixedpars$beta))

  }


  # Convert matrix into list for optimx
  if (variance_targeting==F) {
    initialpars$theta <- as.vector(initialpars$theta)
  }

  # Unlist parameters into a single vector
  # 1) Alpha_Betas
  # 2) Alpha Shares
  # 3) Thetas
  params_trans <- unlist(list(alpha_beta, alphashare, initialpars$omega, initialpars$theta))


  return(list(
    "initialpars" = params_trans,
    "fixedpars"   = fixedpars
    ))


}
