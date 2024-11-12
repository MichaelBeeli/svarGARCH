#' SVAR_Log_Likelihood
#' @rdname SVAR_Log_Likelihood
#' @name Maximize_loglik
#' @title For advanced use only. Calculate and Maximizing Log-likelihood function of Structural VAR with GARCH(1,1) errors. Not intended to be called directly.
#' @description A wrapper function to maximize the log-likelihood function within Estimate_SVAR.
#' @param vec_params Initial guesses for parameters for the log-likelihood function. Must be a vector with n values for the Alpha+Beta, followed by n values of Alpha/(Alpha+Beta), followed by n values of Omega. If variance targeting is true, these must be followed by n(n-1)/2 values for theta. If variance targeting is false, then they must be followed by the column-wise entries of the rotation matrix.
#' @param Residuals The residuals of the VAR/VARMA model.
#' @param variable_names Names of the series. Must follow the same ordering as the residuals.
#' @param fixedpars  Parameters that the user wishes to fix.
#' @param fixvariance Logical, whether to constrain the unconditional variance to 1.
#' @param variance_targeting Whether the asymptotic variance should be constrained to the unconditional variance.
#' @param iter Number of iterations. Default is 1.
#' @param trace_param Trace code for "optimx" package. See "optimx".
#' @importFrom optimx optimx
#' @return GARCH/SVAR Parameters that maximize the log-likelihood.
Maximize_loglik <- function(
    vec_params,
    Residuals,
    variable_names,
    fixedpars=NULL,
    fixvariance,
    variance_targeting = T,
    iter = 1,
    trace_param = 0
    ){


  # Optimize
  #---------
  for(j in 1:iter){

    # Start With nlminb
    maxim_lik_VARMA <- suppressWarnings(optimx(
      par = vec_params,
      fn = GARCH_Loglik_SVAR,
      Residuals = Residuals,
      variable_names = variable_names,
      fixedpars = fixedpars,
      fixvariance = fixvariance,
      variance_targeting = variance_targeting,
      method = "nlminb",
      control = list(trace = trace_param, kkt = F)))

    vec_params <- as.numeric(coef(maxim_lik_VARMA))
    fixedpars <- fixedpars

    # Do Nelder-Mead
    maxim_lik_VARMA <-  suppressWarnings(optimx(
      par = vec_params,
      fn = GARCH_Loglik_SVAR,
      Residuals = Residuals,
      variable_names = variable_names,
      fixedpars = fixedpars,
      fixvariance = fixvariance,
      variance_targeting = variance_targeting,
      method = "Nelder-Mead",
      itnmax = 5000,
      control = list(trace = trace_param, kkt = F)))

    vec_params <- as.numeric(coef(maxim_lik_VARMA))
    fixedpars <- fixedpars
  }

  # Send Back Results
  #-------------------
  estim_results <-
    GARCH_Loglik_SVAR(
      vec_params = vec_params,
      Residuals = Residuals,
      variable_names = variable_names,
      fixedpars = fixedpars,
      fixvariance = fixvariance,
      foroptim = F,
      variance_targeting = variance_targeting
    )

  return(estim_results)

}



#' @rdname SVAR_Log_Likelihood
#' @name GARCH_Loglik_SVAR
#' @description For advanced use only. Builds the rotation matrix, structural residuals and their variance, and calculates the log-likelihood for a given set of parameters.
#' @export
#' @param vec_params The given parameters for the log-likelihood function. Must be a vector with n values for the Alpha+Beta, followed by n values of Alpha/(Alpha+Beta), followed by n values of Omega. If variance targeting is true, these must be followed by n(n-1)/2 values for theta. If variance targeting is false, then they must be followed by the column-wise entries of the rotation matrix.
#' @param Residuals The residuals of the VAR/VARMA model.
#' @param variable_names Names of the series. Must follow the same ordering as the residuals.
#' @param fixedpars Parameters that the user wishes to fix. Has no use when being called directly.
#' @param fixvariance Whether the unconditional variance should be normalized to 1.
#' @param variance_targeting Whether the asymptotic variance should be constrained to the unconditional variance.
#' @param foroptim Logical value, whether the log-likelihood is being calculated as part of an optimization problem. Should be false when the function is being called directly.
#' @param qtransformed To be set false when calling the function directly and the GARCH parameters (Alpha+Beta and the Alpha/(Alpha+Beta)) are entered as their actual values.
#' @importFrom stats var plogis
#' @importFrom steadyICA theta2W
#' @seealso 
#' [steadyICA::theta2W]
#' @references 
#' Killian, L., and Lütkepohl, H. (2017) *Structural Vector Autoregressive Analysis*. Cambridge University Press. 
#' ISBN: 978-1108164818
#' 
#' Lanne, M., and Saikkonen, P. (2007) "A Multivariate Generalized Orthogonal Factor GARCH Model". 
#' *Journal of Business Economics and Statistics, American Statistical Association*, 25 (1): 61–75. doi:10.1198/073500106000000404. 
#' @return Log-likelihood and SVAR/GARCH parameters.
#' @examples
#' Returnsdata <- svarGARCH::Returns
#' companion <- companion_matrix(Returnsdata[,2:4])
#'
#' # Calling the function directly.
#' # Gets the rotation matrix and log-likelihood for a set of specified parameters
#' GARCH_Loglik_SVAR(
#'    vec_params            = c(.9, .98, .95, .2, .1, .15, .1, .02, .05, rep(0, 3*2/2)),
#'    Residuals             = companion$Residuals,
#'    variable_names        = names(Returnsdata[,2:4]),
#'    fixvariance           = TRUE,
#'    variance_targeting    = TRUE,
#'    foroptim              = FALSE,
#'    qtransformed          = FALSE
#' )
#'
#'
GARCH_Loglik_SVAR <-
  function(
    vec_params,
    Residuals,
    variable_names,
    fixedpars=NULL,
    fixvariance=TRUE,
    variance_targeting = T,
    foroptim = T,
    qtransformed=T
    ){

    # Defining the size of the model
    #-------------------------------
    # - Residuals is (T x n)
    nvar        <- ncol(Residuals)
    nvar_plus <- nvar +1
    nvar2 <- 2 * nvar
    nvar2_plus <- nvar2 +1
    nvar_3 <- 3 * nvar
    params_length <- length(vec_params)

    # If alpha and beta parameters are not already in quantile output form, convert
    if(qtransformed==FALSE) {

      vec_params[1:nvar] <- qlogis(vec_params[1:nvar])
      vec_params[nvar_plus:nvar2] <- qlogis(vec_params[nvar_plus:nvar2])

    }


    ## Pass the iterated parameters through
    alpha_beta      <- plogis(vec_params[1:nvar])
    alpha           <- plogis(vec_params[nvar_plus:nvar2]) * alpha_beta
    beta            <- alpha_beta - alpha
    omega           <- vec_params[nvar2_plus:nvar_3]
    theta_index     <- nvar * 3

    # If parameters are fixed, replace with fixed parameters
    if(!is.null(fixedpars)) {

      alpha_beta <- plogis(fixedpars$alphabeta)
      alpha      <- plogis(fixedpars$alphashare) * alpha_beta
      beta       <- alpha_beta - alpha
      omega      <- fixedpars$omega

    }

    # If unconditional variance is normalized to 1, then adjust omega accordingly
    if(fixvariance==TRUE) {

      # Normalizing omega such that unconditional variance of 1
      omega           <- (1 - alpha - beta)

    }


    # Estimating the rotation matrix
    #```````````````````````````````
    if (variance_targeting == F){

      # Fill in the rotation matrix unconstrained
      rotation_matrix <- matrix(
        vec_params[(theta_index + 1):(theta_index + nvar^2)],
        nrow = nvar, ncol = nvar)

    } else if (variance_targeting == T){

      # Fill in the Givens rotation matrix
      vector_rotation <- 2 * pi * plogis(
        vec_params[(theta_index + 1):(theta_index + nvar * (nvar - 1)/2)])
      Givens_mat      <- steadyICA::theta2W(vector_rotation)
      # Rotating the Choleski
      rotation_matrix <- t(chol(var(Residuals, na.rm = T))) %*%
        Givens_mat

    }

    # Initialize the series
    #----------------------
    # Log-likelihood series
    loglik          <- rep(NA, nrow(Residuals))
    # 3 Series of Variance
    variance_series <- NA * Residuals
    # Structural residuals
    Struct_resid    <- t(solve(rotation_matrix) %*% t(Residuals))
    # Last variance object
    last_variance   <- omega/(1-alpha-beta)
    # Last structural residual
    last_struct     <- rep(0,nvar)


    # Loop to form the log-likelihood function
    #----------------------------------------
    for(i in 1:nrow(Residuals)){

      # Calculate the structural residuals variance
      variance_series[i,]   <- omega + alpha * last_struct^2 + beta * last_variance

      # Calculate the log-likelihood
      Condi_var <- rotation_matrix %*%
        diag(c(last_variance)) %*%
        t(rotation_matrix)

      loglik[i] <- -1/2 * (
        nvar * log(2 * pi) +
          log(det(Condi_var)) +
          sum(Struct_resid[i,]^2/last_variance)
      )

      # Update objects for the loop
      last_variance <- variance_series[i,]
      last_struct   <- Struct_resid[i,]
    }

    # return results
    if(foroptim==T){
      return(-sum(loglik))
    } else{
      return(
        list(
          "loglik"          = loglik,
          "variance_series" = variance_series,
          "omega"           = omega,
          "alpha"           = alpha,
          "beta"            = beta,
          "rotation_matrix" = array(rotation_matrix, c(nvar, nvar),
                                    dimnames = list(
                                      variable_names,
                                      paste("shock", 1:nvar)
                                    )),
          "Struct_resid"    = Struct_resid,
          "vec_params"      = vec_params,
          "fixedpars"       = fixedpars
        )
      )
    }

  }
