#' SVAR_Estimation ----
#' @rdname SVAR_Estimation
#' @name Estimate_SVAR
#' @title Estimate Structural VAR by GARCH Errors
#' @export
#' @description Estimates a structural VAR by Maximum likelihood estimation. Identifies as many structural shocks as series using a GARCH(1,1) model for conditionally heteroskedastic residuals.
#' Must provide the residuals for each of the series along with their names. Additional entries are optional.
#' @param Residuals The residuals of the VAR/VARMA model.
#' @param variable_names Names of the series. Must follow the same ordering as the residuals.
#' @param initialpars Optional initial guess for the optimal GARCH(1,1) parameters. If specified, should be a list with an alpha, beta, and omega for each series. Specifying initial guesses for theta is more advanced and optional.
#' @param fixedpars Optional, GARCH(1,1) Parameters that the user wishes to fix. If specified, must be a list with some combination of the alpha, beta, and omega of each series. alpha, beta, and omega must be named vectors in a list.
#' @param iter Number of iterations for the Maximum likelihood estimation. Default is 1.
#' @param fixvariance Logical, whether to normalize the unconditional variance of each structural shock to 1. Default is true.
#' @param variance_targeting Whether the asymptotic variance should be constrained to the unconditional variance. Default is true.
#' @param trace_param Trace code for "optimx" package. See `optimx::optimx`
#' @seealso
#' [optimx::optimx]
#' [rugarch::ugarchfit]
#' [steadyICA::theta2W]
#' @details
#' @section Initial Parameters:
#' `initialpars` are an optional entry for the initial guess for the true parameters. If no guess is provided, the initial GARCH(1, 1) parameters will be taken from `rugarch::ugarchfit`. The "theta" parameters represent the general set of parameters that characterize the rotation matrix for the structural residuals.
#'
#' "alpha", "beta", and "omega" should be named vectors in a list, and should be entered as the actual values for the GARCH parameters.The parameters passed to the log-likelihood function will be transformed to ensure stationarity.
#'
#' "theta" is an optional additional named vector in the list. Let n be the number of series. If the variance targeting option is specified, theta should be  n * (n - 1) / 2 Givens Rotation Matrices inputs (theta). If variance targeting is not specified, then theta should be a column-wise vector of initial guesses for the (n x n) rotation matrix.
#'
#' @section Fixed Parameters:
#' `fixedpars` is an optional input that allows the user to fix combinations of GARCH(1, 1) parameters in the estimation. Not all parameters must be entered, but if one parameter is entered for one series, it must be entered for all series.
#'
#' @section Variance Targeting:
#' `variance_targeting`
#'
#' @references
#' Killian, L., and Lütkepohl, H. (2017) *Structural Vector Autoregressive Analysis*. Cambridge University Press.
#' ISBN: 978-1108164818
#'
#' Lanne, M., and Saikkonen, P. (2007) "A Multivariate Generalized Orthogonal Factor GARCH Model".
#' *Journal of Business Economics and Statistics, American Statistical Association*, 25 (1): 61–75.
#' doi:10.1198/073500106000000404.
#'
#' Hamilton, J. (1994) *Time Series Analysis*. Princeton University Press.
#' ISBN: 9780691042893
#'
#' @return Estimated structural VAR/GARCH parameters, rotation matrix, structural residuals, structural shock variance, and the log-likelihood series.
#' @examples
#'
#' Returnsdata <- svarGARCH::Returns
#' companion <- companion_matrix(Returnsdata[,2:4])
#'
#' # Default estimation method, no fixed parameters, no initial guess for parameters
#' Estimation_results <- Estimate_SVAR(
#'    Residuals                = companion$Residuals,
#'    variable_names           = names(Returnsdata[,2:4]),
#'    fixvariance              = TRUE,
#'    variance_targeting       = TRUE
#' )
#'
#'
Estimate_SVAR <- function(
    Residuals,
    variable_names=NULL,
    initialpars=NULL,
    fixedpars=NULL,
    iter = 1,
    fixvariance = TRUE,
    variance_targeting = T,
    trace_param = 0
){


  #### Run input checks and transform parameters to proper format for estimation
  pars <- Set_GARCH_estim(
    Residuals=Residuals,
    variable_names=variable_names,
    initialpars=initialpars,
    fixedpars=fixedpars,
    iter=iter,
    fixvariance=fixvariance,
    variance_targeting=variance_targeting,
    trace_param=trace_param
  )

  #### Maximize the log-likelihood
  Estimation_results <- Maximize_loglik(
    vec_params = pars$initialpars,
    Residuals=Residuals,
    variable_names=variable_names,
    fixedpars = pars$fixedpars,
    fixvariance = fixvariance,
    variance_targeting = variance_targeting,
    iter=iter,
    trace_param=trace_param
  )

  #### Return results
  return(Estimation_results)

}
