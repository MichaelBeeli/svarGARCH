#' companion_matrix
#' @rdname companion_matrix
#' @name companion_matrix
#' @title Get companion form of VAR/VARMA system
#' @description Estimates a VAR/VARMA system on a demeaned version of the provided data using the MTS package. Then, converts the MTS coefficients into companion form.
#' @export
#' @param Series A multivariate time series, where each series has a column. Can be dataframe columns or a matrix. Do not include periods or dates.
#' @param p Autoregressive order of the VAR/VARMA model. Default is 1.
#' @param q Moving-Average order of the VAR/VARMA model. Default is 0.
#' @param Demean Logical. Whether to demean the series. Default is true.
#' @return Companion Form of VAR/VARMA model, means of each series before demeaning.
#' @examples
#'  # VAR(5) of simulated data
#'   Data <- matrix(cbind(rnorm(30000, mean=0, sd=1)),  nrow=10000, ncol=3)
#'   companion <- companion_matrix(Series=Data, p=5, q=0)
#'
#' # VARMA(1,1) of returns data
#'   Returnsdata <- svarGARCH::Returns
#'   companion <- companion_matrix(Series=Returnsdata[,2:4], p=1, q=1)
#'
#'

companion_matrix <- function(
  Series,
  p=1,
  q=0,
  Demean=TRUE
  ){

  Series <- as.data.frame(Series)

  Series_Means <- sapply(Series, mean, na.rm=T)

  Series_new <- Series

  if(Demean==T) {
    Series_new <- sapply(Series, function(x) x - mean(x, na.rm=T))
  }

  if(q < 0 | p < 0) {
    stop("p and q must be positive integers")
  }
  if(q > 0) {
    Model <- MTS::VARMA(Series_new, p=p, q=q)
  } else {
    Model <- MTS::VAR(Series_new, p=p)
  }

  # Get residuals of model
  residuals <- Model$residuals

  # Initialize
  nshocks <- ncol(Model$data)
  if(is.null(Model$Theta) == T){
    Model$ARorder <- Model$order
    Model$MAorder <- 0
  }

  # AR component
  #-------------
  AR_comp   <- rbind(
    Model$Phi,
    cbind(diag(1, nshocks * (Model$ARorder - 1)),
          matrix(0, nrow = nshocks * (Model$ARorder - 1), ncol = nshocks))
  )

  # MA component
  #-------------
  if(Model$MAorder == 0){
    MA_params <- NULL
  } else {
    MA_params <- list(Model$Theta[,1:nshocks])
    if(Model$MAorder > 1){
      for(i in 2:Model$MAorder){
        MA_params[[i]] <- Model$Theta[,(nshocks + (i-2)*nshocks + 1):(nshocks + (i-1) * nshocks)]
      }
    }
  }

  return(list(
    "Residuals"    = residuals,
    "AR_companion" = AR_comp,
    "MA_companion" = MA_params,
    "Means"        = Series_Means,
    "Intercepts"   = Model$coef[1,]
  ))

}


