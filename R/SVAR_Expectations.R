#' SVAR_Expectations ----
#' @rdname SVAR_Expectations
#' @name SVAR_Expectations
#' @title SVAR forecasts
#' @export
#' @description Compute optimal h-step ahead forecasts for SVAR system at each date. Must provide companion form components, true series, structural residuals, rotation matrix, dates, and a numeric vector of forecast horizons. Returns an array of the optimal forecasts at each date for the specified horizons.
#' @param  AR_companion The autoregressive component of the companion matrix.
#' @param MA_companion The moving average component of the companion matrix.
#' @param Series The underlying data of the SVAR (excluding dates/periods).
#' @param Structural_residuals Structural residuals identified by SVAR.
#' @param rotation_matrix The rotation matrix identified by SVAR.
#' @param Dates_Series The dates to compute forecasts for.
#' @param variable_names Names of the series. Must follow the same ordering as the residuals.
#' @param linear_comb Optional, a named list of new variables that are linear combinations of the original series. The names should be entred as the list names, and the scalars to generate the linear combination should be the elements.
#' @param horizons Forecast horizons to compute.
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
#' Forecasts <- SVAR_Expectations(
#'    AR_companion             =companion$AR_companion,
#'    MA_companion             =companion$MA_companion,
#'    Series                   =Returns[,2:4],
#'    Structural_residuals     =Estimation_results$Struct_resid,
#'    rotation_matrix          =Estimation_results$rotation_matrix,
#'    Dates_Series             =Returns$Date,
#'    variable_names           =c("SPX", "BRK-A", "HPQ"),
#'    linear_comb              =list("BRK-SPX"=c(-1,1,0)),
#'    horizons                 =c(5, 21, 63)
#')
#'
SVAR_Expectations <- function(
    AR_companion,
    MA_companion = NULL,
    Series,
    Intercept = NULL,
    Structural_residuals,
    rotation_matrix,
    Dates_Series,
    variable_names,
    linear_comb,
    horizons
    ){

  # Initialize the Objects
  #-----------------------
  if(is.null(Intercept)) {
    Intercept <- rep(0, ncol(Series))
  }

  nshocks   <- ncol(Structural_residuals)
  totvars <- nshocks + dplyr::coalesce(length(linear_comb), 0)
  allvarnames   <- c(variable_names, names(linear_comb))

  # Calculating the AR and MA orders
  AR_order  <- ncol(AR_companion)/nshocks
  MA_order  <- length(MA_companion)

  # Taking size of the extended system
  nextended <- (AR_order + MA_order) * nshocks

  # Initializing the AR matrix and the rotation matrix in the companion forms
  Total_AR              <- AR_companion
  Total_rotation_matrix <- rbind(
    rotation_matrix,
    matrix(0, nrow = (AR_order - 1) * nshocks, ncol = nshocks)
  )


  # Construct the total AR matrix in case of MA terms
  #- - - - - - - - - - - - - - - - - - - - - - - - -
  if(MA_order > 0){
    Total_AR <- rbind(
      AR_companion,
      matrix(0, nrow = MA_order * nshocks, ncol = ncol(AR_companion))
    )
    # Filling the top-right corner
    top_right <- NULL
    for(i in 1:MA_order){
      top_right <- cbind(top_right,
                         MA_companion[[i]] %*% rotation_matrix)
    }
    top_right <- rbind(
      top_right,
      matrix(0, nrow = AR_order * nshocks, ncol = MA_order * nshocks)
    )

    # Filling up the rotation matrix
    Total_rotation_matrix <- rbind(
      Total_rotation_matrix,
      diag(1, nshocks)
    )

    # Filling up the total AR matrix of the extended model
    if(MA_order == 1){
      Total_AR <- cbind(Total_AR, top_right)
    } else {
      bottom_right <- cbind(
        diag(1, (MA_order - 1) * nshocks),
        matrix(0, nrow = (MA_order - 1) * nshocks, ncol = nshocks)
      )
      Total_AR <- cbind(
        Total_AR,
        rbind(top_right, bottom_right)
      )

      # Filling up the rotation matrix again
      Total_rotation_matrix <- rbind(
        Total_rotation_matrix,
        matrix(0, nrow = (MA_order - 1) * nshocks, ncol = nshocks)
      )
    }
  }

  if (!is.null(linear_comb)) {
    linear_comb_ext <- as.matrix(cbind(
      rbind(
        diag(1, nrow = ncol(rotation_matrix), names = TRUE),
        matrix(unlist(linear_comb), nrow = length(linear_comb), byrow = T)
      ),
      matrix(0, nrow=totvars, ncol(Total_AR) - nshocks)
    ))
  } else {
    linear_comb_ext <-  as.matrix(cbind(
      rbind(
        diag(1, nrow = ncol(rotation_matrix), names = TRUE)
      ),
      matrix(0, nrow=totvars, ncol(Total_AR) - nshocks)
    ))
  }

  dimnames(linear_comb_ext) <- list(c(allvarnames), c(variable_names, rep("", ncol(Total_AR) - nshocks)))

  # Obtaining the forecasts
  #========================
  number_periods    <- nrow(Series)
  Total_forecasts <- array(
    NA, dim = c(number_periods, nshocks * (AR_order + MA_order), max(horizons)),
    dimnames = list(
      as.character(Dates_Series),
      paste0("V",c(1:(nshocks * (AR_order + MA_order)))),
      paste0("H=",c(1:max(horizons))))
  )

  Forecasts_comb  <- array(
    NA, dim = c(number_periods, totvars, max(horizons)),
    dimnames = list(
      as.character(Dates_Series),
      allvarnames,
      paste0("H=",c(1:max(horizons))))
  )

  # Construct matrix of lagged factors
  factors_forecast  <- matrix(
    NA, nrow = nshocks * (AR_order + MA_order), ncol = number_periods
  )

  length_suppressed <- max(AR_order, MA_order)

  Intercept_vec <- c(Intercept, rep(0, 3*(length_suppressed-1)))

  Intercept_matrix <- matrix(rep(Intercept_vec, nrow(Series)), ncol=nshocks * (AR_order + MA_order), nrow = nrow(Series), byrow = T)

  Filled_structural_residuals <- rbind(
    matrix(NA, nrow = length_suppressed, ncol = nshocks),
    Structural_residuals
  )
  # Setting up forecast components
  for(i in length_suppressed:nrow(Series)){
  # Adding residuals if MA component
    if(MA_order>0) {
      factors_forecast[,i] <- c(
        c(t(Series[i:(i - AR_order + 1),])),
        c(t(Filled_structural_residuals[i:(i - MA_order + 1),]))
      )
    } else { # If AR component only, do not add residuals
      factors_forecast[,i] <- c(
        c(t(Series[i:(i - AR_order + 1),])))
    }
  }

  # Initialize first forecast
  Total_forecasts[,,1] <- Intercept_matrix + t((Total_AR %*% factors_forecast))
  Forecasts_comb[,,1]  <- t(linear_comb_ext %*% t(Intercept_matrix + t(Total_AR %*% factors_forecast)))
  # Calculate all other forecasts
  for(i in 2:max(horizons)){
    Total_forecasts[,,i] <-  Intercept_matrix + t((Total_AR %*% t(Total_forecasts[,,i-1])))
    Forecasts_comb[,,i]  <- t(linear_comb_ext %*% t(Total_forecasts[,,i]))
  }


  # Finalize the Object
  #====================
  Forecasts_comb <- Forecasts_comb[, , horizons]

  return(Forecasts_comb)
}

