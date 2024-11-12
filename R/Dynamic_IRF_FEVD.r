#' Dynamic_IRF_FEVD
#' @rdname Dynamic_IRF_FEVD
#' @description Computes Impulse Responses and Forecast Error Variance Decompositions up to specified length for each date in range provided.
#' @title Dynamic Impulse-Response Functions and Forecast Error Variance Decompositions
#' @export
#' @param AR_companion Autoregressive component of the companion matrix.
#' @param MA_companion Moving-average component of the companion matrix.
#' @param GARCH_estimation_results The object returned by the function call to 'Estimate_SVAR'. The estimation extracts the necessary parameters and series from this object.
#' @param Dates_Series Full set of dates the forecast error variance decompositions should be computed for.
#' @param variable_names Names of the series. Must follow the same ordering as the residuals.
#' @param linear_comb Optional, a named list of new variables that are linear combinations of the original series. The names should be entered as the list names, and the scalars to generate the linear combination should be the elements.
#' @param length_IRF For dynamic forecast error variance decompositions only. Maximum horizon to calculate for the impulse responses.
#' @param horizon_plot A list of numeric vectors. Forecast horizons for the forecast error variance decompositions to Show side-by-side at each date. Each list element will generate a panel grid of plots by horizon and variable.
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @references
#' Killian, L., and Lütkepohl, H. (2017) *Structural Vector Autoregressive Analysis*. Cambridge University Press.
#' ISBN: 978-1108164818
#'
#' Lütkepohl, H. (2006) *New Introduction to Multiple Time Series Analysis*.
#' ISBN: 978-3540262398
#' @return A (variable x shock x date x horizon) 4-dimensional array of forecast error variance and decompositions. Also returns dataframe versions and generates plots of each.
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
#' # Show Forecast Error Variance Decompositions at every period, at horizon lengths of 5, 20, 100
#' Dynamic_IRFS_FEVDs <- Dynamic_IRF_FEVD(
#'   AR_companion             = companion$AR_companion,
#'   GARCH_estimation_results = Estimation_results,
#'   Dates_Series             = Returnsdata[,1],
#'   variable_names           = names(Returnsdata[,2:4]),
#'   linear_comb              = list("BRK-SPX"=c(-1, 1, 0)),
#'   length_IRF               = 250,
#'   horizon_plot             = c(5, 20, 100)
#' )
#'
#'
Dynamic_IRF_FEVD <- function(
    AR_companion,
    MA_companion = NULL,
    GARCH_estimation_results,
    Dates_Series,
    variable_names,
    linear_comb = NULL,
    length_IRF,
    horizon_plot = 1
    ){


  # Initialize the Objects
  #-----------------------
  Structural_residuals <- GARCH_estimation_results$Struct_resid
  rotation_matrix      <- GARCH_estimation_results$rotation_matrix
  variance_series      <- GARCH_estimation_results$variance_series
  nshocks              <- ncol(Structural_residuals)
  nperiods             <- dim(Structural_residuals)[1]
  totvars <- nshocks + length(linear_comb)
  allvarnames <- c(variable_names, names(linear_comb))

  # Defining the dates
  Dates_residuals <- Dates_Series[(length(Dates_Series) - nperiods + 1):length(Dates_Series)]

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
  #- - - - - - - - - - - - - - - - - - - - - - - - -
  # Construct the total variance covariance matrix
  extended_condi_variance <- tcrossprod(Total_rotation_matrix)

  Omega     <- GARCH_estimation_results$omega
  Alpha     <- GARCH_estimation_results$alpha
  Beta      <- GARCH_estimation_results$beta

  # Compute the variance forecasts
  #-------------------------------
  Variance_forecasts      <- array(
    NA,
    dim = c(nperiods, nshocks, length_IRF),
    dimnames = list(
      as.character(Dates_residuals),
      variable_names,
      paste("horizon", 1:length_IRF)
    ))

  # Leading to get the right timing for the variance forecasts
  Variance_forecasts[,,1] <- dplyr::lead(variance_series, 1)

  AR_times_rotation       <- array(
    NA, dim = c(dim(Total_AR)[1], dim(Total_rotation_matrix)[2], length_IRF))
  AR_times_rotation[,,1]  <- Total_rotation_matrix


  # Loop over the horizon
  for(i in 2:length_IRF){
    Variance_forecasts[,,i] <- matrix(
      Omega * (1 - (Alpha + Beta)^(i-1)) / (1 - Alpha - Beta),
      nrow = nperiods, ncol = length(Omega), byrow = T) +
      t(diag((Alpha + Beta)^(i-1)) %*% t(Variance_forecasts[,,1]))

    AR_times_rotation[,,i]  <- Total_AR %*% AR_times_rotation[,,i-1]
  }


  # Setting up a multiplying matrix to generate the linear combinations
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


  # Coercing the forecasting matrix and computing the large Vectorized matrix
  # - This selects the interesting components for variance computation
  Constrained_AR_times_rotation <- array(
    NA, dim = c(totvars, dim(AR_times_rotation)[2], length_IRF))

  # - This stores the vectorized matrix
  Large_vec_matrix <- array(
    NA, dim = c(totvars, nshocks, length_IRF),
    dimnames = list(
      allvarnames,
      paste("shock", 1:nshocks),
      paste("horizon", 1:length_IRF)
    ))
  select_diag_mat  <- (
    diag(1, totvars) %x% diag(1, totvars)
  )[seq(from = 1, to = totvars^2, by = totvars+1),]
  ones             <- matrix(1, nrow = totvars, ncol = 1)

  # - This loops to construct the objects recursively
  for(i in 1:length_IRF){
    Constrained_AR_times_rotation[,,i] <- linear_comb_ext %*% AR_times_rotation[,,i]

    Large_vec_matrix[,,i] <- select_diag_mat %*% ((Constrained_AR_times_rotation[,,i] %x% ones) *
                                                    (ones %x% Constrained_AR_times_rotation[,,i]))
  }

  # Building the time series
  #-------------------------
  Variance_total_time_series <- array(
    NA,
    dim = c(nperiods, nshocks, totvars, length_IRF),
    dimnames = list(
      as.character(Dates_residuals),
      paste("shock", 1:nshocks),
      allvarnames,
      paste("horizon", 1:length_IRF)
    )
  )

  pb = txtProgressBar(min = 0, max = length_IRF, initial = 0)

  for(horizon in 1:length_IRF){

    # Print Progress
    setTxtProgressBar(pb,horizon)

    # Looping over shocks
    for(shock in 1:nshocks){
      Term_sum <- matrix(0, nrow = nperiods, ncol = totvars)

      for(i in 1:horizon){
        Term_sum <- Term_sum +
          matrix(Variance_forecasts[, shock, horizon - i + 1], ncol = 1) %*%
          matrix(Large_vec_matrix[, shock, i], nrow = 1)
      }

      Variance_total_time_series[, shock, , horizon] <- Term_sum

    }
  }
  close(pb)


  # Normalizing:
  #-------------
  Variance_decomp_time_series <- apply(Variance_total_time_series,
                                       c(1,3,4), function(x){x/sum(x)})

  # Transforming into data frames
  #------------------------------
  Variance_decomp_time_series_df <-
    tibble::as_tibble(as.data.frame.table(Variance_decomp_time_series))

  Variance_decomp_time_series_df <- Variance_decomp_time_series_df %>%
    rename(shock    = !!quote(Var1),
           date     = !!quote(Var2),
           variable = !!quote(Var3),
           horizon  = !!quote(Var4),
           value    = !!quote(Freq))
  Variance_decomp_time_series_df <- Variance_decomp_time_series_df %>%
    mutate(horizon = as.numeric(horizon),
           date = as.Date(lubridate::fast_strptime(as.character(date), "%Y-%m-%d")))
  change_variable_names          <- Variance_decomp_time_series_df$variable
  levels(change_variable_names)  <- allvarnames
  Variance_decomp_time_series_df$variable <- change_variable_names

  # Plotting
  #---------
  plots <- lapply(horizon_plot, Dynm_Vdecomp_Plot, df=Variance_decomp_time_series_df)


  return(list("detailed_variance_timeSeries"      = Variance_total_time_series,
              "Variance_decomposed_timeSeries"    = Variance_decomp_time_series,
              "Variance_decomposed_timeSeries_df" = Variance_decomp_time_series_df,
              "Variance_Decomposition_plots"      = plots
  ))
}




#' @rdname Dynamic_IRF_FEVD
#' @name Dynamic_FEVD_Plot
#' @export
#' @param horizon_plot A numeric vector of forecast horizons.
#' @param df A dataframe of the dynamic forecast error variance decompositions. Recommended to use the 'Variance_decomposed_timeSeries_df' component of the object returned by 'Get_Dynamic_IRF_FEVD'.
#' @importFrom ggplot2 ggplot aes theme theme_bw geom_hline ggtitle geom_col element_text labs facet_wrap vars geom_line scale_color_manual scale_linetype_manual element_blank element_text facet_grid geom_area scale_fill_viridis_d labeller label_both
#' @return A plot of dynamic forecast error variance decompositions at the given horizons.
Dynm_Vdecomp_Plot <- function(horizon_plot, df) {
  p <- ggplot2::ggplot(df %>% dplyr::filter(!!quote(horizon) %in% horizon_plot),
                       aes(x = !!quote(date), y = !!quote(value), fill = !!quote(shock))) +
    geom_area(alpha = 1) +
    scale_fill_viridis_d() + theme_bw() + labs(y="", x="") +
    facet_grid(horizon ~ variable, labeller = labeller(.rows = label_both))
  return(p)
}
