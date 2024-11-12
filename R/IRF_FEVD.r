#' IRF_FEVD
#' @rdname IRF_FEVD
#' @description Computes Impulse Responses and Forecast Error Variance Decompositions up until specified horizon.
#' @title Impulse-Response Functions and Forecast Error Variance Decompositions
#' @export
#' @param AR_companion Autoregressive component of the companion matrix.
#' @param MA_companion Moving-average component of the companion matrix.
#' @param GARCH_estimation_results The object returned by the function call to 'Estimate_SVAR'. The structural residuals and rotation matrix are extracted and used in the estimation.
#' @param variable_names Names of the series. Must follow the same ordering as the residuals.
#' @param linear_comb Optional, a named list of new variables that are linear combinations of the original series. The names should be entered as the list names, and the scalars to generate the linear combination should be the elements.
#' @param horizon Horizon that the standard impulse responses and forecast error variance decompositions should be calculated for.
#' @importFrom viridis viridis
#' @importFrom ggplot2 ggplot aes theme theme_bw geom_hline ggtitle geom_col element_text labs facet_wrap vars geom_line scale_color_manual scale_linetype_manual element_blank element_text facet_grid scale_color_viridis_d scale_alpha_manual
#' @references
#' Killian, L., and Lütkepohl, H. (2017) *Structural Vector Autoregressive Analysis*. Cambridge University Press.
#' ISBN: 978-1108164818
#'
#' Lütkepohl, H. (2006) *New Introduction to Multiple Time Series Analysis*.
#' ISBN: 978-3540262398
#'  @return A (variable x shock x horizon) 3-dimensional array of impulse responses and forecast error variance decompositions. Also returns dataframe versions and generates plots of each.
#' @examples
#' Returnsdata <- svarGARCH::Returns
#' companion <- companion_matrix(Returnsdata[,2:4])
#'
#' # Default estimation method, no fixed parameters, no initial guess for parameters
#' Estimation_results <- Estimate_SVAR(
#'    Residuals                = companion$Residuals,
#'    variable_names           = names(Returnsdata[,2:4]),
#'    fixvariance              = TRUE,
#'    variance_targeting = TRUE
#' )
#'
#' IRFS_FEVDs <- IRF_FEVD(
#'   AR_companion             = companion$AR_companion,
#'   MA_companion             = companion$MA_companion,
#'   GARCH_estimation_results = Estimation_results,
#'   variable_names           = names(Returnsdata[,2:4]),
#'   linear_comb              = list("BRK-SPX"=c(-1, 1, 0)),
#'   horizon                  = 250
#' )
#'
#'
IRF_FEVD <- function(
    AR_companion,
    MA_companion = NULL,
    GARCH_estimation_results,
    variable_names,
    linear_comb = NULL,
    horizon = 251
    ){

  # Initialize the Objects
  #-----------------------
  Structural_residuals <- GARCH_estimation_results$Struct_resid
  rotation_matrix      <- GARCH_estimation_results$rotation_matrix
  nshocks       <- ncol(Structural_residuals)
  ncombs        <- dplyr::coalesce(length(linear_comb), 0)
  totvars       <- ncol(rotation_matrix) + ncombs
  allvarnames   <- c(variable_names, names(linear_comb))

  # Calculating the AR and MA orders
  AR_order  <- ncol(AR_companion)/nshocks
  MA_order  <- length(MA_companion)
  TotOrder <- AR_order + MA_order

  # Taking size of the extended system
  nextended <- (AR_order + MA_order) * totvars

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


  # Initialize the objects to retrieve
  IRF         <- array(NA, c(totvars, nshocks, horizon),
                       dimnames = list(
                         allvarnames,
                         paste("shock", 1:nshocks),
                         paste("horizon", 1:horizon)
                       ))

  Variance    <- array(NA, c(nshocks*(TotOrder), nshocks*(TotOrder), horizon),
                       dimnames = list(
                         rep(paste("shock", 1:nshocks), TotOrder),
                         rep(paste("shock", 1:nshocks), TotOrder),
                         paste("horizon", 1:horizon)
                       ))
  Vdecomp     <- array(NA, c(totvars, nshocks, horizon),
                       dimnames = list(
                         allvarnames,
                         paste("shock", 1:nshocks),
                         paste("horizon", 1:horizon)
                       ))

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

  # Compute the quantities of interest
  #-----------------------------------
  # IRFS
  Term_AR <- diag(1, nrow(Total_AR))
  for(i in 1:horizon){
    IRF[,,i] <- (linear_comb_ext %*% (Term_AR) %*% Total_rotation_matrix)
    Term_AR <- Term_AR %*% Total_AR
  }

  # Total Variance of Forecast Errors
  Variance[,,1] <- extended_condi_variance
  for(i in 2:horizon){
    Variance[,,i] <- Total_AR %*% Variance[,,i-1] %*% t(Total_AR) +
      extended_condi_variance
  }

  Variance_linear_comb  <-
    array(
      apply(Variance, 3, function(x){linear_comb_ext %*% x %*% t(linear_comb_ext)}),
      c(nrow(linear_comb_ext), nrow(linear_comb_ext), dim(Variance)[3])
    )

  E   <- diag(1, totvars)
  E2  <- diag(1, nshocks)

  Vdecomp[,,1] <-  (linear_comb_ext %*%  Total_rotation_matrix)^2
  Term_AR <- Total_AR
  for(k in 2:horizon){
    for(i in 1:nrow(linear_comb_ext)){
      for(j in 1:ncol(Total_rotation_matrix)){
        Vdecomp[i,j,k] <-
          (linear_comb_ext[i,] %*% Term_AR %*%
             Total_rotation_matrix %*% E2[,j])^2 +
          Vdecomp[i,j,k-1]
      }
    }
    Term_AR <- Term_AR %*% Total_AR
  }


  # Normalizing by total variance
  for(k in 1:horizon){
    Vdecomp[,,k] <- (diag(1/diag(Variance_linear_comb[,,k]))) %*% Vdecomp[,,k]
  }


  # Generate final arrays

   ## Create data frame of IRFs for linear combinations of original variables
    IRF_data_frame <- data.frame(
      "horizon" = c(0:(horizon-1) %x% matrix(1,nshocks*nrow(linear_comb_ext),1)),
      "variable" = factor(rownames(linear_comb_ext), levels = rownames(linear_comb_ext)),
      "shock" = paste("shock", c(c(1:nshocks) %x% matrix(1,nrow(linear_comb_ext),1))),
      "value" = c(IRF[1:nrow(linear_comb_ext),,])
    )

  # Generate Variance Decomposition Dataframe
    Vdecomp_data_frame <- data.frame(
      "horizon"  = c(0:(horizon-1) %x% matrix(1,nshocks*nrow(linear_comb_ext),1)),
      "variable" = rownames(Vdecomp),
      "shock"    = paste("shock", c(c(1:nshocks) %x% matrix(1,nrow(linear_comb_ext),1))),
      "value"    = c(Vdecomp[1:totvars,,])
    ) %>% dplyr::group_by(horizon, !!quote(variable)) %>%
      mutate(VarShare = cumsum(!!quote(value)))


  # Plotting
  #---------
    colors_used <- viridis(nshocks)

    nplot <- nrow(linear_comb_ext)

    plot_vdecomp <- ggplot(Vdecomp_data_frame, aes(x=!!quote(horizon), y=!!quote(VarShare), color=!!quote(shock), alpha=!!quote(shock))) +
      geom_line() +
      facet_wrap(facets = vars(!!quote(variable)), nrow = 1) + theme_bw() + scale_color_viridis_d() + scale_alpha_manual(values = c(rep(1, nshocks-1), 0)) +labs( y="")

    p <- ggplot(IRF_data_frame, aes(x = !!quote(horizon), y = !!quote(value)))+
      geom_line(color="royalblue", lty=1, linewidth=1.2, alpha=0.4) + facet_grid(variable~shock, scales  ="free") + theme_bw() + geom_hline(yintercept = 0, linetype=3) +labs(y="")


  # RETURN RESULTS
    return(list(
      "IRF" = IRF,
      "IRF_DataFrame" = IRF_data_frame,
      "Vdecomp" = Vdecomp,
      "Vdecomp_DataFrame" = Vdecomp_data_frame,
      "IRF_plot" = p,
      "Vdecomp_plot" = plot_vdecomp
    ))


}

