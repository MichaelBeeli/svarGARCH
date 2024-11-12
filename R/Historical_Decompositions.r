#' Build_Historical_Decomp
#' @rdname Historical_Decompositions
#' @name Build_Historical_Decomp
#' @title Build and Summarize Historical Decompositions
#' @export
#' @description Computes the historical decomposition of the time series from the structural shocks identified by GARCH.
#' @param AR_companion Autoregressive component of the companion matrix
#' @param MA_companion Moving average component of the companion matrix.
#' @param Structural_residuals Structural residuals identified by SVAR
#' @param rotation_matrix Rotation matrix identified by SVAR
#' @param Dates_Series Dates for the historical decomposition
#' @param variable_names Names of the series. Must follow the same ordering as the residuals.
#' @param linear_comb Optional, a named list of new variables that are linear combinations of the original series. The names should be entred as the list names, and the scalars to generate the linear combination should be the elements.
#' @importFrom ggplot2 ggplot aes theme theme_bw geom_hline ggtitle geom_col element_text labs
#' @importFrom dplyr select mutate group_by summarise
#' @importFrom magrittr %>%
#' @importFrom tidyselect any_of
#' @references
#' Killian, L., and Lütkepohl, H. (2017) *Structural Vector Autoregressive Analysis*. Cambridge University Press.
#' ISBN: 978-1108164818
#'
#' Lütkepohl, H. (2006) *New Introduction to Multiple Time Series Analysis*.
#' ISBN: 978-3540262398
#' @return Historical decompositions of Time Series by Structural Shocks. Returns an array and dataframe version. Also returns the provided specification for the linear combinations
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
#' # Historical Decompositions with a linear combination of the original series
#' Historical_Decompositions <- Build_Historical_Decomp(
#'   AR_companion = companion$AR_companion,
#'   MA_companion = companion$MA_companion,
#'   Structural_residuals = Estimation_results$Struct_resid,
#'   rotation_matrix      = Estimation_results$rotation_matrix,
#'   Dates_Series         = Returnsdata$Date,
#'   variable_names       = names(Returnsdata[,2:4]),
#'   linear_comb          = list("BRK-SPX"=c(-1, 1, 0))
#' )
#'
#' # Summarizing Historical Decompositions
#' Hist_decomp_plots <- Summary_Historical_Decomp(
#'   historical_df  = Historical_Decompositions$Historical_DataFrame,
#'   true_df        = Returnsdata,
#'   variable_names = names(Returnsdata[,2:4]),
#'   linear_comb    = Historical_Decompositions$Linear_comb
#' )
#'

# HISTORICAL DECOMPOSITION ----
#=========================
Build_Historical_Decomp <- function(
    AR_companion,
    MA_companion = NULL,
    Structural_residuals,
    rotation_matrix,
    Dates_Series,
    variable_names,
    linear_comb = NULL
){

  # Initialize the final matrix object
  #-----------------------------------
  nshocks           <- ncol(Structural_residuals)
  decomposed_series <- matrix(
    0,
    nrow = nrow(Structural_residuals) + 1,
    ncol = ncol(Structural_residuals)^2
  ) # We implicitly assume as many shocks as variables

  # Initialize useful objects for the computation
  #----------------------------------------------

  # Pre-multiplying matrix extended
  rotation_matrix_extended <- NULL
  for(i in 1:nshocks){
    rotation_matrix_extended <- rbind(
      rotation_matrix_extended,
      diag(rotation_matrix[i,])
    )
  }

  # Calculating the AR and MA orders
  AR_order <- ncol(AR_companion)/nshocks
  MA_order <- length(MA_companion)

  # Coercing AR Companion matrix to first n rows
  AR_first_rows <- AR_companion[1:nshocks,]

  # Defining the dates
  Dates_residuals <- Dates_Series[(max(AR_order, MA_order)+1):length(Dates_Series)]

  # Performing loop to reconstruct the historical series
  #-----------------------------------------------------
  if(MA_order == 0){
    for(i in (AR_order + 1):nrow(decomposed_series)){
      # We first focus on the autoregressive term
      # Initialize the term
      Term_AR <- rep(0, ncol(decomposed_series))

      # We add successively the different lag orders
      for(j in 1:AR_order){
        Term_AR <- Term_AR +
          c(
            (AR_first_rows[, (nshocks * (j-1) + 1):(nshocks * j)] %x% diag(1, nshocks)) %*%
              decomposed_series[i-j,]
          )
      }


      # Adding the shocks
      decomposed_series[i,] <- c(
        Term_AR + rotation_matrix_extended %*% Structural_residuals[i-1,]
      )
    }

    # BELOW IS THE MOVING AVERAGE CASE
    #- - - - - - - - - - - - - - - - -
  } else {
    # Structural shocks to use in loop in case MA terms
    last_structural_shock <- matrix(0, nrow = nshocks, ncol = 1)

    # MA matrix to use in loop
    matrix_MA_extended <- NULL
    for(j in 1:MA_order){
      temp_column_matrix <- NULL

      for(i in 1:nshocks){
        temp_column_matrix <- rbind(
          temp_column_matrix,
          diag((MA_companion[[j]] %*% rotation_matrix)[i,])
        )
      }

      matrix_MA_extended <- cbind(
        matrix_MA_extended,
        temp_column_matrix
      )
      # Resulting matrix is (n^2 x nq) collapsing together the different
      # MA terms from 1 to q.
    }

    # Loop
    #-----
    for(i in (max(AR_order, MA_order) + 1):nrow(decomposed_series)){
      # Initialize the first terms
      Term_AR <- rep(0, ncol(decomposed_series))

      # We add successively the different lag orders
      for(j in 1:AR_order){
        Term_AR <- Term_AR +
          c(
            (AR_first_rows[, (nshocks * (j-1) + 1):(nshocks * j)] %x% diag(1, nshocks)) %*%
              decomposed_series[i-j,]
          )
      }

      # MA term
      residuals_for_MA <- NULL
      if(i == 2){
        residuals_for_MA <- last_structural_shock
      } else if(i == MA_order + 1){
        residuals_for_MA <- NULL
        for(j in 1:(MA_order - 1)){
          residuals_for_MA <- c(
            residuals_for_MA,
            Structural_residuals[i - 1 - j,])
        }
        residuals_for_MA <- c(residuals_for_MA, last_structural_shock)
      } else {
        residuals_for_MA <- NULL
        for(j in 1:MA_order){
          residuals_for_MA <- c(
            residuals_for_MA,
            Structural_residuals[i - 1 - j,])
        }
      }
      Term_MA <- matrix_MA_extended %*% residuals_for_MA

      # Adding the shocks
      decomposed_series[i,] <- c(
        Term_AR + Term_MA +
          rotation_matrix_extended %*% Structural_residuals[i-1,]
      )
    }
  }

  # Get linear combinations
  #=========
  if(!is.null(linear_comb)) {

      lcombs <- matrix(0, nrow = nrow(Structural_residuals) + 1, ncol = length(linear_comb)*nshocks)
      array_decomp <- array(decomposed_series, dim = list(nrow(Structural_residuals) + 1, nshocks,nshocks))

      for (i in 1:length(linear_comb)) {

        start <- (i-1)*nshocks+1
        stop <- i*nshocks

        lcombs[, start:stop] <- t(apply(array_decomp, MARGIN = 1, function(x) x %*% linear_comb[[i]]))

      }
      decomposed_series <- cbind(decomposed_series, lcombs)



  }
  totvars <- nshocks + length(linear_comb)
  allvarnames <- c(variable_names, names(linear_comb))
  # Reconstruct a data frame
  #=========
  var_name <- NULL
  shock_name <- NULL

  for (i in 1:totvars) {

    var_name <- c(
      var_name,
      rep(allvarnames[i], nshocks * nrow(Structural_residuals))
    )

    for (j in 1:nshocks) {
      shock_name <- c(
        shock_name,
        rep(paste("shock",j), nrow(Structural_residuals))
      )
    }

  }
  # Form the data frame
  #=========
  decomposed_series_DataFrame <- data.frame(
    "Date"     = Dates_residuals,
    "variable" = factor(var_name, levels = unique(var_name)),
    "shock"    = shock_name,
    "value"    = c(decomposed_series[2:nrow(decomposed_series),])
  )

  # Form the full array
  #=========
  Historical_array <- array(decomposed_series, dim = list(nrow(Structural_residuals) + 1, nshocks,totvars),
                            dimnames = list(
                              paste("date", 1:(nrow(Structural_residuals) + 1)),
                              paste("shock", 1:nshocks),
                              allvarnames
                            ))

  return(
    list(
      "Historical_array"     = Historical_array,
      "Historical_DataFrame" = decomposed_series_DataFrame,
      "Linear_comb"          = linear_comb
    )
  )
}

#' Summary_Historical_Decomp
#' @rdname Historical_Decompositions
#' @name Summary_Historical_Decomp
#' @export
#' @description Convenience function to visualize the historical decompositions and the extent of estimation error in the historical decomposition
#' @param historical_df The "Historical_Dataframe" object returned by the "Build_Historical_Decomp" function.
#' @param true_df A dataframe of the true series. In order, should include a column named "Date", followed by each of the columns of the original series in order, excluding linear combinations.
#' @param linear_comb Optional, a named list of new variables that are linear combinations of the original series. The names should be entred as the list names, and the scalars to generate the linear combination should be the elements.
#' @param demeaned Logical value for whether the true_df has been demeaned. If false, it will be demeaned before the comparison is shown. Default is true.
#' @importFrom ggplot2 ggplot aes theme theme_bw geom_hline ggtitle geom_col element_text labs facet_wrap vars geom_line scale_color_manual scale_linetype_manual element_blank element_text facet_grid
#' @importFrom dplyr select mutate group_by summarise rename full_join pull across
#' @importFrom magrittr %>%  %<>%
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom ggplot2 .data
#' @return Plots of Historical Decomposition Time Series
Summary_Historical_Decomp <- function(
    historical_df,
    true_df,
    variable_names,
    linear_comb=NULL,
    demeaned=T
) {

  #### Ensuring comparison can be performed properly
  truecols <- ncol(true_df)
  names <- historical_df %>% dplyr::pull("variable")

  hist_base_cols <- setdiff(unique(names), names(linear_comb))

  if (names(true_df)[1]!="Date") {
    stop("Please ensure that the column 'Date' is the first column of the true dataframe.")
  }
  if(!(truecols-1 == length(hist_base_cols))) {
    stop("Historical and True Dataframes are incompatible. Ensure that column names and order are identical")
  }
  if(length(setdiff(names(true_df[,-1]), hist_base_cols)) + length(setdiff(hist_base_cols, names(true_df[,-1])) > 0)) {
     stop("Historical and True Dataframes are incompatible. Ensure that column names and order are identical")
  }

  #### Demeaning true dataframe if necessary
  if(demeaned==F) {
    true_df[,2:ncol(true_df)] <- sapply(true_df[,2:ncol(true_df)], function(x) x - mean(x))
  }


  #### Get linear combinations of True Data
  for (i in seq_len(length(linear_comb))) {

    name <- names(linear_comb)[i]

    true_df[name] <- as.matrix(true_df %>% select(tidyselect::any_of(variable_names))) %*% linear_comb[[i]]

  }

  #### Prepare clean data for merge
  true_df_clean <- true_df   %>% pivot_longer(
    cols = -c(!!quote(Date)),
    names_to = "variable",
    values_to = "True"
  )

  #### Merge the two dataframes together to allow for comparisons
  estdf <- historical_df %>% group_by(across(tidyselect::all_of(c("Date", "variable")))) %>%
    summarise(Estimated=sum(!!quote(value), na.rm=T))

  compdf <- full_join(estdf, true_df_clean, by=c("Date", "variable")) %>%
    pivot_longer(
      cols = c("Estimated", "True"),
      names_to = "Source",
      values_to = "value"
    ) %>% mutate(variable=toupper(!!quote(variable)))

  diffdf <- compdf %>% pivot_wider(
    id_cols = c("Date", "variable"),
    names_from = "Source",
    values_from = "value"
  ) %>% mutate(Error=!!quote(True)-!!quote(Estimated))

  # PLOTTING
  #=========


  p <- ggplot(historical_df, aes(x = .data$Date, y = .data$value)) +
    geom_line() +
    facet_grid(.data$variable~.data$shock, scales = "free")  + theme_bw() + theme(
      legend.position = "bottom",
      strip.background = element_blank()
    ) + labs(title="Historical Decompositions")


  cplot <- ggplot(compdf, aes(x= .data$Date, y= .data$value, linetype= .data$Source, color= .data$Source)) +
    facet_wrap(facets=vars(.data$variable), scales="free") +
    geom_line(alpha=0.7) + scale_linetype_manual(values = c(2,1)) + scale_color_manual(values=c("red", "black")) + theme_bw() + theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_text(size=12)
    ) + labs(title="Historical Reconstruction Comparison", y="")

  eplot <- ggplot(diffdf, aes(x= .data$Date, y= .data$Error)) +
    facet_wrap(facets=vars( .data$variable), scales="free") +
    geom_line(alpha=0.7, color="red")  + theme_bw() + theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size=12)
    ) + labs(title="Historical Reconstruction Error", y="") +
    geom_hline(yintercept = 0, linetype=2)


  # SHOW PRIMARY PLOT
  #=========
  suppressWarnings(print(p))


  # Send plots to user
  #---------------------
  return(
    list(
      "Historical_Decompositions"            = p,
      "Historical_Reconstruction_Comparison" = cplot,
      "Historical_Reconstruction_Error"      = eplot
    )
  )
}
