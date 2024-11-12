#' summary_heteroskedasticity
#' @rdname summary_heteroskedasticity
#' @export
#' @name summary_heteroskedasticity
#' @title Summary of Heteroskedasticity present in the data
#' @description This function provides a simple analysis of heteroskedasticity of errors. It can be used to quickly check for positive evidence of heteroskedasticity. Generates ACFs and PACFs of squared errors, and the results of the Multivariate ARCH test from the MTS package.
#' @param Residuals Residuals from a time series system
#' @param lagsize The lag length for the ACFs
#' @importFrom ggplot2 ggplot aes theme theme_bw geom_hline ggtitle geom_col element_text labs
#' @importFrom ggpubr annotate_figure text_grob ggarrange
#' @importFrom purrr map2
#' @importFrom stats acf pacf
#' @importFrom MTS MarchTest
#' @seealso  [MTS::MarchTest] For suite of tests to detect conditional heteroskedasticity in Multivariate Time Series. 
#' @references 
#' Lütkepohl, H. (2006) *New Introduction to Multiple Time Series Analysis*. 
#' ISBN: 978-3540262398
#' @return Multivariate ARCH Test Results, Auto-Correlation Functions of Squares
#' @examples
#'
#' Returnsdata <- svarGARCH::Returns 
#' companion <- companion_matrix(Returnsdata[,2:4], p=1, q=1)
#' summary_heteroskedasticity(companion$Residuals, lagsize=12)
#'
#'
#' Series <- matrix(cbind(rnorm(20000, mean=1, sd=1)), nrow=10000, ncol=2)
#' Model <- MTS::VAR(Series, p=1)
#' summary_heteroskedasticity(Model$residuals, lagsize=12)

summary_heteroskedasticity <- function(Residuals, lagsize=12) {



  lagplus <- lagsize+1
  # Function to generate pacfs --
  genpacfs <- function(x, y, partial=FALSE) {
    series <- stats::acf(x^2, lag.max = lagsize, plot = FALSE)$acf[2:lagplus]

    if (partial == TRUE) {
      series <- stats::pacf(x^2, lag.max = lagsize, plot = FALSE)$acf
    }

    sig <- 1.96/sqrt((length(x)))

    df <- data.frame(cbind(seq(1,lagsize), series))

    plot <- ggplot(data=df, aes(x=factor(!!quote(V1)), y= series)) +
      theme_bw() +
      geom_col(width = 0.05, color="black") +
      geom_hline(yintercept = sig, linetype="dashed", color="blue") + geom_hline(yintercept = -1*sig, linetype="dashed", color="blue")  + labs(x="", y="") + ggtitle(paste0("Series", y)) + theme(
        plot.title = element_text(size=16, hjust=0.5),
        axis.text = element_text(color="black", size=12, angle=90)
      )

    return(plot)
  }



  res <- data.frame(Residuals)

  acfs <- map2(res, seq(1, ncol(Residuals)), genpacfs, partial=FALSE)
  pacfs <- map2(res, seq(1, ncol(Residuals)), genpacfs, partial=TRUE)

  plot <- ggarrange(plotlist = acfs)
  plot_p <- ggarrange(plotlist = pacfs)


  acfplot <- annotate_figure(plot, top = text_grob(paste0("ACF Squared Errors"), size = 22, face="bold"))
  pacfplot <- annotate_figure(plot_p, top = text_grob(paste0("PACF Squared Errors"), size = 22, face="bold"))


  #### Next need to implement ARCH tests --> see Lutkepohl & MTS package

  MTS::MarchTest(Residuals, lag = lagsize)

  readline(prompt="Press [enter] to see ACFS")

  print(acfplot)

  readline(prompt="Press [enter] to see PACFS")

  print(pacfplot)

  return(list("ACFs"=acfs, "PACFs"=pacfs))

}

