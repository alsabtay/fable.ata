#' @docType package
#' @keywords package
"_PACKAGE"

#' @importFrom ATAforecasting ATA.SeasAttr ATA.BoxCoxAttr ATA
#' @importFrom tsibble measured_vars as_tibble tibble
#' @importFrom stats frequency ts
train_ata <- function(.data, specials, ...){
  if(length(tsibble::measured_vars(.data)) > 1){
    stop("Only univariate responses are supported by ATA Method")
  }
  level <- specials$level[[1]]
  trend <- specials$trend[[1]]
  season <- specials$season[[1]]
  accuracy <- specials$accuracy[[1]]
  transform <- specials$transform[[1]]
  holdout <- specials$holdout[[1]]
  # Prepare data for modelling
    model_data <- tsibble::as_tibble(.data)[c(rlang::expr_text(index(.data)), tsibble::measured_vars(.data))]
    colnames(model_data) <- c("ds", "y")
    if (any(is.na(model_data$y))) {
      stop("ATA method does not support missing values.")
    }
    period <- fabletools::get_frequencies(season$period, .data, .auto = "smallest")
    pre_data <- quietly(tsbox::ts_ts)(.data)
    model_data <- stats::ts(model_data$y, frequency = unname(period), start=start(pre_data))
  if (holdout$holdout == TRUE & accuracy$criteria == "AMSE") {
    stop("ATA Method does not support 'AMSE' for 'holdout' forecasting.")
  }
  seas_opt_crit <- ATAforecasting::ATA.SeasAttr(suroot.test=season$suroot_test, corrgram.tcrit = season$suroot_tcrit, suroot.uroot = season$suroot_uroot, uroot.pkg="tseries",
                                uroot.test=trend$uroot_test, uroot.type=trend$uroot_type, uroot.alpha=trend$uroot_alpha, uroot.maxd=1)
  bc_opt_crit <- ATAforecasting::ATA.BoxCoxAttr(bcMethod = transform$bcMethod, bcLower = transform$bcLower, bcUpper = transform$bcUpper)

  # Build and train model
  mdl_ATA <- safely(quietly(ATAforecasting::ATA))(
                  X = model_data,
                  Y = NULL,
                  parP = level$parP,
                  parQ = trend$parQ,
                  parPHI = trend$parPHI,
                  model.type = ifelse(trend$type=="N","A",ifelse(trend$type=="Ad", "A", ifelse(trend$type=="Md","M",trend$type))),
                  seasonal.test = ifelse(season$type=="N", FALSE, season$test),
                  seasonal.model = ifelse(season$type=="N", "none", season$method),
                  seasonal.period = unname(period),
                  seasonal.type = ifelse(season$type=="N", "M", season$type),
                  seasonal.test.attr = seas_opt_crit,
                  find.period = NULL,
                  accuracy.type = accuracy$criteria,
                  nmse = accuracy$nmse,
                  level.fixed = level$level_fixed,
                  trend.fixed = ifelse(trend$opt_trend=="fixed", TRUE, FALSE),
                  trend.search = ifelse(trend$opt_trend=="search", TRUE, FALSE),
                  h = 1,
                  partition.h = ifelse(holdout$holdout == TRUE, holdout$set_size, NULL),
                  holdout = holdout$holdout,
                  holdout.adjustedP = holdout$adjustedP,
                  holdin = FALSE,
                  transform.order = ifelse(transform$order=="none", "before", transform$order),
                  transform.method = ifelse(transform$method=="none", NULL, transform$method),
                  transform.attr = bc_opt_crit,
                  lambda = transform$lambda,
                  shift = transform$shift,
                  initial.level = level$initial_level,
                  initial.trend = trend$initial_trend,
                  ci.level = 95,
                  start.phi = trend$parPHI_range[1],
                  end.phi = trend$parPHI_range[2],
                  size.phi = trend$parPHI_increment,
                  negative.forecast = TRUE,
                  print.out = FALSE,
                  plot.out = FALSE)

  # Return model
  ata_out <-  list(
        "par" = tsibble::tibble(term = names(mdl_ATA$par.specs), estimate = unname(mdl_ATA$par.specs)),
        "est" = list(".fitted" = mdl_ATA$fitted,
                   ".resid" = mdl_ATA$residuals),
        "fit" = mdl_ATA$accuracy$fits,
        "components" = list("coefp" = mdl_ATA$coefp,
                          "coefq" = mdl_ATA$coefp,
                          "P" = mdl_ATA$parP,
                          "Q" = mdl_ATA$parQ,
                          "PHI" = mdl_ATA$parPHI,
                          "response" = mdl_ATA$actual,
                          "level" = mdl_ATA$level,
                          "trend" = mdl_ATA$trend,
                          "seasonal" = mdl_ATA$seasonal,
                          "seasindex" = mdl_ATA$seasonal.index,
                          "seasadj" = mdl_ATA$seasonal.adjusted,
                          "remainder" = mdl_ATA$residuals),
        "transform" = list("method" = transform$method,
                         "order" = transform$order,
                         "lambda" = mdl_ATA$lambda,
                         "shift" = mdl_ATA$shift),
        "holdout" = list("holdout" = mdl_ATA$holdout,
                       "adjustment" = mdl_ATA$holdout.adjustedP),
        "spec" = list("errortype" = "A",
                 "trendtype" = mdl_ATA$model.type,
                 "seasontype" = mdl_ATA$seasonal.model,
                 "damped" = ifelse(mdl_ATA$parPHI==1, FALSE, TRUE),
                 "period" = mdl_ATA$seasonal.period,
                 "method" = mdl_ATA$method),
        "model_output" = mdl_ATA)
        structure(ata_out, class = "fable_ata")
}

specials_ata <- fabletools::new_specials(
   level = function(parP= NULL, level_fixed = FALSE, initial_level = FALSE){
                   list("parP" = parP, "level_fixed" = level_fixed, "initial_level" = initial_level)
                  },
   trend = function(type = "A", parQ = NULL, initial_trend = FALSE, opt_trend = "none",
                   parPHI = NULL, parPHI_range = c(0.8, 1.0), parPHI_increment = 0.01,
                   uroot_test = "adf", uroot_alpha = 0.05, uroot_type = "level")
                   {

                     if (type == "N"){
                       parQ = 0
                       warning("Q has been set 0 because of no trend option.")
                     }
                     if (type == "Ad" | type == "Md"){
                       parPHI = NULL
                       warning("PHI has been set NULL as damped trend is choosen.")
                     }
                     list("type" = type, "parQ" = parQ, "initial_trend" = initial_trend,
                          "parPHI" = parPHI, "parPHI_range" = parPHI_range, "parPHI_increment" = parPHI_increment,
                          "uroot_test" = uroot_test, "uroot_alpha" = uroot_alpha, "uroot_type" = uroot_type)
                  },
   season = function(type = "A", test = TRUE, period = NULL, method = "decomp",
                  suroot_test = "correlogram", suroot_tcrit = 1.28, suroot_uroot = TRUE)
                  {
                    if (type == "N") {
                      period <- 1
                    }
                    list("type" = type, "test" = test, "period" = period, "method" = method,
                          "suroot_test" = suroot_test, "suroot_tcrit" = suroot_tcrit, "suroot_uroot" = suroot_uroot)
                  },
   accuracy = function(criteria="sMAPE", nmse = 3, ic = "AIC"){
                     if (nmse > 30 & criteria == "AMSE") {
                       nmse <- 30
                       warning("'nmse' must be less than 30. 'nmse' is set to 30.")
                     }else if ((is.null(nmse) | nmse <= 1) & criteria == "AMSE") {
                       nmse <- 3
                       warning("'nmse' must be greater than 1. 'nmse' is set to 3.")
                     }else{
                     }
                      list("criteria" = criteria, "nmse" = nmse, "ic" = ic)
                    },
   transform = function(method="none", order = "none", lambda = NULL, shift = 0,
                      bcMethod = "guerrero", bcLower = 0, bcUpper = 5){
                        list("method" = method, "order" = order, "lambda" = lambda, "shift" = shift,
                        "bcMethod" = bcMethod, "bcLower" = bcLower, "bcUpper" = bcUpper)
                    },
   holdout = function(holdout = FALSE, adjustment = TRUE, set_size = NULL ){
                        list("holdout" = holdout, "adjustment" = adjustment, "set_size" = set_size)
                    }
  .required_specials = c("level", "trend", "season", "accuracy", "transform", "holdout")

)

#' ATAforecasting: Automatic Time Series Analysis and Forecasting using Ata Method with Box-Cox Power Transformations Family and Seasonal Decomposition Techniques
#'
#' Returns ATA(p,q,phi) applied to \code{X}.
#' The Ata method based on the modified simple exponential smoothing as described in Yapar, G. (2016) <doi:10.15672/HJMS.201614320580> ,
#' Yapar G., Capar, S., Selamlar, H. T., Yavuz, I. (2017) <doi:10.15672/HJMS.2017.493> and Yapar G., Selamlar, H. T., Capar, S., Yavuz, I. (2019)
#' <doi:10.15672/hujms.461032> is a new univariate time series forecasting method which provides innovative solutions to issues faced during
#' the initialization and optimization stages of existing methods.
#' Forecasting performance of the Ata method is superior to existing methods both in terms of easy implementation and accurate forecasting.
#' It can be applied to non-seasonal or seasonal time series which can be decomposed into four components (remainder, level, trend and seasonal).
#' This methodology performed well on the M3 and M4-competition data.
#' Returns ATA(p,q,phi) applied to \code{X}.
#'
#' @param formula Model specification (see "Specials" section).
#' @param ... Other arguments
#'
#' @section Specials:
#'
#' The _specials_ define the methods and parameters for the components (level, trend, seasonality, accuracy, transform, holdout) of an ATA method.
#'
#' There are a couple of limitations to note about ATA method:
#'
#' - It supports only additive error term.
#' - It does not support exogenous regressors.
#' - It does not support missing values. You can complete missing values in the data with imputed values (e.g. with [tsibble::fill_gaps()], [tidyr::fill()], or by fitting a different model type and then calling [fabletools::interpolate()]) before fitting the model.
#'
#' \subsection{level}{
#' The `level` special is used to specify the form of the level term.
#' \preformatted{
#' level(parP = NULL, level_fixed = TRUE, initial_level = FALSE)
#' }
#'
#' \tabular{ll}{
#'   `parP`     \tab The value of the smoothing parameter for the level. If `p = 0`, the level will not change over time. Conversely, if `p = 1` the level will update similarly to a random walk process. If NULL or "opt", it is estimated. \code{p} has all integer values from 1 to \code{length(data)}. \cr
#'   `level_fixed`      \tab If TRUE, "pStarQ"  --> First, fits ATA(p,0) where p = p* is optimized for q=0. Then, fits ATA(p*,q) where q is optimized for p = p*. \cr
#'   `initial_level`     \tab If NULL, FALSE is default. If FALSE, ATA Method calculates the pth observation in \code{data} for level. If TRUE, ATA Method calculates average of first p value in \code{data}for level. \cr
#'
#' }
#' }
#'
#' \subsection{trend}{
#' The `trend` special is used to specify the form of the trend term and associated parameters.
#' \preformatted{
#' trend(type = "A", parQ = NULL, initial_trend = FALSE, opt_trend = "none",
#'        parPHI = NULL, parPHI_range = c(0.8, 1.0), parPHI_increment = 0.01,
#'        uroot_test = "adf", uroot_alpha = 0.05, uroot_type = "level")
#' }
#'
#' \tabular{ll}{
#'   `type`     \tab The form of the trend term: either none ("N"), additive ("A"), multiplicative ("M") or damped variants ("Ad", "Md"). \cr
#'   `parQ`      \tab The value of the smoothing parameter for the slope. If `Q = 0`, the slope will not change over time. Conversely, if `Q = 1` the slope will have no memory of past slopes. \cr
#'   `parPHI` \tab The value of the dampening parameter for the slope. If `phi = 0`, the slope will be dampened immediately (no slope). Conversely, if `phi = 1` the slope will not be dampened. \cr
#'   `parPHI_range`       \tab If `phi=NULL`, `phi_range` provides bounds for the optimised value of `phi`.\cr
#'   `parPHI_increment`  \tab If `phi=NULL`, `parPHI_increment` provides increment step for searching `phi`. If NULL, `parPHI_increment` will be determined as the value that allows the `parPHI_range` to be divided into 20 equal parts. \cr
#'   `initial_trend`        \tab If NULL, FALSE is default. If FALSE, ATA Method calculates the qth observation in \code{X(T)-X(T-1)} for trend. If TRUE, ATA Method calculates average of first q value in \code{X(T)-X(T-1)} for trend. \cr
#'   `opt_trend`        \tab The value of the dampening parameter for the slope. If `phi = 0`, the slope will be dampened immediately (no slope). Conversely, if `phi = 1` the slope will not be dampened. \cr
#'   `uroot_test`        \tab Type of unit root test before all type seasonality test. Possible values are "adf", "pp" and "kpss". \cr
#'   `uroot_alpha`   \tab Significant level of the unit root test, possible values range from 0.01 to 0.1. \cr
#'   `uroot_type`        \tab Specification of the deterministic component in the regression for unit root test. Possible values are "level" and "trend". \cr
#' }
#' }
#'
#' \subsection{season}{
#' The `season` special is used to specify the form of the seasonal term and associated parameters. To specify a nonseasonal model you would include `season(method = "N")`.
#' \preformatted{
#' season(type = "A", test = TRUE, period = NULL, method = "decomp",
#'        suroot_test = "correlogram", suroot_tcrit = 1.28, suroot_uroot = TRUE)
#' }
#' \tabular{ll}{
#'   `type`     \tab The form of the seasonal term: either none ("N"), additive ("A") or multiplicative ("M"). \cr
#'   `test`     \tab Testing for stationary and seasonality. If TRUE, the method firstly uses \code{test="adf"}, Augmented Dickey-Fuller, unit-root test then the test returns the least number of differences required to pass the test at level \code{alpha}. After the unit-root test, seasonal test applies on the stationary \code{data}. \cr
#'   `period`     \tab The periodic nature of the seasonality. This can be a number indicating the number of observations in each seasonal period (for example, annual seasonality would be "1").  \cr
#'   `method`      \tab A string identifying method for seasonal decomposition. If NULL, "decomp" method is default. Possible values are c("none", "decomp", "stl", "stlplus", "tbats", "stR") phrases of methods denote. \cr
#'   `suroot_test`     \tab Type of seasonal unit root test to use. Possible values are "correlogram", "seas", "hegy", "ch" and "ocsb". \cr
#'   `suroot_tcrit`     \tab t-value for autocorrelogram.  \cr
#'   `suroot_uroot`     \tab If TRUE, unit root test for stationary before seasonal unit root test is allowed. \cr
#' }
#' }
#'
#' @return A model specification.
#'
#' @importFrom fabletools new_model_class new_model_definition
#' @importFrom rlang enquo
#'
#' @export
ATAM <- function(formula, ...){
        atam_model <- fabletools::new_model_class("ATAM", train = train_ata, specials = specials_ata)
        fabletools::new_model_definition(atam_model, !!rlang::enquo(formula))
}

#' Forecast a model from the fable ATA model
#'
#'
#' @param object The time series model used to produce the forecasts
#' @param new_data A `tsibble` containing future information used to forecast.
#' @param h The forecast horison (can be used instead of `new_data` for regular time series with no exogenous regressors).
#' @param ci_level Confidence Interval levels for forecasting. Default value is 95.
#' @param negative_forecast Negative values are allowed for forecasting. Default value is TRUE. If FALSE, all negative values for forecasting are set to 0.
#' @param ... Other arguments
#'
#' @return A vector of fitted residuals.
#'
#' @importFrom tsibble is_tsibble as_tibble measured_vars
#' @importFrom rlang enquo expr_text
#' @importFrom tsbox ts_ts
#' @importFrom stats tsp ts
#' @importFrom ATAforecasting ATA.Forecast
#' @importFrom distributional dist_degenerate
#'
#' @export
forecast.fable_ata <- function(object, new_data, h=NULL, ci_level=95, negative_forecast=TRUE, ...){
  mdl <- object$model_output
  spec_mdl <- object$spec
  if(!is.null(h) && !is.null(new_data)){
      warning("Input forecast horizon 'h' will be ignored as 'new_data/test_set' has been provided.")
      h <- nrow(new_data)
  }
  if (is.null(h)){
    if (spec_mdl$period==4){
      h <- 8
    }else if (spec_mdl$period==5){
      h <- 10
    }else if (spec_mdl$period==7){
      h <- 14
    }else if (spec_mdl$period==12){
      h <- 24
    }else if (spec_mdl$period==24){
      h <- 48
    }else {
      h <- 6
    }
    warning(paste("Input forecast horizon has been changed with ", h))
  }

  # Prepare data
  if (tsibble::is_tsibble(new_data)){
      test_set <- tsibble::as_tibble(new_data)[c(rlang::expr_text(index(new_data)), tsibble::measured_vars(new_data))]
      colnames(test_set) <- c("ds", "yh")
      pre_data <- ts_ts(new_data)
      test_set <- stats::ts(pre_data, frequency = spec_mdl$period, start=tsp(pre_data)[1])
  }else{
    test_set = new_data
  }

  fc <- safely(quietly(ATAforecasting::ATA.Forecast))(mdl, h, test_set, ci_level, negative_forecast)
  # Return forecasts
  distributional::dist_degenerate(fc$forecast)
}

#' Extract fitted values
#'
#' Extracts the fitted values from an estimated ATA model.
#'
#' @param object An estimated model.
#' @param ... Unused.
#'
#' @return A vector of fitted values.
#'
#' @export
fitted.fable_ata <- function(object, ...){
  object$est[[".fitted"]]
}

#' Extract model residuals
#'
#' Extracts the residuals from an estimated ATA model.
#'
#' @param object An estimated model.
#' @param ... Unused.
#'
#' @return A vector of residuals.
#'
#' @export
residuals.fable_ata <- function(object, ...){
  object$est[[".resid"]]
}


#' Extract estimated states from an ATA model.
#'
#' @param object An estimated model.
#' @param ... Unused.
#'
#' @return A [fabletools::dable()] containing estimated states.
#'
#' @importFrom fabletools as_dable
#' @importFrom rlang sym
#'
#' @export
components.fable_ata <- function(object, ...){
  cmp <- object$components
  three_cmp <- list("level" = cmp$level, "trend" = cmp$trend, "remainder" = cmp$remainder)
  seas_cmp <- cmp$seasonal
  fabletools::as_dable(cmp,
            resp = !!sym(response),
            method = model_sum(object),
            seasons = seas_cmp,
            aliases = three_cmp
           )
}

#' Glance an ATA model
#'
#' @param x An estimated model.
#' @param ... Unused.
#'
#' @return A one row tibble summarising the model's fit.
#'
#' @export
glance.fable_ata <- function(x, ...){
  x$fit
}

#' Tidy a ATA model
#'
#' @param x An estimated model.
#' @param ... Unused.
#'
#' @return The model's coefficients in a `tibble`.
#'
#' @export
tidy.fable_ata <- function(x, ...){
  x$par
}

#' Summary of ATA model
#'
#' @param x An estimated model.
#' @param ... Unused.
#'
#' @return The model's summary specs in a `tibble`.
#'
#' @export
model_sum.fable_ata <- function(x, ...){
  x$spec[["method"]]
}

#' Format of ATA model
#'
#' @param x An estimated model.
#' @param ... Unused.
#'
#' @export
format.fable_ata <- function(x, ...){
  "ATA"
}

#' Specialized Screen Print Function of ATA model
#'
#' @param object An estimated model.
#' @param ... Unused.
#'
#' @export
report.fable_ata <- function(object, ...) {
  x <- object$model_output
  opscipen <- options("scipen" = 100, "digits"=7)
    on.exit(options(opscipen))
	x <- object
    cat(x$method,"\n\n")
    if (x$level.fixed==TRUE){
      cat("   level.fixed: TRUE","\n\n")
    }
    if (x$trend.opt!="none"){
      cat(paste("   trend optimization method: trend.", x$trend.opt, "\n\n", sep=""))
    }
    if(!is.null(x$transform.method)){
      cat(paste("   '",x$transform.method, "' transformation method was selected.","\n\n", sep=""))
    }
    if(!is.null(x$lambda)){
      cat("   Box-Cox transformation: lambda=",round(x$lambda,4), "\n\n")
    }
    cat(paste("   model.type:",x$model.type, "\n\n"))
    if (x$is.season==FALSE){
      cat("   seasonal.model: no seasonality","\n\n")
    }else {
      cat(paste("   seasonal.model:",x$seasonal.model, "\n\n"))
    }
    if (x$is.season==TRUE){
      cat(paste("   seasonal.type:",x$seasonal.type, "\n\n"))
    }
    cat(paste("   forecast horizon:",x$h, "\n\n"))
    cat(paste("   accuracy.type:",x$accuracy.type, "\n\n"))

    cat("In-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$fits$sigma2, x$accuracy$fits$loglik, x$accuracy$MAE$inSample$MAE, x$accuracy$MSE$inSample$MSE, x$accuracy$MSE$inSample$RMSE, x$accuracy$MPE$inSample$MPE, x$accuracy$MAPE$inSample$MAPE, x$accuracy$sMAPE$inSample$sMAPE, x$accuracy$MASE$inSample$MASE, x$accuracy$OWA$inSample$OWA)
    names(stats) <- c("sigma2", "loglik", "MAE", "MSE", "RMSE", "MPE", "MAPE", "sMAPE", "MASE", "OWA")
    cat("\n")
    print(stats)
    cat("\n")

    cat("In-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$MAE$inSample$MdAE, x$accuracy$MSE$inSample$MdSE, x$accuracy$MSE$inSample$RMdSE, x$accuracy$MPE$inSample$MdPE, x$accuracy$MAPE$inSample$MdAPE, x$accuracy$sMAPE$inSample$sMdAPE)
    names(stats) <- c("MdAE", "MdSE", "RMdSE", "MdPE", "MdAPE", "sMdAPE")
    cat("\n")
    print(stats)
    cat("\n")


    cat("Out-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$MAE$outSample$MAE, x$accuracy$MSE$outSample$MSE, x$accuracy$MSE$outSample$RMSE, x$accuracy$MPE$outSample$MPE, x$accuracy$MAPE$outSample$MAPE, x$accuracy$sMAPE$outSample$sMAPE, x$accuracy$MASE$outSample$MASE, x$accuracy$OWA$outSample$OWA)
    names(stats) <- c("MAE", "MSE", "RMSE", "MPE", "MAPE", "sMAPE", "MASE",  "OWA")
    cat("\n")
    print(stats)
    cat("\n")

    cat("Out-Sample Accuracy Measures:","\n")
    stats <- c(x$accuracy$MAE$outSample$MdAE, x$accuracy$MSE$outSample$MdSE, x$accuracy$MSE$outSample$RMdSE, x$accuracy$MPE$outSample$MdPE, x$accuracy$MAPE$outSample$MdAPE, x$accuracy$sMAPE$outSample$sMdAPE)
    names(stats) <- c("MdAE", "MdSE", "RMdSE", "MdPE", "MdAPE", "sMdAPE")
    cat("\n")
    print(stats)
    cat("\n")

    cat("Information Criteria:","\n")
    stats <- c(x$accuracy$fits$AIC, x$accuracy$fits$AICc, x$accuracy$fits$BIC)
    names(stats) <- c("AIC", "AICc", "BIC")
    cat("\n")
    print(stats)
    cat("\n")

    stats <- c(x$execution.time[1], x$execution.time[2], x$execution.time[3])
    names(stats) <- c("user","system","elapsed")
    cat("\n")
    print(stats)
    cat("\n")
    cat(paste("calculation.time:",x$calculation.time, "\n\n"))
    cat("\n")

    cat("Forecasts:","\n")
    print(x$forecast)
    cat("\n\n")
}
