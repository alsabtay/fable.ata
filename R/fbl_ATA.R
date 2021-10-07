#' @docType package
#' @keywords package
"_PACKAGE"

globalVariables(c(".","self", "origin"))


no_xreg <- function(...) {
  abort("Exogenous regressors are not supported for ATA method.")
}

train_ATA <- function(.data, specials, opt_crit, nmse, ic, ...){
  if(length(tsibble::measured_vars(.data)) > 1){
    abort("Only univariate responses are supported by ATA Method")
  }
  # if (class(.data)[1]!="ts" & class(.data)[1]!="msts" & !is_tsibble(.data)){
  #   abort("Class of data must be ts/msts/tsibble object with single/multiple period seasonality. ATA Method was terminated!")
  # }
  # Prepare data for modelling
  if (is.ts(.data)){
    model_data = .data
    period <- stats::frequency(model_data)
    if (any(is.na(model_data))) {
      abort("Ata method does not support missing values.")
    }
  }else{
    model_data <- tsibble::as_tibble(.data)[c(rlang::expr_text(index(.data)), measured_vars(.data))]
    colnames(model_data) <- c("ds", "y")
    if (any(is.na(model_data$y))) {
      abort("Ata method does not support missing values.")
    }
    period <- fabletools::get_frequencies(period, .data, .auto = "smallest")
    pre_data <- tsbox::ts_ts(.data)
    model_data <- stats::ts(model_data$y, frequency = unname(period), start=tsp(pre_data)[1])
  }

  level <- specials$level
  trend <- specials$trend
  season <- specials$season
  transform <- specials$transform
  holdout <- specials$holdout

  seas_opt_crit <- ATAforecasting::ATA.SeasAttr(suroot.test=season$suroot_test, corrgram.tcrit = season$suroot_tcrit, uroot.pkg="tseries",
                                uroot.test=trend$uroot_test, uroot.type=trend$uroot_type, uroot.alpha=trend$uroot_alpha, uroot.maxd=1)
  bc_opt_crit <- ATAforecasting::ATA.BoxCoxAttr(bcMethod = transform$bcMethod, bcLower = transform$bcLower, bcUpper = transform$bcUpper)

  # Build and train model
  mdl_ATA <- ATAforecasting::ATA(X = model_data,
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
                  accuracy.type = opt_crit,
                  nmse = nmse,
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
  structure(
    list(
        par = tsibble::tibble(term = names(mdl_ATA$par.specs), estimate = unname(mdl_ATA$par.specs)),
        est = list(.fitted = mdl_ATA$fitted,
                   .resid = mdl_ATA$residuals),
        fit = mdl_ATA$accuracy$fits,
        components = list(coefp = mdl_ATA$coefp,
                          coefq = mdl_ATA$coefp,
                          P = mdl_ATA$parP,
                          Q = mdl_ATA$parQ,
                          PHI = mdl_ATA$parPHI,
                          response = mdl_ATA$actual,
                          level = mdl_ATA$level,
                          trend = mdl_ATA$trend,
                          seasonal = mdl_ATA$seasonal,
                          seasindex = mdl_ATA$seasonal.index,
                          seasadj = mdl_ATA$seasonal.adjusted,
                          remainder = mdl_ATA$residuals),
        transform = list(method = transform$method,
                         order = transform$order,
                         lambda = mdl_ATA$lambda,
                         shift = mdl_ATA$shift),
        holdout = list(holdout = mdl_ATA$holdout,
                       adjustment = mdl_ATA$holdout.adjustedP),
        spec = c(errortype = "A",
                 trendtype = mdl_ATA$model.type,
                 seasontype = mdl_ATA$seasonal.model,
                 damped = ifelse(mdl_ATA$parPHI==1, FALSE, TRUE),
                 period = mdl_ATA$seasonal.period,
                 method = mdl_ATA$method),
        model_output = mdl_ATA,
        class = "fbl_ATA")
      )
}

specials_ATA <- fabletools::new_specials(
   level = function(parP= NULL, level_fixed = FALSE, initial_level = FALSE){
                   as.list(environment())
                  },
   trend = function(type = c("N", "A", "M", "Ad", "Md"), parQ = NULL,
                   parPHI = NULL, parPHI_range = c(0.8, 1.0), parPHI_increment = 0.02,
                   uroot_test = "adf", uroot_alpha = 0.05, uroot_type = "level",
                   initial_trend = FALSE, opt_trend = "none"){

                     if (!all(is.element(type, c("N", "A", "M", "Ad", "Md")))) {
                       stop("Invalid trend type")
                     }
                     if (!all(is.element(uroot_test, c("adf", "pp", "kpss")))) {
                       stop("Invalid unit root test type")
                     }
                     if (!all(is.element(uroot_type, c("level", "trend")))) {
                       stop("Invalid specification of the deterministic component in the regression for unit root test")
                     }
                     if (!all(is.element(opt_trend, c("none","fixed","search")))) {
                       stop("Invalid trend optimization method")
                     }
                     if (parPHI_range[1] > parPHI_range[2]) {
                       abort("Lower PHI limits must be less than upper limits")
                     }
                     if (parPHI_increment<=0 | parPHI_increment>1){
                       abort("Increment size must be higher than 0 and less than 1")
                     }
                     type <- match.arg(type)
                     if (type == "N"){
                       parQ = 0
                       warn("Q has been set 0 because of no trend option.")
                     }
                     if ((type == "Ad" | type == "Md") & parPHI == 1){
                       parPHI = NULL
                       warn("PHI has been set NULL as damped trend is choosen.")
                     }
                     as.list(environment())
                  },
   season = function(type = c("N", "A", "M"), test = TRUE, period = NULL,
                  suroot_test = "correlogram", suroot_tcrit = 1.28,
                  method = "decomp"){
                    if (!all(is.element(type, c("N", "A", "M")))) {
                      stop("Invalid season type")
                    }
                    if (!all(is.element(suroot_test, c("correlogram", "seas", "hegy", "ch", "ocsb")))) {
                      stop("Invalid seasonal unit root test type")
                    }
                    if (!all(is.element(method, c("decomp", "stl", "stlplus", "stR", "tbats")))) {
                      stop("Invalid seasonal decompostion method")
                    }
                    type <- match.arg(type)
                    if (type == "N") {
                      period <- 1
                    }
                    as.list(environment())
                  },
   transform = function(method="none", order = "none", lambda = NULL, shift = 0,
                      bcMethod = "guerrero", bcLower = 0, bcUpper = 5){
                        if (!all(is.element(bcMethod, c("guerrero", "loglik")))) {
                          stop("Invalid method for calculating Box-Cox's lambda")
                        }
                        if (!all(is.element(order, c("none","before","after")))) {
                          stop("Invalid transformation order type")
                        }
                        if (!all(is.element(method, c("Box_Cox", "Sqrt", "Reciprocal", "Log", "NegLog", "Modulus", "BickelDoksum", "Manly", "Dual", "YeoJohnson", "GPower", "GLog")))) {
                          stop("Invalid power transformation method")
                        }
                        as.list(environment())
                    },
   holdout = function(holdout = FALSE, adjustment = TRUE, set_size = NULL ){
                      if (min(set_size)<=0){
                        abort("Percentage of holdout set size must be higher than 0.")
                      }
                      as.list(environment())
                    },
   xreg_specials = NULL ,
  .required_specials = c("trend", "season")
)


#' @export
ATA <- function(formula,
                    opt_crit = "sMAPE",
                    nmse = 3,
                    ic = "AIC", ...){
        ATA_model <- fabletools::new_model_class(model = "ATA",
                                                     train = train_ATA,
                                                     specials = specials_ATA,
                                                     check = all_tsbl_checks_ata
                                                     )
        fabletools::new_model_definition(ATA_model,
                                         !!rlang::enquo(formula),
                                         opt_crit = opt_crit,
                                         nmse = nmse,
                                         ic = ic, ...)
}


#' @export
forecast.fbl_ATA <- function(object, new_data, h=NULL, ci_level=95, negative_forecast=TRUE, ...){
  mdl <- object$model_output
  spec_mdl <- object$spec
  if(!is.null(h) && !is.null(new_data)){
      warn("Input forecast horizon 'h' will be ignored as 'new_data/test_set' has been provided.")
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
    warn(paste("Input forecast horizon has been changed with ", h))
  }

  # Prepare data
  if (is_tsibble(new_data)){
      test_set <- as_tibble(new_data)[c(expr_text(index(new_data)), measured_vars(new_data))]
      colnames(test_set) <- c("ds", "yh")
      pre_data <- ts_ts(new_data)
      test_set <- ts(pre_data, frequency = spec_mdl$period, start=tsp(pre_data)[1])
  }else{
    test_set = new_data
  }

  fc <- ATAforecasting::ATA.Forecast(mdl, h, test_set, ci_level, negative_forecast)$forecast
  # Return forecasts
  distributional::dist_degenerate(fc$forecast)
}


#' @export
fitted.fbl_ATA <- function(object, ...){
  object$est[[".fitted"]]
}


#' @export
residuals.fbl_ATA <- function(object, ...){
  object$est[[".resid"]]
}


#' @export
components.fbl_ATA <- function(object, ...){
  cmp <- object$spec
  fabletools::as_dable(cmp,
            resp = !!sym("response"),
            method = "ATA",
            seasons = !!sym("seasonal"),
            aliases = list(level = level, trend = trend, remainder = remainder)
           )
}


#' @export
glance.fbl_ATA <- function(x, ...){
  x$fit
}

#' @export
tidy.fbl_ATA <- function(x, ...){
  x$par
}

#' @export
model_sum.fbl_ATA <- function(x){
  x$spec[["method"]]
}

#' @export
format.fbl_ATA <- function(x, ...){
  "Ata Method"
}


#' @export
report.fbl_ATA <- function(object, ...) {
  ATAforecasting::ATA.print(object$model_output)
}
