

#' List wrappers in survSuperLearner
#'
#' This function lists all prediction algorithms and/or screening algorithms built in to the \code{survSuperLearner} package.
#'
#' @usage \code{survlistWrappers()}
#'
#' @param what If \code{"both"}, lists both prediction and screening functions; if \code{"survSL"}, lists only prediction functions; if \code{"survscreen"}, lists only screening functions; otherwise all functions in the package.
#' @return Invisible character vector with the requested functions.

survlistWrappers <- function (what = "both") {
  everything <- sort(getNamespaceExports("survSuperLearner"))
  if (what == "both") {
    message("All prediction algorithm wrappers in survSuperLearner:\n")
    print(everything[grepl(pattern = "^surv[S]L", everything)])
    message("\nAll screening algorithm wrappers in survSuperLearner:\n")
    print("All")
    print(everything[grepl(pattern = "survscreen", everything)])
  }
  else if (what == "survSL") {
    message("All prediction algorithm wrappers in survSuperLearner:\n")
    print(everything[grepl(pattern = "^surv[S]L", everything)])
  }
  else if (what == "survscreen") {
    message("All screening algorithm wrappers in survSuperLearner:\n")
    print("All")
    print(everything[grepl(pattern = "survscreen", everything)])
  }
  else {
    message("All functions in survSuperLearner:\n")
    print(everything)
  }
  invisible(everything)
}

#' Wrapper functions for prediction algorithms in survSuperLearner
#'
#' This is a template function for a \code{survSuperLearner} prediction algorithm. You can use this to write your own prediction algorithms.
#'
#' @param time Observed follow-up time; i.e. minimum of the event and censoring times.
#' @param event Observed event indicator; i.e, whether the follow-up time corresponds to an event or censoring.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction. Should have the same variable names and structure as \code{X}.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Observation clusters.
#' @param ... Additional ignored arguments.
#' @return \item{pred}{Matrix of predictions, with the same number of rows as \code{newX} and number of columns equal to the length of \code{new.times}. Rows index new observations, and columns index new times at which the survival was computed.}
#' \item{fit}{Elements of the model fit. Typically includes \code{fit$object}, the fitted object, but may include additional elements obtained in the fitting process that would be needed for prediction.}

survSL.template <- function(time, event, X, newX, new.times, obsWeights, id, ...) {
  pred <- numeric()
  fit <- vector("list", length = 0)
  class(fit) <- c("survSL.template")
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Wrapper functions for screening algorithms in survSuperLearner
#'
#' This is a template function for a \code{survSuperLearner} screening algorithm. You can use this to write your own screening algorithms.
#'
#' @param time Observed follow-up time; i.e. minimum of the event and censoring times.
#' @param event Observed event indicator; i.e, whether the follow-up time corresponds to an event or censoring.
#' @param X Training covariate data.frame.
#' @param obsWeights Observation weights.
#' @param id Observation clusters.
#' @param ... Additional ignored arguments.
#' @return Logical vector of the same length as the number of columns as \code{X} indicating which variables should be included.

survscreen.template <- function(time, event, X, obsWeights, id, ...) {
  whichScreen <- rep(TRUE, ncol(X))
  return(whichScreen)
}

#' Wrapper function for Kaplan-Meier prediction algorithm
#'
#' This prediciton algorithm ignores all covariates and simply computes the Kaplan-Meier estimator of the marginal survival function of the event as indicated by the right-censored data \code{time} and \code{event} using the \code{\link[survival]{survfit}} function.
#'
#' @param time Observed follow-up time; i.e. minimum of the event and censoring times.
#' @param event Observed event indicator; i.e, whether the follow-up time corresponds to an event or censoring.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction. Should have the same variable names and structure as \code{X}.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param ... Additional ignored arguments.
#' @return \item{pred}{Matrix of predictions, with the same number of rows as \code{newX} and number of columns equal to the length of \code{new.times}. Rows index new observations, and columns index new times at which the survival was computed.}
#' \item{fit}{One-element list including \code{object}, the fitted \code{\link[survival]{survfit}} object.}


survSL.km <- function(time, event, X, newX, new.times, obsWeights, ...) {
  fit.km <- survival::survfit(
    survival::Surv(time, event)~1,
    weights = obsWeights
  )
  pred <- matrix(stats::stepfun(fit.km$time, c(1,fit.km$surv), right = FALSE)(new.times),
                 nrow=nrow(newX), ncol = length(new.times), byrow=TRUE)

  fit <- list(object = fit.km)
  class(fit) <- c("survSL.km")
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Prediction functions for Kaplan-Meier prediction algorithm
#'
#' Obtains predicted survivals from a fitted \code{survSL.km} object.
#'
#' @param object Fitted \code{survSL.km} object.
#' @param newX New covariate data.frame for which to obtain predictions.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param ... Additional ignored arguments.
#' @return Matrix of predictions, with the same number of rows as \code{newX} and number of columns equal to the length of \code{new.times}. Rows index new observations, and columns index new times at which the survival was computed.

predict.survSL.km <- function(object, newX, new.times, ...) {
  matrix(stats::stepfun(object$object$time, c(1,object$object$surv), right = FALSE)(new.times),
         nrow=nrow(newX), ncol = length(new.times), byrow=TRUE)
}

# est <- survSL.km(time, 1-event, X = data.frame(W), newX = data.frame(W), new.times = 0:15, obsWeights=NULL, id=NULL)
# pred <- predict.survSL.km(est$fit, newX = data.frame(W)[1:3,], new.times = c(1.5,2.5))

#' Wrapper function for Cox proportional hazards regression prediction algorithm
#'
#' This prediciton algorithm uses the partial maximum likelihood estimator of the coefficients and Breslow estimator of the baseline cumulative hazard in the Cox proportional hazards model using the \code{\link[survival]{coxph}} and \code{\link[survival]{survfit}} functions.
#'
#' @param time Observed follow-up time; i.e. minimum of the event and censoring times.
#' @param event Observed event indicator; i.e, whether the follow-up time corresponds to an event or censoring.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction. Should have the same variable names and structure as \code{X}.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param ... Additional ignored arguments.
#' @return \item{pred}{Matrix of predictions, with the same number of rows as \code{newX} and number of columns equal to the length of \code{new.times}. Rows index new observations, and columns index new times at which the survival was computed.}
#' \item{fit}{One-element list including \code{object}, the fitted \code{\link[survival]{coxph}} object.}


survSL.coxph <- function(time, event, X, newX, new.times, obsWeights, ...) {
  fit.coxph <- survival::coxph(
    survival::Surv(time, event) ~ .,
    data = as.data.frame(cbind(time=time, event=event, X)),
    weights = obsWeights
  )
  pred <- t(summary(survival::survfit(fit.coxph,
                                      newdata=newX,
                                      se.fit = FALSE,
                                      conf.int = FALSE),
                    times=new.times)$surv)
  if(ncol(pred) < length(new.times)) {
    pred <- cbind(pred, matrix(pred[,ncol(pred)], nrow=nrow(pred), ncol=length(new.times) - ncol(pred)))
  }

  fit <- list(object = fit.coxph)
  class(fit) <- c("survSL.coxph")
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Prediction functions for Cox regression prediction algorithm
#'
#' Obtains predicted survivals from a fitted \code{survSL.coxph} object.
#'
#' @param object Fitted \code{survSL.coxph} object.
#' @param newX New covariate data.frame for which to obtain predictions.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param ... Additional ignored arguments.
#' @return Matrix of predictions, with the same number of rows as \code{newX} and number of columns equal to the length of \code{new.times}. Rows index new observations, and columns index new times at which the survival was computed.

predict.survSL.coxph <- function(object, newX, new.times, ...) {
  pred <- t(summary(survival::survfit(formula = object$object,
                                      newdata = newX,
                                      se.fit = FALSE,
                                      conf.int = FALSE),
                    times=new.times)$surv)
  if(ncol(pred) < length(new.times)) {
    pred <- cbind(pred, matrix(pred[,ncol(pred)], nrow=nrow(pred), ncol=length(new.times) - ncol(pred)))
  }
  return(pred)
}

# est <- survSL.coxph(time, 1-event, X = data.frame(W), newX = data.frame(W), new.times = 0:15, obsWeights=NULL, id=NULL)
# pred <- predict.survSL.coxph(est$fit, newX = data.frame(W)[1:3,], new.times = c(1.5,2.5))

#' Wrapper function for Random Survival Forests prediction algorithm
#'
#' This prediciton algorithm uses the \code{\link[randomForestSRC]{rfsrc}} function from the \code{randomForestSRC} package to estimate a survival random forest.
#'
#' @param time Observed follow-up time; i.e. minimum of the event and censoring times.
#' @param event Observed event indicator; i.e, whether the follow-up time corresponds to an event or censoring.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction. Should have the same variable names and structure as \code{X}.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param id Currently ignored
#' @param ... Additional arguments passed on to \code{\link[randomForestSRC]{rfsrc}}.
#' @return \item{pred}{Matrix of predictions, with the same number of rows as \code{newX} and number of columns equal to the length of \code{new.times}. Rows index new observations, and columns index new times at which the survival was computed.}
#' \item{fit}{One-element list including \code{object}, the fitted \code{\link[randomForestSRC]{rfsrc}}. object.}
#' @references Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. \emph{The Annals of Applied Statistics}, 2(3), 841-860.

survSL.rfsrc <- function(time, event, X, newX, new.times, obsWeights, id, ...) {
  data <- data.frame(time, event)
  data <- cbind(data, X)
  fit.rfsrc <- randomForestSRC::rfsrc(Surv(time, event) ~ ., data=data, importance = FALSE, case.wt = obsWeights, ...)
  survs <- predict(fit.rfsrc, newdata=newX, importance='none')$survival
  pred <- t(sapply(1:nrow(survs), function(i) {
    stats::approx(c(0,fit.rfsrc$time.interest), c(1,survs[i,]), method='constant', xout = new.times, rule = 2)$y
  }))

  fit <- list(object = fit.rfsrc)
  class(fit) <- c("survSL.rfsrc")
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Prediction functions for survival random forest prediction algorithm
#'
#' Obtains predicted survivals from a fitted \code{survSL.rfsrc} object.
#'
#' @param object Fitted \code{survSL.rfsrc} object.
#' @param newX New covariate data.frame for which to obtain predictions.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param ... Additional ignored arguments.
#' @return Matrix of predictions, with the same number of rows as \code{newX} and number of columns equal to the length of \code{new.times}. Rows index new observations, and columns index new times at which the survival was computed.

predict.survSL.rfsrc <- function(object, newX, new.times, ...) {
  survs <- predict(object$object, newdata=newX, importance='none')$survival
  pred <- t(sapply(1:nrow(survs), function(i) {
    stats::approx(c(0,object$object$time.interest), c(1,survs[i,]), method='constant', xout = new.times, rule = 2)$y
  }))
  return(pred)
}

# est <- survSL.rfsrc(time, 1-event, X = data.frame(W), newX = data.frame(W), new.times = c(0:14, 14.999), obsWeights=rep(1, length(time)), id=NULL)
# pred <- predict.survSL.rfsrc(est$fit, newX = data.frame(W)[1:3,], new.times = c(1.5,2.5))

#' Wrapper function for generalized additive Cox regression prediction algorithm
#'
#' This prediciton algorithm uses the \code{\link[mgcv]{gam}} function from the \code{mgcv} package to estimate a generalized additive Cox proportional hazards regression model. This model generalizes the usual Cox proportional hazards model to allow for an additive combination of smooth and possibly non-linear functions of the continuous covariates.
#'
#' @param time Observed follow-up time; i.e. minimum of the event and censoring times.
#' @param event Observed event indicator; i.e, whether the follow-up time corresponds to an event or censoring.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction. Should have the same variable names and structure as \code{X}.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param cts.num The lower cutoff of unique values at which a covariate should be treated as continuous. Any covariate with number of unique values strictly larger than \code{cts.num} will be treated as continuous and receive a smooth term in the GAM. If \code{X} contains unordered factors with large numbers of unique values, \code{cts.num} should be set to larger than these numbers of unique values, otherwise an error may be thrown.
#' @param ... Additional ignored arguments.
#' @return \item{pred}{Matrix of predictions, with the same number of rows as \code{newX} and number of columns equal to the length of \code{new.times}. Rows index new observations, and columns index new times at which the survival was computed.}
#' \item{fit}{One-element list including \code{object}, the fitted \code{\link[mgcv]{gam}}. object.}


survSL.gam <- function(time, event, X, newX, new.times, cts.num = 5, ...) {

  if (any(round(list(...)[["obsWeights"]],8)!=1))
    warning("Argument 'obsWeights' is ignored by the 'gam' library algorithm")
  if ("gam" %in% loadedNamespaces())
    warning("mgcv and gam packages are both in use. You might see an error because both packages use the same function names.")
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  if (sum(!cts.x) > 0) {
    gam.model <- as.formula(paste("time~", paste(paste("s(",
                                                    colnames(X[, cts.x, drop = FALSE]),
                                                    ")", sep = ""), collapse = "+"), "+", paste(colnames(X[,
                                                                                                           !cts.x, drop = FALSE]), collapse = "+")))
  } else {
    gam.model <- as.formula(paste("time~", paste("s(",
                                                    colnames(X[, cts.x, drop = FALSE]),
                                                    ")", sep = "", collapse = "+")))
  }
  if (sum(!cts.x) == length(cts.x)) {
    gam.model <- as.formula(paste("time~", paste(colnames(X),
                                              collapse = "+"), sep = ""))
  }

  fit.gam <- mgcv::gam(gam.model, family=mgcv::cox.ph(), data = X, weights=event)

  new.data <- data.frame(time=rep(new.times, each = nrow(newX)))
  for (col in names(newX)) new.data[[col]] <- rep(newX[[col]], length(new.times))
  pred <- predict(fit.gam, newdata=new.data, type="response", se=FALSE)
  pred <- matrix(pred, nrow = nrow(newX), ncol = length(new.times))

  fit <- list(object = fit.gam)
  class(fit) <- c("survSL.gam")
  out <- list(pred = pred, fit = fit)
  return(out)
}


#' Prediction functions for generalized additive Cox regression prediction algorithm
#'
#' Obtains predicted survivals from a fitted \code{survSL.gam} object.
#'
#' @param object Fitted \code{survSL.gam} object.
#' @param newX New covariate data.frame for which to obtain predictions.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param ... Additional ignored arguments.
#' @return Matrix of predictions, with the same number of rows as \code{newX} and number of columns equal to the length of \code{new.times}. Rows index new observations, and columns index new times at which the survival was computed.

predict.survSL.gam <- function(object, newX, new.times, ...) {
  new.data <- data.frame(time=rep(new.times, each = nrow(newX)))
  for (col in names(newX)) new.data[[col]] <- rep(newX[[col]], length(new.times))
  pred <- predict(object$object, newdata=new.data, type="response", se=FALSE)
  pred <- matrix(pred, nrow = nrow(newX), ncol = length(new.times))
  return(pred)
}

# est <- survSL.gam(time, event, X = X, newX = X, new.times = c(0:14, 14.999), obsWeights=rep(1, length(time)), id=NULL)
# pred <- predict.survSL.gam(est$fit, newX = data.frame(W)[1:3,], new.times = c(0:14, 14.999))

#' Wrapper function for zero-inflated parametric survival regression prediction algorithms
#'
#' These prediciton algorithm use the \code{\link[survival]{survreg}} function from the \code{survival} package to estimate parametric survival regressions. See details for specific parametric models.
#'
#' @param time Observed follow-up time; i.e. minimum of the event and censoring times.
#' @param event Observed event indicator; i.e, whether the follow-up time corresponds to an event or censoring.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction. Should have the same variable names and structure as \code{X}.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param ... Additional ignored arguments.
#' @details Parametric regression models assume a particular parametric form for the distribution of the event given covariates, where the covariates contribute in a linear way to one of the parameters of the distribution. Implemented here are currently exponential (\code{survSL.expreg}), Weibull (\code{survSL.weibreg}), and log-logistic (\code{survSL.loglogreg}) regressions.
#'
#' Note that survival regressions typically assume that the distribution of the event is continuous and strictly positive. Therefore, they will throw errors if there are large discrete components of the observed distribution of the event. Since some survival outcomes have positive mass at zero, we have amended the standard survival regression to include a component at zero, making this a \emph{zero-inflated} regression model. Specifically, if there are observed event times equal to zero, then a preliminary logistic regression is fit to predict the probability that the event time is exactly zero, and the survival regression is fit to the strictly positive part of the distribution.
#' @return \item{pred}{Matrix of predictions, with the same number of rows as \code{newX} and number of columns equal to the length of \code{new.times}. Rows index new observations, and columns index new times at which the survival was computed.}
#' \item{fit}{Two-element list including \code{reg.object}, the fitted \code{\link[survival]{survreg}} object, and \code{pos.object}, the fitted \code{\link[stats]{glm}} object for the probability that the event was positive (or 1 if no zeroes were detected).}

survSL.expreg <- function(time, event, X, newX, new.times, obsWeights, ...) {

  if(any(time == 0 & event == 1)) {
    timepos <- as.numeric(time > 0 & event == 1)
    fit.pos <- stats::glm(timepos ~ ., data=cbind(timepos, X)[event == 1,], family='binomial', weights = obsWeights[event == 1])
    pos.pred <- predict(fit.pos, newdata = newX, type = 'response')
  } else {
    fit.pos <- 1
    pos.pred <- rep(1, nrow(newX))
  }

  fit.expreg <- survival::survreg(survival::Surv(time[time > 0], event[time > 0]) ~ .,
                                  data = X[time > 0,,drop=FALSE],
                                  weights = obsWeights[time > 0], dist = 'exponential')
  pred <- predict(fit.expreg, newdata = newX, type = 'quantile', p = seq(0, .999, by=.001))
  pred <- try(t(sapply(1:nrow(pred), function(j) {
    pos.pred[j] * (1-stats::approx(pred[j,], seq(0, .999, by=.001), xout = new.times, method = 'linear', rule = 2)$y)
  })), silent = TRUE)
  if(inherits(pred, "try-error")) stop("Survival regression failed to produce predictions.")

  fit <- list(reg.object = fit.expreg, pos.object = fit.pos)
  class(fit) <- c("survSL.expreg")
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Prediction functions for parametric survival model prediction algorithms
#'
#' Obtains predicted survivals from a fitted parametric survival regression object.
#'
#' @param object Fitted \code{survSL.expreg}, \code{survSL.weibreg}, etc. object.
#' @param newX New covariate data.frame for which to obtain predictions.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param ... Additional ignored arguments.
#' @return Matrix of predictions, with the same number of rows as \code{newX} and number of columns equal to the length of \code{new.times}. Rows index new observations, and columns index new times at which the survival was computed.

predict.survSL.expreg <- function(object, newX, new.times, ...) {

  if(inherits(object$pos.object, "glm")) {
    pos.pred <- predict(object$pos.object, newdata = newX, type = 'response')
  } else {
    pos.pred <- rep(1, nrow(newX))
  }

  pred <- predict(object$reg.object, newdata = newX, type = 'quantile', p = seq(0, .999, by=.001))
  pred <- t(sapply(1:nrow(pred), function(j) {
    pos.pred[j] * (1-stats::approx(pred[j,], seq(0, .999, by=.001), xout = new.times, method = 'linear', rule = 2)$y)
  }))

  return(pred)
}

# est <- survSL.expreg(time, 1-event, X = data.frame(W), newX = data.frame(W), new.times = c(0:14, 14.999), obsWeights=rep(1, length(time)), id=NULL)
# pred <- predict.survSL.expreg(est$fit, newX = data.frame(W)[1:3,], new.times = c(0:14, 14.999))

#' @rdname survSL.expreg
survSL.weibreg <- function(time, event, X, newX, new.times, obsWeights, id, ...) {

  if(any(time == 0 & event == 1)) {
    timepos <- as.numeric(time > 0 & event == 1)
    fit.pos <- stats::glm(timepos ~ ., data=cbind(timepos, X)[event == 1,], family='binomial', weights = obsWeights[event == 1])
    pos.pred <- predict(fit.pos, newdata = newX, type = 'response')
  } else {
    fit.pos <- 1
    pos.pred <- rep(1, nrow(newX))
  }

  fit.weibreg <- survival::survreg(survival::Surv(time[time > 0], event[time > 0]) ~ .,
                                  data = X[time > 0,,drop=FALSE],
                                  weights = obsWeights[time > 0], dist = 'weibull')
  pred <- predict(fit.weibreg, newdata = newX, type = 'quantile', p = seq(0, .999, by=.001))
  pred <- try(t(sapply(1:nrow(pred), function(j) {
    pos.pred[j] * (1-stats::approx(pred[j,], seq(0, .999, by=.001), xout = new.times, method = 'linear', rule = 2)$y)
  })), silent=TRUE)
  if(inherits(pred, "try-error")) stop("Survival regression failed to produce predictions.")

  fit <- list(reg.object = fit.weibreg, pos.object = fit.pos)
  class(fit) <- c("survSL.weibreg")
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' @rdname predict.survSL.expreg
predict.survSL.weibreg <- function(object, newX, new.times, ...) {

  if(inherits(object$pos.object, "glm")) {
    pos.pred <- predict(object$pos.object, newdata = newX, type = 'response')
  } else {
    pos.pred <- rep(1, nrow(newX))
  }

  pred <- predict(object$reg.object, newdata = newX, type = 'quantile', p = seq(0, .999, by=.001))
  pred <- t(sapply(1:nrow(pred), function(j) {
    pos.pred[j] * (1-stats::approx(pred[j,], seq(0, .999, by=.001), xout = new.times, method = 'linear', rule = 2)$y)
  }))

  return(pred)
}

# est <- survSL.weibreg(time, 1-event, X = data.frame(W), newX = data.frame(W), new.times = c(0:14, 14.999), obsWeights=rep(1, length(time)), id=NULL)
# pred <- predict.survSL.weibreg(est$fit, newX = data.frame(W)[1:3,], new.times = c(0:14, 14.999))

#' @rdname survSL.expreg
survSL.loglogreg <- function(time, event, X, newX, new.times, obsWeights, id, ...) {

  if(any(time == 0 & event == 1)) {
    timepos <- as.numeric(time > 0 & event == 1)
    fit.pos <- stats::glm(timepos ~ ., data=cbind(timepos, X)[event == 1,], family='binomial', weights = obsWeights[event == 1])
    pos.pred <- predict(fit.pos, newdata = newX, type = 'response')
  } else {
    fit.pos <- 1
    pos.pred <- rep(1, nrow(newX))
  }

  fit.loglogreg <- survival::survreg(survival::Surv(time[time > 0], event[time > 0]) ~ .,
                                   data = X[time > 0,,drop=FALSE],
                                   weights = obsWeights[time > 0], dist = 'loglogistic')
  pred <- predict(fit.loglogreg, newdata = newX, type = 'quantile', p = seq(0, .999, by=.001))
  pred <- try(t(sapply(1:nrow(pred), function(j) {
    pos.pred[j] * (1-stats::approx(pred[j,], seq(0, .999, by=.001), xout = new.times, method = 'linear', rule = 2)$y)
  })), silent=TRUE)
  if(inherits(pred, "try-error")) stop("Survival regression failed to produce predictions.")

  fit <- list(reg.object = fit.loglogreg, pos.object = fit.pos)
  class(fit) <- c("survSL.loglogreg")
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' @rdname predict.survSL.expreg
predict.survSL.loglogreg <- function(object, newX, new.times, ...) {

  if(inherits(object$pos.object, "glm")) {
    pos.pred <- predict(object$pos.object, newdata = newX, type = 'response')
  } else {
    pos.pred <- rep(1, nrow(newX))
  }

  pred <- predict(object$reg.object, newdata = newX, type = 'quantile', p = seq(0, .999, by=.001))
  pred <- t(sapply(1:nrow(pred), function(j) {
    pos.pred[j] * (1-stats::approx(pred[j,], seq(0, .999, by=.001), xout = new.times, method = 'linear', rule = 2)$y)
  }))

  return(pred)
}

# est <- survSL.loglogreg(time, 1-event, X = data.frame(W), newX = data.frame(W), new.times = c(0:14, 14.999), obsWeights=rep(1, length(time)), id=NULL)
# pred <- predict.survSL.weibreg(est$fit, newX = data.frame(W)[1:3,], new.times = c(0:14, 14.999))

#' Wrapper function for piecewise constant hazard regression
#'
#' This prediciton algorithm uses the \code{\link[pch]{pchreg}} function from the \code{pch} package to estimate piecewise constant hazard regressions.
#'
#' @param time Observed follow-up time; i.e. minimum of the event and censoring times.
#' @param event Observed event indicator; i.e, whether the follow-up time corresponds to an event or censoring.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction. Should have the same variable names and structure as \code{X}.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param breaks Number or numeric vector of breaks to be passed to \code{\link[pch]{pchreg}}.
#' @param ... Additional ignored arguments.

survSL.pchreg <- function(time, event, X, newX, new.times, obsWeights, breaks = 4, ...) {

  fit.pchreg <- pch::pchreg(survival::Surv(time, event) ~ .,
                                  data = X,
                                  weights = obsWeights, breaks = breaks)
  pred <- try(sapply(new.times, function(t0) {
    predict(fit.pchreg, newdata = cbind(newX, time = t0), type = 'distr')$Surv
  }), silent = TRUE)
  if(inherits(pred, "try-error")) stop("PCH regression failed to produce predictions.")

  fit <- list(reg.object = fit.pchreg)
  class(fit) <- c("survSL.pchreg")
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Prediction functions for piecewise constant hazard regression
#'
#' Obtains predicted survivals from a fitted piecewise constant hazard regression object.
#'
#' @param object Fitted \code{survSL.pchreg} object.
#' @param newX New covariate data.frame for which to obtain predictions.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param ... Additional ignored arguments.
#' @return Matrix of predictions, with the same number of rows as \code{newX} and number of columns equal to the length of \code{new.times}. Rows index new observations, and columns index new times at which the survival was computed.

predict.survSL.pchreg <- function(object, newX, new.times, ...) {

  pred <- sapply(new.times, function(t0) {
    predict(object$fit$reg.object, newdata = cbind(newX, time = t0), type = 'distr')$Surv
  })
  return(pred)
}


#' Wrapper function for piecewise constant hazard SuperLearner
#'
#' This prediciton algorithm assumes that the hazard is constant on intervals, and estimates the hazard
#' within each interval conditional on covariates flexibly using SuperLearner.
#'
#' @param time Observed follow-up time; i.e. minimum of the event and censoring times.
#' @param event Observed event indicator; i.e, whether the follow-up time corresponds to an event or censoring.
#' @param X Training covariate data.frame.
#' @param newX Test covariate data.frame to use for prediction. Should have the same variable names and structure as \code{X}.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param obsWeights Observation weights.
#' @param breaks Either the number of intervals to use or the endpoints of the intervals.
#' @param SL.library Library to use for SuperLearning of individual bins.
#' @param ... Additional ignored arguments.
#'
survSL.pchSL <- function(time, event, X, newX, new.times, obsWeights, breaks, SL.library, ...) {
  if(length(breaks) == 1) {
    n.intervals <- breaks
    breaks <- c(0,as.numeric(stats::quantile(time[event == 1], probs = seq(0,1,by=1/n.intervals)[-1])))
  } else {
    n.intervals <- length(breaks) - 1
  }
  #breaks[length(breaks)] <- Inf
  intervals <- cbind(breaks[-length(breaks)], breaks[-1])

  sl.fits <- lapply(1:n.intervals, function(j) {
    if(j == 1) samp <- time >= 0
    else samp <- time > intervals[j,1]
    if(j < n.intervals) outcome <- as.numeric(time <= intervals[j,2] & event == 1)
    else outcome <- event
    fit <- try(SuperLearner::SuperLearner(Y = outcome[samp], X = X[samp,,drop=FALSE], newX = newX, SL.library = SL.library, family = 'binomial', method = 'method.NNloglik', obsWeights = obsWeights[samp]), silent=TRUE)
    if(inherits(fit, "try-error")) {
      fit <- try(SuperLearner::SuperLearner(Y = outcome[samp], X = X[samp,,drop=FALSE], newX = newX, SL.library = SL.library, family = 'binomial', method = 'method.NNLS', obsWeights = obsWeights[samp]), silent=TRUE)
      if(inherits(fit, "try-error")) {
        fit <- try(SuperLearner::SuperLearner(Y = outcome[samp], X = X[samp,,drop=FALSE], newX = newX, SL.library = SL.library, family = 'binomial', method = 'method.NNLS', obsWeights = obsWeights[samp]), silent=TRUE)
        if(inherits(fit, "try-error")) stop("Error in computing bin SuperLearner.")
      }
    }
    fit
  })

  hazard.ests <- sapply(1:n.intervals, function(j) {
    sl.fits[[j]]$SL.predict / (intervals[j,2] - intervals[j,1])
  })

  cum.hazard.ests <- t(apply(hazard.ests, 1, function(row) {
    cumsum(row * diff(breaks))
  }))

  new.time.bins <- findInterval(new.times, breaks, all.inside = TRUE)

  pred <- sapply(1:length(new.times), function(j) {
    bin <- new.time.bins[j]
    if(bin > 1) base <- cum.hazard.ests[,bin-1]
    else base <- 0
    exp(-(base + hazard.ests[,bin] * (new.times[j] - breaks[bin])))
  })

  fit <- list(sl.fits = sl.fits, breaks=breaks)
  class(fit) <- c("survSL.pchSL")
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Prediction functions for piecewise constant hazard SuperLearner prediction algorithm
#'
#' Obtains predicted survivals from a fitted \code{survSL.pchSL} object.
#'
#' @param object Fitted \code{survSL.pchSL} object.
#' @param newX New covariate data.frame for which to obtain predictions.
#' @param new.times Times at which to obtain to obtain the predicted survivals.
#' @param ... Additional ignored arguments.
#' @return Matrix of predictions, with the same number of rows as \code{newX} and number of columns equal to the length of \code{new.times}. Rows index new observations, and columns index new times at which the survival was computed.

predict.survSL.pchSL <- function(object, newX, new.times, ...) {
  breaks <- object$fit$breaks
  n.intervals <- length(breaks) - 1
  intervals <- cbind(breaks[-length(breaks)], breaks[-1])

  hazard.ests <- sapply(1:n.intervals, function(j) {
    predict(object$fit$sl.fits[[j]], newdata=newX)$pred / (intervals[j,2] - intervals[j,1])
  })

  cum.hazard.ests <- t(apply(hazard.ests, 1, function(row) {
    cumsum(row * diff(breaks))
  }))

  new.time.bins <- findInterval(new.times, breaks, all.inside = TRUE)

  pred <- sapply(1:length(new.times), function(j) {
    bin <- new.time.bins[j]
    if(bin > 1) base <- cum.hazard.ests[,bin-1]
    else base <- 0
    exp(-(base + hazard.ests[,bin] * (new.times[j] - breaks[bin])))
  })

  return(pred)
}


#' Wrapper function for glmnet screening algorithm
#'
#' This screening algorithm uses the \code{\link[glmnet]{glmnet}} function from the \code{glmnet} package to select covariates.
#'
#' @param time Observed follow-up time; i.e. minimum of the event and censoring times.
#' @param event Observed event indicator; i.e, whether the follow-up time corresponds to an event or censoring.
#' @param X Training covariate data.frame.
#' @param obsWeights Observation weights.
#' @param alpha Penalty exponent for \code{glmnet}. Defaults to 1 (lasso penalty).
#' @param minscreen Minimum number of covariates to return. Defaults to 2.
#' @param nfolds Number of folds for cross-validation selection of penalty parameter. Defaults to 10.
#' @param nlambda Number of penalty parameters to search over. Defaults to 100.
#' @param ... Additional ignored arguments.
#' @details The penalty parameter is selected using cross-validation via \code{\link[glmnet]{cv.glmnet}}. If this results in fewer than \code{minscreen} covariates, the penalty is increased to include \code{minscreen} covariates.
#' @return Logical vector of the same length as the number of columns of \code{X} indicating which variables were included.

survscreen.glmnet <- function(time, event, X, obsWeights, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, ...) {
  if (!is.matrix(X)) {
    X <- stats::model.matrix(~-1 + ., X)
  }
  time[event == 0] <- time[event == 0] + min(diff(sort(unique(time)))) / 2
  if(any(time == 0)) time[time == 0] <- min(time[time > 0]) / 2
  fit.glmnet <- glmnet::cv.glmnet(y = survival::Surv(time, event), x = X,
                                  weights = obsWeights, family = 'cox', alpha = alpha, nfolds = nfolds, nlambda = nlambda)

  whichVariable <- (as.numeric(coef(fit.glmnet$glmnet.fit, s = fit.glmnet$lambda.min)) != 0)
  if (sum(whichVariable) < minscreen) {
    warning("fewer than minscreen variables passed the glmnet screen, increased lambda to allow minscreen variables")
    sumCoef <- apply(as.matrix(fit.glmnet$glmnet.fit$beta), 2, function(x) sum((x != 0)))
    newCut <- which.max(sumCoef >= minscreen)
    whichVariable <- (as.matrix(fit.glmnet$glmnet.fit$beta)[,newCut] != 0)
  }
  return(whichVariable)
}

# survscreen.glmnet(time, event, X = data.frame(W), obsWeights=rep(1,length(time)), id=NULL)

#' Wrapper function for marginal Cox regression screening algorithm
#'
#' This screening algorithm uses marginal \code{\link[survival]{coxph}} regressions to select covariates that have significant marginal relationships with the event.
#'
#' @param time Observed follow-up time; i.e. minimum of the event and censoring times.
#' @param event Observed event indicator; i.e, whether the follow-up time corresponds to an event or censoring.
#' @param X Training covariate data.frame.
#' @param obsWeights Observation weights.
#' @param alpha Penalty exponent for \code{glmnet}. Defaults to 1 (lasso penalty).
#' @param minscreen Minimum number of covariates to return. Defaults to 2.
#' @param min.p Threshold p-value used to decide if a covariate is included. Defaults to 0.1
#' @param ... Additional ignored arguments.
#' @details A univariate Cox regression is run for each covariate; those with p-values less than \code{min.p} are included.
#' @return Logical vector of the same length as the number of columns of \code{X} indicating which variables were included.

survscreen.marg <- function(time, event, X, obsWeights, minscreen = 2, min.p = 0.1, ...) {
  pvals <- apply(X, 2, function(col) {
    est <- survival::coxph(survival::Surv(time, event) ~ ., data =  as.data.frame(cbind(time=time, event=event, col)),
                           weights = obsWeights)
    summary(est)$waldtest['pvalue']
  })
  whichVariable <- pvals <= .1
  if(sum(whichVariable) < minscreen) {
    whichVariable <- rep(FALSE, ncol(X))
    whichVariable[order(pvals)[1:minscreen]] <- TRUE
  }
  return(whichVariable)
}



#survscreen.marg.10(time, 1-event, X = data.frame(W), obsWeights=rep(1,length(time)), id=NULL)

All <- function(X, ...) {
  rep(TRUE, ncol(X))
}
