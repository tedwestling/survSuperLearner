% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SL_wrappers.R
\name{survSL.rfsrc}
\alias{survSL.rfsrc}
\title{Wrapper function for Random Survival Forests prediction algorithm}
\usage{
survSL.rfsrc(time, event, X, newX, new.times, obsWeights, id, ...)
}
\arguments{
\item{time}{Observed follow-up time; i.e. minimum of the event and censoring times.}

\item{event}{Observed event indicator; i.e, whether the follow-up time corresponds to an event or censoring.}

\item{X}{Training covariate data.frame.}

\item{newX}{Test covariate data.frame to use for prediction. Should have the same variable names and structure as \code{X}.}

\item{new.times}{Times at which to obtain to obtain the predicted survivals.}

\item{obsWeights}{Observation weights.}

\item{id}{Currently ignored}

\item{...}{Additional arguments passed on to \code{\link[randomForestSRC]{rfsrc}}.}
}
\value{
\item{pred}{Matrix of predictions, with the same number of rows as \code{newX} and number of columns equal to the length of \code{new.times}. Rows index new observations, and columns index new times at which the survival was computed.}
\item{fit}{One-element list including \code{object}, the fitted \code{\link[randomForestSRC]{rfsrc}}. object.}
}
\description{
This prediciton algorithm uses the \code{\link[randomForestSRC]{rfsrc}} function from the \code{randomForestSRC} package to estimate a survival random forest.
}
\references{
Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. \emph{The Annals of Applied Statistics}, 2(3), 841-860.
}
