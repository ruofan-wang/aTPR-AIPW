\name{atpr_fit_methods}
\alias{atpr_fit_methods}
\alias{print.atpr_fit}
\alias{coef.atpr_fit}
\alias{vcov.atpr_fit}
\alias{confint.atpr_fit}
\alias{summary.atpr_fit}
\alias{print.summary.atpr_fit}
\title{Methods for fitted adjusted-TPR contrasts}
\usage{
\method{print}{atpr_fit}(x, digits = 4, ...)

\method{coef}{atpr_fit}(object, ...)

\method{vcov}{atpr_fit}(object, ...)

\method{confint}{atpr_fit}(object, parm = "delta", level = 0.95, ...)

\method{summary}{atpr_fit}(object, level = 0.95, ...)

\method{print}{summary.atpr_fit}(x, digits = 4, ...)
}
\arguments{
\item{x, object}{An \code{atpr_fit} object returned by \code{\link[=fit_atpr]{fit_atpr()}}.}

\item{digits}{Number of digits printed.}

\item{...}{Additional arguments; currently unused.}

\item{parm}{Unused; the object contains the single contrast delta.}

\item{level}{Confidence level for a working normal-approximation interval.}
}
\value{
\code{coef()} returns the named contrast; \code{vcov()} returns its estimated
variance; \code{confint()} returns a two-column matrix; and \code{summary()} returns a
list with estimates, diagnostics, and a working normal interval.
Print methods invisibly return their input.
}
\description{
Print, summarize, or extract a point estimate and its working
influence-function uncertainty.
}
