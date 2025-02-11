\name{plot.cv.grpnet}
\alias{plot.cv.grpnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Cross-Validation Curve for cv.grpnet Fits
}
\description{
Plots the mean cross-validation error, along with lower and upper standard deviation curves, as a function of \code{log(lambda)}.
}
\usage{
\method{plot}{cv.grpnet}(x, sign.lambda = 1, nzero = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class "cv.grpnet"
}
  \item{sign.lambda}{
Default plots \code{log(lambda)} on the x-axis. Set to -1 to plot \code{-1*log(lambda)} on the x-axis instead.
}
  \item{nzero}{
Should the number of non-zero groups be printed on the top of the x-axis?
}
  \item{...}{
Additional arguments passed to the \code{\link{plot}} function.
}
}
\details{
Produces cross-validation plot only (i.e., nothing is returned).
}
\value{
No return value (produces a plot)
}
\references{
Friedman, J.,  Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. \emph{Journal of Statistical Software, 33}(1), 1-22. \doi{10.18637/jss.v033.i01}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Syntax and functionality were modeled after the \code{plot.cv.glmnet} function in the \bold{glmnet} package (Friedman, Hastie, & Tibshirani, 2010).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cv.grpnet}} for k-fold cross-validation of \code{lambda}

\code{\link{plot.grpnet}} for plotting the regularization path
}
\examples{
# see 'cv.grpnet' for plotting examples
?cv.grpnet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ graphs }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
