\name{print}
\alias{print.coef.grpnet}
\alias{print.cv.grpnet}
\alias{print.grpnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
S3 'print' Methods for grpnet
}
\description{
Prints some basic information about the coefficients (for \code{\link{coef.grpnet}} objects), observed cross-validation error (for \code{\link{cv.grpnet}} objects), or the computed regularization path (for \code{\link{grpnet}} objects). 
}
\usage{
\method{print}{coef.grpnet}(x, ...)

\method{print}{cv.grpnet}(x, digits = max(3, getOption("digits") - 3), ...)

\method{print}{grpnet}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{\link{coef.grpnet}}, \code{\link{cv.grpnet}}, or \code{\link{grpnet}}
}
  \item{digits}{
the number of digits to print (must be a positive integer)  
}
  \item{...}{
additional arguments for \code{\link{print}} (currently ignored)  
}
}
\details{
For \code{\link{coef.grpnet}} objects, prints the non-zero coefficients and uses "." for coefficients shrunk to zero.

For \code{\link{cv.grpnet}} objects, prints the function \code{call}, the cross-validation \code{type.measure}, and a two-row table with information about the \code{min} and \code{1se} solutions.

For \code{\link{grpnet}} objects, prints a data frame with columns \cr
* nGrp: number of non-zero groups for each \code{lambda} \cr
* Df: effective degrees of freedom for each \code{lambda} \cr
* \%Dev: percentage of null deviance explained for each \code{lambda} \cr
* Lambda: the values of \code{lambda} \cr
}

\value{
No return value (produces a printout)
}

\references{
Friedman, J.,  Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. \emph{Journal of Statistical Software, 33}(1), 1-22. \doi{10.18637/jss.v033.i01}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Some syntax and functionality were modeled after the \code{print} functions in the \bold{glmnet} package (Friedman, Hastie, & Tibshirani, 2010).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{coef.grpnet}} for extracting coefficients

\code{\link{cv.grpnet}} for k-fold cross-validation of \code{lambda}

\code{\link{grpnet}} for fitting grpnet regularization paths
}
\examples{
# see 'coef.grpnet' for coefficient printing examples
?coef.grpnet

# see 'cv.grpnet' for cross-validation error printing examples
?cv.grpnet

# see 'grpnet' for regularization path printing examples
?grpnet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ print }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
