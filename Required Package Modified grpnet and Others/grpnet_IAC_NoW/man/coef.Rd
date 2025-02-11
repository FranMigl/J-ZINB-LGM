\name{coef}
\alias{coef.cv.grpnet}
\alias{coef.grpnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract Coefficients for cv.grpnet and grpnet Fits
}
\description{
Obtain coefficients from a cross-validated group elastic net regularized GLM (cv.grpnet) or a group elastic net regularized GLM (grpnet) object.
}
\usage{
\method{coef}{cv.grpnet}(object, 
     s = c("lambda.min", "lambda.1se"),
     ...)
     
\method{coef}{grpnet}(object, 
     s = NULL,
     ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class "cv.grpnet" or "grpnet"
}
\item{s}{
Lambda value(s) at which predictions should be obtained. For "cv.grpnet" objects, default uses the lambda that minimizes the cross-validation loss function. For "grpnet" objects, default uses \code{s = object$lambda}. Interpolation is used for \code{s} values that are not included in \code{object$lambda}.  
}
\item{...}{
Additional arguments (ignored)
}
}
\details{
\emph{coef.cv.grpnet}: \cr
Returns the coefficients that are used by the \code{\link{predict.cv.grpnet}} function to form predictions from a fit \code{\link{cv.grpnet}} object. 

\emph{coef.grpnet}: \cr
Returns the coefficients that are used by the \code{\link{predict.grpnet}} function to form predictions from a fit \code{\link{grpnet}} object. 
}
\value{
For multinomial response variables, returns a list of length \code{length(object$ylev)}, where the \code{j}-th element is a matrix of dimension \code{c(ncoef, length(s))} giving the coefficients for \code{object$ylev[j]}.

For other response variables, returns a matrix of dimension \code{c(ncoef, length(s))}, where the \code{i}-th column gives the coefficients for \code{s[i]}.
}
\references{
Friedman, J.,  Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. \emph{Journal of Statistical Software, 33}(1), 1-22. \doi{10.18637/jss.v033.i01}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The syntax of these functions closely mimics that of the \code{coef.cv.glmnet} and \code{coef.glmnet} functions in the \bold{glmnet} package (Friedman, Hastie, & Tibshirani, 2010).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{print.coef.grpnet}} for printing \code{\link{coef.grpnet}} objects

\code{\link{predict.cv.grpnet}} for predicting from \code{\link{cv.grpnet}} objects

\code{\link{predict.grpnet}} for predicting from \code{\link{grpnet}} objects
}
\examples{
######***######   grpnet   ######***######

# load data
data(auto)

# fit model (formula method, response = mpg)
mod <- grpnet(mpg ~ ., data = auto)

# extract coefs for regularization path (output = 12 x 100 matrix)
coef(mod)

# extract coefs at 3 particular points (output = 12 x 3 matrix)
coef(mod, s = c(1.5, 1, 0.5))


######***######   cv.grpnet   ######***######

# load data
data(auto)

# 10-fold cv (formula method, response = mpg)
set.seed(1)
mod <- cv.grpnet(mpg ~ ., data = auto, alpha = 1)

# extract coefs for "min" solution (output = 12 x 1 matrix)
coef(mod)

# extract coefs for "1se" solution (output = 12 x 1 matrix)
coef(mod, s = "lambda.1se")

# extract coefs at 3 particular points (output = 12 x 3 matrix)
coef(mod, s = c(1.5, 1, 0.5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ regression }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
