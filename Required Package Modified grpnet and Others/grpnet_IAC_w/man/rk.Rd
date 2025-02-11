\name{rk}
\alias{rk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reproducing Kernel Basis
}
\description{
Generate a reproducing kernel basis matrix for a nominal, ordinal, or polynomial smoothing spline.
}
\usage{
rk(x, df = NULL, knots = NULL, m = NULL, intercept = FALSE, 
   Boundary.knots = NULL, warn.outside = TRUE, 
   periodic = FALSE, xlev = levels(x))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  the predictor vector of length \code{n}. Can be a factor, integer, or numeric, see Note.
  }
  \item{df}{
  the degrees of freedom, i.e., number of knots to place at quantiles of \code{x}. Defaults to 5 but ignored if \code{knots} are provided.
  }
  \item{knots}{
  the breakpoints (knots) defining the spline. If \code{knots} are provided, the \code{df} is defined as \code{length(unique(c(knots, Boundary.knots)))}.
  }
  \item{m}{
  the derivative penalty order: 0 = ordinal spline, 1 = linear spline, 2 = cubic spline, 3 = quintic spline
  }
  \item{intercept}{
  should an intercept be included in the basis?  
  }
  \item{Boundary.knots}{
  the boundary points for spline basis. Defaults to \code{range(x)}.
  }
  \item{warn.outside}{
  if \code{TRUE}, a warning is provided when \code{x} values are outside of the \code{Boundary.knots}
  }
  \item{periodic}{
  should the spline basis functions be constrained to be periodic with respect to the \code{Boundary.knots}?
  }
  \item{xlev}{
  levels of \code{x} (only applicable if \code{x} is a \code{\link{factor}})
  }
}
\details{
Given a vector of function realizations \eqn{f}, suppose that \eqn{f = X \beta}, where \eqn{X} is the (unregularized) spline basis and \eqn{\beta} is the coefficient vector. Let \eqn{Q} denote the postive semi-definite penalty matrix, such that \eqn{\beta^\top Q \beta} defines the roughness penalty for the spline. See Helwig (2017) for the form of \eqn{X} and \eqn{Q} for the various types of splines.

Consider the spectral parameterization of the form \eqn{f = Z \alpha} where
\deqn{Z = X Q^{-1/2}}
is the regularized spline basis (that is returned by this function), and \eqn{\alpha = Q^{1/2} \beta} are the reparameterized coefficients. Note that \eqn{X \beta = Z \alpha} and \eqn{\beta^\top Q \beta = \alpha^\top \alpha}, so the spectral parameterization absorbs the penalty into the coefficients (see Helwig, 2021).

Syntax of this function is designed to mimic the syntax of the \code{\link[splines]{bs}} function.
}
\value{
Returns a basis function matrix of dimension \code{n} by \code{df} (plus 1 if an \code{intercept} is included) with the following attributes:
\item{df }{degrees of freedom}
\item{knots }{knots for spline basis}
\item{m }{derivative penalty order}
\item{intercept }{was an intercept included?}
\item{Boundary.knots }{boundary points of \code{x}}
\item{periodic }{is the basis periodic?}
\item{xlev }{factor levels (if applicable)}
}
\references{
Helwig, N. E. (2017). Regression with ordered predictors via ordinal smoothing splines. Frontiers in Applied Mathematics and Statistics, 3(15), 1-13. \doi{10.3389/fams.2017.00015}

Helwig, N. E. (2021). Spectrally sparse nonparametric regression via elastic net regularized smoothers. Journal of Computational and Graphical Statistics, 30(1), 182-191. \doi{10.1080/10618600.2020.1806855}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The (default) type of spline basis depends on the \code{\link{class}} of the input \code{x} object:

* If \code{x} is an unordered factor, then a nominal spline basis is used

* If \code{x} is an ordered factor (and \code{m = NULL}), then an ordinal spline basis is used

* If \code{x} is an integer or numeric (and \code{m = NULL}), then a cubic spline basis is used

Note that you can override the default behavior by specifying the \code{m} argument.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

######***######   NOMINAL SPLINE BASIS   ######***######

x <- as.factor(LETTERS[1:5])
basis <- rk(x)
plot(1:5, basis[,1], t = "l", ylim = extendrange(basis))
for(j in 2:ncol(basis)){
  lines(1:5, basis[,j], col = j)
}


######***######   ORDINAL SPLINE BASIS   ######***######

x <- as.ordered(LETTERS[1:5])
basis <- rk(x)
plot(1:5, basis[,1], t = "l", ylim = extendrange(basis))
for(j in 2:ncol(basis)){
  lines(1:5, basis[,j], col = j)
}


######***######   LINEAR SPLINE BASIS   ######***######

x <- seq(0, 1, length.out = 101)
basis <- rk(x, m = 1)
plot(x, basis[,1], t = "l", ylim = extendrange(basis))
for(j in 2:ncol(basis)){
  lines(x, basis[,j], col = j)
}


######***######   CUBIC SPLINE BASIS   ######***######

x <- seq(0, 1, length.out = 101)
basis <- rk(x)
basis <- scale(basis)  # for visualization only!
plot(x, basis[,1], t = "l", ylim = extendrange(basis))
for(j in 2:ncol(basis)){
  lines(x, basis[,j], col = j)
}


######***######   QUINTIC SPLINE BASIS   ######***######

x <- seq(0, 1, length.out = 101)
basis <- rk(x, m = 3)
basis <- scale(basis)  # for visualization only!
plot(x, basis[,1], t = "l", ylim = extendrange(basis))
for(j in 2:ncol(basis)){
  lines(x, basis[,j], col = j)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ smooth }% __ONLY ONE__ keyword per line
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
