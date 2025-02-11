\name{rk.model.matrix}
\alias{rk.model.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct Design Matrices via Reproducing Kernels
}
\description{
Creates a design (or model) matrix using the \code{\link{rk}} function to expand variables via a reproducing kernel basis. 
}
\usage{
rk.model.matrix(object, data = environment(object), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a \code{\link{formula}} or \code{\link{terms}} object describing the fit model
}
  \item{data}{
a data frame containing the variables referenced in \code{object}
}
  \item{...}{
additional arguments passed to the \code{\link{rk}} function, e.g., \code{df}, \code{knots}, \code{m}, etc. Arguments must be passed as a named list, see Examples.
}
}
\details{
Designed to be a more flexible alternative to the \code{\link{model.matrix}} function. The \code{\link{rk}} function is used to construct a marginal basis for each variable that appears in the input \code{object}. Tensor product interactions are formed by taking a \code{\link{row.kronecker}} product of marginal basis matrices. Interactions of any order are supported using standard formulaic conventions, see Note.
}
\value{
The design matrix corresponding to the input formula and data, which has the following attributes:

\item{assign}{an integer vector with an entry for each column in the matrix giving the term in the formula which gave rise to the column}
\item{term.labels}{a character vector containing the labels for each of the terms in the model}
\item{knots}{a named list giving the knots used for each variable in the formula}
\item{m}{a named list giving the penalty order used for each variable in the formula}
\item{periodic}{a named list giving the periodicity used for each variable in the formula}
\item{xlev}{a named list giving the factor levels used for each variable in the formula}
}
\references{
Helwig, N. E. (2017). Regression with ordered predictors via ordinal smoothing splines. Frontiers in Applied Mathematics and Statistics, 3(15), 1-13. \doi{10.3389/fams.2017.00015}

Helwig, N. E. (2021). Spectrally sparse nonparametric regression via elastic net regularized smoothers. Journal of Computational and Graphical Statistics, 30(1), 182-191. \doi{10.1080/10618600.2020.1806855}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
For formulas of the form \code{y ~ x + z}, the constructed model matrix has the form \code{cbind(rk(x), rk(z))}, which simply concatenates the two marginal basis matrices. For formulas of the form \code{y ~ x : z}, the constructed model matrix has the form \code{row.kronecker(rk(x), rk(z))}, where \code{\link{row.kronecker}} denotes the row-wise kronecker product. The formula \code{y ~ x * z} is a shorthand for \code{y ~ x + z + x : z}, which concatenates the two previous results. Unless it is suppressed (using \code{0+}), the first column of the basis will be a column of ones named \code{(Intercept)}.
}

\seealso{
See \code{\link{rk}} for details on the reproducing kernel basis
}
\examples{
# load auto data
data(auto)

# additive effects
x <- rk.model.matrix(mpg ~ ., data = auto)
dim(x)                      # check dimensions
attr(x, "assign")           # check group assignments
attr(x, "term.labels")      # check term labels

# two-way interactions
x <- rk.model.matrix(mpg ~ . * ., data = auto)
dim(x)                      # check dimensions
attr(x, "assign")           # check group assignments
attr(x, "term.labels")      # check term labels

# specify df for horsepower, weight, and acceleration
# note: default df = 5 is used for displacement and model.year
df <- list(horsepower = 6, weight = 7, acceleration = 8)
x <- rk.model.matrix(mpg ~ ., data = auto, df = df)
sapply(attr(x, "knots"), length)   # check df

# specify knots for model.year
# note: default knots are selected for other variables
knots <- list(model.year = c(1970, 1974, 1978, 1982))
x <- rk.model.matrix(mpg ~ ., data = auto, knots = knots)
sapply(attr(x, "knots"), length)   # check df

}
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ smooth }% __ONLY ONE__ keyword per line
