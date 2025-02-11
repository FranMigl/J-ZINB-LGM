\name{row.kronecker}
\alias{row.kronecker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Row-Wise Kronecker Product
}
\description{
Calculates the row-wise Kronecker product between two matrices with the same number of rows.
}
\usage{
row.kronecker(X, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
matrix of dimension \eqn{n \times p}
}
  \item{Y}{
matrix of dimension \eqn{n \times q}
}
}
\details{
Given \code{X} of dimension \code{c(n, p)} and \code{Y} of dimension \code{c(n, q)}, this function returns

\code{cbind(x[,1] * Y, x[,2] * Y, ..., x[,p] * Y)}

which is a matrix of dimension \code{c(n, p*q)}
}
\value{
Matrix of dimension \eqn{n \times pq} where each row contains the Kronecker product between the corresponding rows of \code{X} and \code{Y}.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Used by the \code{\link{rk.model.matrix}} to construct basis functions for interaction terms

See \code{\link{kronecker}} for the regular kronecker product
}
\examples{
X <- matrix(c(1, 1, 2, 2), nrow = 2, ncol = 2)
Y <- matrix(1:6, nrow = 2, ncol = 3)
row.kronecker(X, Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ array }
\keyword{ algebra }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
