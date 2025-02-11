\name{predict.grpnet}
\alias{predict.grpnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict Method for grpnet Fits
}
\description{
Obtain predictions from a fit group elastic net regularized GLM (grpnet) object.
}
\usage{
\method{predict}{grpnet}(object, 
        newx,
        newdata,
        s = NULL,
        type = c("link", "response", "class", "terms", "importance", 
                 "coefficients", "nonzero", "groups", "ncoefs", "ngroups", "norm"),
        ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object of class "grpnet"
}
  \item{newx}{
Matrix of new \code{x} scores for prediction (default S3 method). Must have \eqn{p} columns arranged in the same order as the \code{x} matrix used to fit the model. Ignored for the last six types of predictions.
}
  \item{newdata}{
Data frame of new \code{data} scores for prediction (S3 "formula" method). Must contain all variables in the \code{formula} used to fit the model. Ignored for the last six types of predictions.
}
  \item{s}{
Lambda value(s) at which predictions should be obtained. Default uses \code{s = object$lambda}. Interpolation is used for \code{s} values that are not included in \code{object$lambda}.  
}
  \item{type}{
Type of prediction to return. "link" gives predictions on the link scale (\eqn{\eta}). "response" gives predictions on the mean scale (\eqn{\mu}). "class" gives predicted class labels (for "binomial" and "multinomial" families). "terms" gives the predictions for each term (group) in the model (\eqn{\eta_k}). "importance" gives the variable importance index for each term (group) in the model. "coefficients" returns the coefficients used for predictions. "nonzero" returns a list giving the indices of non-zero coefficients for each \code{s}. "groups" returns a list giving the labels of non-zero groups for each \code{s}. "ncoefs" returns the number of non-zero coefficients for each \code{s}. "ngroups" returns the number of non-zero groups for each \code{s}. "norm" returns the L2 norm of each group's coefficients for each \code{s}.
}
  \item{...}{
Additional arguments (ignored)
}
}
\details{
When \code{type == "link"}, the predictions for each \eqn{\lambda} have the form
\deqn{\boldsymbol\eta_\lambda = \mathbf{X}_{\mathrm{new}} \boldsymbol\beta_\lambda}
where \eqn{\mathbf{X}_{\mathrm{new}}} is the argument \code{newx} (or the design matrix created from \code{newdata} by applying \code{object$formula}) and \eqn{\boldsymbol\beta_\lambda} is the coefficient vector corresponding to \eqn{\lambda}.

When \code{type == "response"}, the predictions for each \eqn{\lambda} have the form
\deqn{\boldsymbol\mu_\lambda = g^{-1}(\boldsymbol\eta_\lambda) }
where \eqn{g^{-1}(\cdot)} is the inverse link function stored in \code{object$family$linkinv}.

When \code{type == "class"}, the predictions for each \eqn{\lambda} have the form
\deqn{\mathbf{y}_\lambda = \arg\max_l \boldsymbol\mu_\lambda(l) }
where \eqn{\boldsymbol\mu_\lambda(l)} gives the predicted probability that each observation belongs to the \eqn{l}-th category (for \eqn{l = 1,\ldots,m}) using the regularization parameter \eqn{\lambda}. 

When \code{type == "terms"}, the groupwise predictions for each \eqn{\lambda} have the form
\deqn{\boldsymbol\eta_{k\lambda} = \mathbf{X}_k^{\mathrm{(new)}} \boldsymbol\beta_{k\lambda}} where \eqn{\mathbf{X}_k^{\mathrm{(new)}}} is the portion of the argument \code{newx} (or the design matrix created from \code{newdata} by applying \code{object$formula}) that corresponds to the \eqn{k}-th term/group, and \eqn{\boldsymbol\beta_{k\lambda}} are the corresponding coefficients.

When \code{type == "importance"}, the variable importance indices are defined as 
\deqn{\pi_k = \left( \boldsymbol\eta_{k\lambda}^\top \mathbf{C} \boldsymbol\eta_{0\lambda} \right) \left( \boldsymbol\eta_{0\lambda}^\top \mathbf{C} \boldsymbol\eta_{0\lambda} \right)^{-1}}
where \eqn{\mathbf{C} = (\mathbf{I}_n - \frac{1}{n} \mathbf{1}_n \mathbf{1}_n^\top)} denotes a centering matrix, and \eqn{\boldsymbol\eta_{0\lambda} = \sum_{k=1}^K \boldsymbol\eta_{k\lambda}}. Note that \eqn{\sum_{k=1}^K \pi_k = 1}, but some \eqn{\pi_k} could be negative. When they are positive, \eqn{\pi_k} gives the approximate proportion of model (explained) variation that is attributed to the \eqn{k}-th term.
}
\value{
Depends on three factors... \cr
1. the exponential family distribution \cr
2. the length of the input \code{s} \cr
3. the \code{type} of prediction requested

For most response variables, the typical output will be...
\item{*}{a matrix of dimension \code{c(newnobs, length(s))} if \code{length(s) > 1}}
\item{*}{a vector of length \code{newnobs} if length(s) == 1}

For multinomial response variables, the typical output will be... 
\item{*}{an array of dimension \code{c(newnobs, length(object$ylev), length(s))} if \code{type \%in\% c("link", "response")}}
\item{*}{a matrix of dimension \code{c(newobs, length(s))} if \code{type == "class"}}

Note: if \code{type == "class"}, then the output will be the same class as \code{object$ylev}. Otherwise, the output will be real-valued (or integer for the counts).

If \code{type == "terms"} and \code{family != "multinomial"}, the output will be...
\item{*}{an array of dimension \code{c(newnobs, nterms, length(s))} if \code{length(s) > 1}}
\item{*}{a matrix of dimension \code{c(newnobs, nterms)} if \code{length(s) == 1}}

If \code{type == "terms"} and \code{family == "multinomial"}, the output will be a list of length \code{length(object$ylev)} where each element gives the terms for the corresponding response class.

If \code{type == "importance"} and \code{family != "multinomial"}, the output will be...
\item{*}{a matrix of dimension \code{c(nterms, length(s))} if \code{length(s) > 1}}
\item{*}{a vector of length \code{nterms} if \code{length(s) == 1}}

If \code{type == "importance"} and \code{family == "multinomial"}, the output will be a list of length \code{length(object$ylev)} where each element gives the importance for the corresponding response class. If \code{length(s) == 1}, the output will be simplified to matrix.

If \code{type == "coefficients"}, the output will be the same as that produced by \code{\link{coef.grpnet}}.

If \code{type == "nonzero"}, the output will be a list of length \code{length(s)} where each element is a vector of integers (indices).

If \code{type == "groups"}, the output will be a list of length \code{length(s)} where each element is a vector of characters (\code{term.labels}).

If \code{type \%in\% c("ncoefs", "ngroups")}, the output will be a vector of length \code{length(s)} where each element is an integer.

If \code{type == "norm"}, the output will be a matrix of dimension \code{c(K, length(s))}, where each cell gives the L2 norm for the corresponding group and smoothing parameter. Note that \code{K} denotes the number of groups.

}
\references{
Friedman, J.,  Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. \emph{Journal of Statistical Software, 33}(1), 1-22. \doi{10.18637/jss.v033.i01}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Some internal code (e.g., used for the interpolation) is borrowed from the \code{predict.glmnet} function in the \bold{glmnet} package (Friedman, Hastie, & Tibshirani, 2010).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{grpnet}} for fitting grpnet regularization paths

\code{\link{predict.cv.grpnet}} for predicting from \code{\link{cv.grpnet}} objects
}
\examples{
######***######   family = "gaussian"   ######***######

# load data
data(auto)

# fit model (formula method, response = mpg)
mod <- grpnet(mpg ~ ., data = auto)

# get fitted values for regularization path (output = 392 x 100 matrix)
fit.path <- predict(mod, newdata = auto)

# get fitted values at 3 particular points (output = 392 x 3 matrix)
fit.some <- predict(mod, newdata = auto, s = c(1.5, 1, 0.5))

# compare rmse for solutions
rmse.path <- sqrt(colMeans((auto$mpg - fit.path)^2))
rmse.some <- sqrt(colMeans((auto$mpg - fit.some)^2))
plot(log(mod$lambda), rmse.path, cex = 0.5)
points(log(c(1.5, 1, 0.5)), rmse.some, pch = 0, col = "red")


\donttest{
######***######   family = "binomial"   ######***######

# load data
data(auto)

# define response (1 = American, 0 = other)
y <- ifelse(auto$origin == "American", 1, 0)

# define predictors
x <- model.matrix(~ ., data = auto[,1:7])[,-1]

# define group (according to colnames(x))
g <- c(1, 2, 2, 2, 2, 3, 4, 5, 6, 7)

# fit model (default method, response = y)
mod <- grpnet(x, y, g, family = "binomial")

# get predicted classes for regularization path (output = 392 x 100 matrix)
fit.path <- predict(mod, newx = x, type = "class")

# get predicted classes at 3 particular points (output = 392 x 3 matrix)
fit.some <- predict(mod, newx = x, type = "class", s = c(.15, .1, .05))

# compare misclassification rate for solutions
miss.path <- 1 - colMeans(y == fit.path)
miss.some <- 1 - colMeans(y == fit.some)
plot(log(mod$lambda), miss.path, cex = 0.5)
points(log(c(.15, .1, .05)), miss.some, pch = 0, col = "red")



######***######   family = "poisson"   ######***######

# load data
data(auto)

# fit model (formula method, response = horsepower)
mod <- grpnet(horsepower ~ ., data = auto, family = "poisson")

# get fitted values for regularization path (output = 392 x 100 matrix)
fit.path <- predict(mod, newdata = auto, type = "response")

# get fitted values at 3 particular points (output = 392 x 3 matrix)
fit.some <- predict(mod, newdata = auto, type = "response", s = c(15, 10, 5))

# compare rmse for solutions
rmse.path <- sqrt(colMeans((auto$horsepower - fit.path)^2))
rmse.some <- sqrt(colMeans((auto$horsepower - fit.some)^2))
plot(log(mod$lambda), rmse.path, cex = 0.5)
points(log(c(15, 10, 5)), rmse.some, pch = 0, col = "red")



######***######   family = "negative.binomial"   ######***######

# load data
data(auto)

# fit model (formula method, response = horsepower)
mod <- grpnet(horsepower ~ ., data = auto, family = "negative.binomial", theta = 100)

# get fitted values for regularization path (output = 392 x 100 matrix)
fit.path <- predict(mod, newdata = auto, type = "response")

# get fitted values at 3 particular points (output = 392 x 3 matrix)
fit.some <- predict(mod, newdata = auto, type = "response", s = c(15, 10, 5))

# compare rmse for solutions
rmse.path <- sqrt(colMeans((auto$horsepower - fit.path)^2))
rmse.some <- sqrt(colMeans((auto$horsepower - fit.some)^2))
plot(log(mod$lambda), rmse.path, cex = 0.5)
points(log(c(15, 10, 5)), rmse.some, pch = 0, col = "red")



######***######   family = "multinomial"   ######***######

# load data
data(auto)

# fit model (formula method, response = origin)
mod <- grpnet(origin ~ ., data = auto, family = "multinomial")

# get predicted classes for regularization path (output = 392 x 100 matrix)
fit.path <- predict(mod, newdata = auto, type = "class")

# get predicted classes at 3 particular points (output = 392 x 3 matrix)
fit.some <- predict(mod, newdata = auto, type = "class", s = c(.1, .01, .001))

# compare misclassification rate for solutions
miss.path <- 1 - colMeans(auto$origin == fit.path)
miss.some <- 1 - colMeans(auto$origin == fit.some)
plot(log(mod$lambda), miss.path, cex = 0.5)
points(log(c(.1, .01, .001)), miss.some, pch = 0, col = "red")



######***######   family = "Gamma"   ######***######

# load data
data(auto)

# fit model (formula method, response = mpg)
mod <- grpnet(mpg ~ ., data = auto, family = "Gamma")

# get fitted values for regularization path (output = 392 x 100 matrix)
fit.path <- predict(mod, newdata = auto, type = "response")

# get fitted values at 3 particular points (output = 392 x 3 matrix)
fit.some <- predict(mod, newdata = auto, type = "response", s = c(0.2, 0.1, 0.05))

# compare rmse for solutions
rmse.path <- sqrt(colMeans((auto$mpg - fit.path)^2))
rmse.some <- sqrt(colMeans((auto$mpg - fit.some)^2))
plot(log(mod$lambda), rmse.path, cex = 0.5)
points(log(c(0.2, 0.1, 0.05)), rmse.some, pch = 0, col = "red")



######***######   family = "inverse.gaussian"   ######***######

# load data
data(auto)

# fit model (formula method, response = mpg)
mod <- grpnet(mpg ~ ., data = auto, family = "inverse.gaussian")

# get fitted values for regularization path (output = 392 x 100 matrix)
fit.path <- predict(mod, newdata = auto, type = "response")

# get fitted values at 3 particular points (output = 392 x 3 matrix)
fit.some <- predict(mod, newdata = auto, type = "response", s = c(0.2, 0.1, 0.05))

# compare rmse for solutions
rmse.path <- sqrt(colMeans((auto$mpg - fit.path)^2))
rmse.some <- sqrt(colMeans((auto$mpg - fit.some)^2))
plot(log(mod$lambda), rmse.path, cex = 0.5)
points(log(c(0.2, 0.1, 0.05)), rmse.some, pch = 0, col = "red")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ regression }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
