\name{grpnet}
\alias{grpnet}
\alias{grpnet.default}
\alias{grpnet.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a Group Elastic Net Regularized GLM
}
\description{
Fits generalized linear models with a group elastic net penalty. Predictor groups can be manually input (default S3 method) or inferred from the model (S3 "formula" method). The regularization path is computed at a data-generated (default) or user-provided sequence of lambda values. }
\usage{
grpnet(x, \dots)

\method{grpnet}{default}(x, 
       y, 
       group, 
       family = c("gaussian", "binomial", "multinomial", "poisson", 
                  "negative.binomial", "Gamma", "inverse.gaussian"),
       weights = NULL, 
       offset = NULL, 
       alpha = 1, 
       nlambda = 100,
       lambda.min.ratio = ifelse(nobs < nvars, 0.05, 0.0001), 
       lambda = NULL, 
       penalty.factor = NULL,
       penalty = c("LASSO", "MCP", "SCAD"),
       gamma = ifelse(penalty == "MCP", 3, 4),
       theta = 1,
       standardize = TRUE, 
       intercept = TRUE, 
       thresh = 1e-04, 
       maxit = 1e05, 
       ...)
      
\method{grpnet}{formula}(formula,
       data, 
       use.rk = TRUE,
       family = c("gaussian", "binomial", "multinomial", "poisson", 
                  "negative.binomial", "Gamma", "inverse.gaussian"),
       weights = NULL,
       offset = NULL,
       alpha = 1,
       nlambda = 100,
       lambda.min.ratio = ifelse(nobs < nvars, 0.05, 0.0001),
       lambda = NULL,
       penalty.factor = NULL,
       penalty = c("LASSO", "MCP", "SCAD"),
       gamma = ifelse(penalty == "MCP", 3, 4),
       theta = 1,
       standardize = TRUE,
       thresh = 1e-04,
       maxit = 1e05,
       ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Model (design) matrix of dimension \code{nobs} by \code{nvars} (\eqn{n \times p}).
}
  \item{y}{
Response vector of length \eqn{n}. Matrix inputs are allowed for binomial and multinomial families (see "Binomial and multinomial" section).
}
  \item{group}{
Group label vector (factor, character, or integer) of length \eqn{p}. Predictors with the same label are grouped together for regularization.
}
  \item{formula}{
Model formula: a symbolic description of the model to be fitted. Uses the same syntax as \code{\link{lm}} and \code{\link{glm}}. 
}
  \item{data}{
Optional data frame containing the variables referenced in \code{formula}.
}
  \item{use.rk}{
If \code{TRUE} (default), the \code{\link{rk.model.matrix}} function is used to build the model matrix. Otherwise, the \code{\link{model.matrix}} function is used to build the model matrix. Additional arguments to the \code{\link{rk.model.matrix}} function can be passed via the \code{...} argument.
}
  \item{family}{
Character specifying the assumed distribution for the response variable. Partial matching is allowed. Options include \code{"gaussian"} (real-valued response), \code{"binomial"} (binary response), \code{"multinomial"} (multi-class response), \code{"poisson"} (count response), \code{"negative.binomial"} (count response), \code{"Gamma"} (positive real-valued), or \code{"inverse.gaussian"} (positive real-valued). 
}
  \item{weights}{
Optional vector of length \eqn{n} with non-negative weights to use for weighted (penalized) likelihood estimation. Defaults to a vector of ones.
}
  \item{offset}{
Optional vector of length \eqn{n} with an a priori known term to be included in the model's linear predictor. Defaults to a vector of zeros.  
}
  \item{alpha}{
Regularization hyperparameter satisfying \eqn{0 \leq \alpha \leq 1} that gives the balance between the group L1 (lasso) and group L2 (ridge) penalty. Setting \eqn{\alpha = 1} uses a group lasso penalty, setting \eqn{\alpha = 0} uses a group ridge penalty, and setting \eqn{0 < \alpha < 1} uses a group elastic net group penalty.  
}
  \item{nlambda}{
Number of \eqn{\lambda} values to use in the regularization path. Ignored if \code{lambda} is provided.
}
  \item{lambda.min.ratio}{
The proportion \eqn{0 < \pi < 1} that defines the minimum regularization parameter \eqn{\lambda_{\mathrm{min}}} as a fraction of the maximum regularization parameter \eqn{\lambda_{\mathrm{max}}} via the relationship \eqn{\lambda_{\mathrm{min}} = \pi \lambda_{\mathrm{max}}}. Ignored if \code{lambda} is provided. Note that \eqn{\lambda_{\mathrm{max}}} is defined such that all penalized effects are shrunk to zero.
}
  \item{lambda}{
Optional vector of user-supplied regularization parameter values. 
}
  \item{penalty.factor}{
Default S3 method: vector of length \eqn{K} giving the non-negative penalty weight for each predictor group. The order of the weights should correspond to the order of \code{levels(as.factor(group))}. Defaults to \eqn{\sqrt{p_k}} for all \eqn{k = 1,\ldots,K}, where \eqn{p_k} is the number of coefficients in the \eqn{k}-th group. If \code{penalty.factor[k] = 0}, then the \eqn{k}-th group is unpenalized, and the corresponding term is always included in the model. 

S3 "formula" method: named list giving the non-negative penalty weight for terms specified in the formula. Incomplete lists are allowed. Any term that is specified in \code{formula} but not in \code{penalty.factor} will be assigned the default penalty weight of \eqn{\sqrt{p_k}}. If \code{penalty.factor$z = 0}, then the variable \code{z} is unpenalized and always included in the model.
}
  \item{penalty}{
Character specifying which (group) penalty to use: LASSO , MCP, or SCAD.
}
  \item{gamma}{
Penalty hyperparameter that satisfies \eqn{\gamma > 1} for MCP and \eqn{\gamma > 2} for SCAD. Ignored for LASSO penalty.
} 
  \item{theta}{
Additional ("size") parameter for negative binomial responses, where the variance function is defined as \eqn{V(\mu) = \mu + \mu^2/ \theta}
}
  \item{standardize}{
Logical indicating whether the predictors should be groupwise standardized. If \code{TRUE} (default), each column of \code{x} is mean-centered and each predictor group is scaled to have a mean-square of one before fitting the model. Regardless of whether standardization is used, the coefficients are always returned on the original data scale.
}
  \item{intercept}{
Logical indicating whether an intercept term should be included in the model. Note that the intercept is always unpenalized.   
}
  \item{thresh}{
Convergence threshold (tolerance). The algorithm is determined to have converged once the maximum relative change in the coefficients is below this threshold. See "Convergence" section.
}
  \item{maxit}{
Maximum number of iterations to allow.
}
  \item{\dots}{
Additional arguments used by the default or formula method.
}
}
\details{
Consider a generalized linear model of the form
\deqn{
g(\mu) = \mathbf{X}^\top \boldsymbol\beta
}
where \eqn{\mu = E(Y | \mathbf{X})} is the conditional expectation of the response \eqn{Y} given the predictor vector \eqn{\mathbf{X}}, the function \eqn{g(\cdot)} is a user-specified (invertible) link function, and \eqn{\boldsymbol\beta} are the unknown regression coefficients. Furthermore, suppose that the predictors are grouped, such as
\deqn{
\mathbf{X}^\top \boldsymbol\beta = \sum_{k=1}^K \mathbf{X}_k^\top \boldsymbol\beta_k
}
where \eqn{\mathbf{X} = (\mathbf{X}_1, \ldots, \mathbf{X}_K)} is the grouped predictor vector, and \eqn{\boldsymbol\beta = (\boldsymbol\beta_1, \ldots, \boldsymbol\beta_K)} is the grouped coefficient vector. 

Given \eqn{n} observations, this function finds the \eqn{\boldsymbol\beta} that minimizes
\deqn{
L(\boldsymbol\beta | \mathbf{D}) + \lambda P_\alpha(\boldsymbol\beta)
}
where \eqn{L(\boldsymbol\beta | \mathbf{D})} is the loss function with \eqn{\mathbf{D} = \{\mathbf{y}, \mathbf{X}\}} denoting the observed data, \eqn{P_\alpha(\boldsymbol\beta)} is the group elastic net penalty, and \eqn{\lambda \geq 0} is the regularization parameter. 

The loss function has the form 
\deqn{
L(\boldsymbol\beta | \mathbf{D}) = \frac{1}{n} \sum_{i=1}^n w_i \ell_i(\boldsymbol\beta | \mathbf{D}_i)
}
where \eqn{w_i > 0} are the user-supplied \code{weights}, and \eqn{\ell_i(\boldsymbol\beta | \mathbf{D}_i)} is the \eqn{i}-th observation's contribution to the loss function. Note that \eqn{\ell(\cdot) = -\log(f_Y(\cdot))} denotes the negative log-likelihood function for the given \code{family}.

The group elastic net penalty function has the form
\deqn{
P_\alpha(\boldsymbol\beta) = \alpha P_1(\boldsymbol\beta) + (1 - \alpha) P_2(\boldsymbol\beta)
}
where \eqn{\alpha \in [0,1]} is the user-specified \code{alpha} value, 
\deqn{
P_1(\boldsymbol\beta) = \sum_{k=1}^K \omega_k \| \boldsymbol\beta_k \|
}
is the group lasso penalty with \eqn{\omega_k \geq 0} denoting the \eqn{k}-th group's \code{penalty.factor}, and 
\deqn{
P_2(\boldsymbol\beta) = \frac{1}{2} \sum_{k=1}^K \omega_k \| \boldsymbol\beta_k \|^2
}
is the group ridge penalty. Note that \eqn{\| \boldsymbol\beta_k \|^2 = \boldsymbol\beta_k^\top \boldsymbol\beta_k} denotes the squared Euclidean norm. When \code{penalty \%in\% c("MCP", "SCAD")}, the group L1 penalty \eqn{P_1(\boldsymbol\beta)} is replaced by the group MCP or group SCAD penalty.

}
\value{
An object of class \code{"grpnet"} with the following elements:
\item{call }{matched call}
\item{a0 }{intercept sequence of length \code{nlambda}}
\item{beta }{coefficient matrix of dimension \code{nvars} by \code{nlambda}}
\item{alpha }{balance between the group L1 (lasso) and group L2 (ridge) penalty}
\item{lambda}{sequence of regularization parameter values}
\item{family }{exponential family defining the loss function}
\item{dev.ratio }{proportion of (null) deviance explained for each \code{lambda} (= 1 - dev / nulldev)}
\item{nulldev }{null deviance for each \code{lambda}}
\item{df }{effective degrees of freedom for each \code{lambda}}
\item{nzgrp }{number of non-zero groups for each \code{lambda}}
\item{nzcoef }{number of non-zero coefficients for each \code{lambda}}
\item{ylev }{levels of response variable (only for binomial and multinomial families)}
\item{nobs }{number of observations}
\item{group}{group label vector}
\item{ngroups }{number of groups \eqn{K}}
\item{npasses }{number of iterations for each \code{lambda}}
\item{offset }{logical indicating if an offset was included}
\item{args }{list of input argument values}
\item{formula }{input formula (possibly after expansion)}
\item{term.labels }{terms that appear in formula (if applicable)}
\item{rk.args }{arguments for rk.model.matrix function (if applicable)}
}
\references{
Breheny, P., & Huang, J. (2015). Group descent algorithms for nonconvex penalized linear and logistic regression models with grouped predictors. \emph{Statistics and Computing, 25}, 173-187. \doi{10.1007/s11222-013-9424-2}

Friedman, J.,  Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. \emph{Journal of Statistical Software, 33}(1), 1-22. \doi{10.18637/jss.v033.i01}

Yang, Y., & Zou, H. (2015). A fast unified algorithm for solving group-lasso penalize learning problems. \emph{Statistics and Computing, 25}, 1129-1141. \doi{10.1007/s11222-014-9498-5}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The syntax of (the default S3 method for) this function closely mimics that of the \code{glmnet} function in the \bold{glmnet} package (Friedman, Hastie, & Tibshirani, 2010). 
}


\section{S3 "formula" method}{
\bold{Important}: When using the S3 "formula" method, the S3 "predict" method forms the model matrix for the predictions by applying the model formula to the new data. As a result, to ensure that the corresponding S3 "predict" method works correctly, some formulaic features should be avoided. 

\bold{Polynomials}: When including polynomial terms, the \code{\link{poly}} function should be used with option \code{raw = TRUE}. Default use of the \code{\link{poly}} function (with \code{raw = FALSE}) will work for fitting the model, but will result in invalid predictions for new data. Polynomials can also be included via the \code{\link{I}} function, but this isn't recommended because the polynomials terms wouldn't be grouped together, i.e., the terms \code{x} and \code{I(x^2)} would be treated as two separate groups of size one instead of a single group of size two.

\bold{Splines}: B-splines (and other spline bases) can be included via the S3 "formula" method. However, to ensure reasonable predictions for new data, it is necessary to specify the knots directly. For example, if \code{x} is a vector with entries between zero and one, the code \code{bs(x, df = 5)} will *not* produce valid predictions for new data, but the code \code{bs(x, knots = c(0.25, 0.5, 0.75), Boundary.knots = c(0, 1))} will work as intended. Instead of attempting to integrate a call to \code{bs()} or \code{rk()} into the model formula, it is recommended that splines be included via the \code{use.rk = TRUE} argument.
}

\section{Family argument and link functions}{
Unlike the \code{\link{glm}} function, the \code{family} argument of the \code{\link{grpnet}} function \cr
* should be a character vector (not a \code{\link{family}} object) \cr
* does not allow for specification of a link function 

Currently, there is only one available link function for each \code{family}: \cr
* gaussian (identity): \eqn{\mu = \mathbf{X}^\top \boldsymbol\beta} \cr
* binomial (logit): \eqn{\log(\frac{\pi}{1 - \pi}) = \mathbf{X}^\top \boldsymbol\beta} \cr
* multinomial (symmetric): \eqn{\pi_\ell = \frac{\exp(\mathbf{X}^\top \boldsymbol\beta_\ell)}{\sum_{l = 1}^m \exp(\mathbf{X}^\top \boldsymbol\beta_l)}} \cr
* poisson (log): \eqn{\log(\mu) = \mathbf{X}^\top \boldsymbol\beta} \cr
* negative.binomial (log): \eqn{\log(\mu) = \mathbf{X}^\top \boldsymbol\beta} \cr
* Gamma (log): \eqn{\log(\mu) = \mathbf{X}^\top \boldsymbol\beta} \cr
* inverse.gaussian (log): \eqn{\log(\mu) = \mathbf{X}^\top \boldsymbol\beta} \cr
}

\section{Binomial and multinomial}{
For \code{"binomial"} responses, three different possibilities exist for the input response: \cr
1. vector coercible into a factor with two levels \cr
2. matrix with two columns (# successes, # failures) \cr
3. numeric vector with entries between 0 and 1 \cr
In this case, the \code{weights} argument should be used specify the total number of trials.

For \code{"multinomial"} responses, two different possibilities exist for the input reponse: \cr
1. vector coercible into a factor with more than two levels \cr 
2. matrix of integers (counts) for each category level \cr
}

\section{Convergence}{
The algorithm is determined to have converged once

\eqn{\max_j \frac{| \beta_j - \beta_j^{\mathrm{old}} |}{1 + |\beta_j^{\mathrm{old}}|}  < \epsilon }

where \eqn{j \in \{1,\ldots,p\}} and \eqn{\epsilon} is the \code{thresh} argument. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.grpnet}} for plotting the regularization path

\code{\link{predict.grpnet}} for predicting from \code{\link{grpnet}} objects

\code{\link{cv.grpnet}} for k-fold cross-validation of \code{lambda}
}
\examples{
######***######   family = "gaussian"   ######***######

# load data
data(auto)

# fit model (formula method, response = mpg)
mod <- grpnet(mpg ~ ., data = auto)

# print regularization path info
mod

# plot coefficient paths
plot(mod)


\donttest{
######***######   family = "binomial"   ######***######

# load data
data(auto)

# define response (1 = American, 0 = other)
y <- ifelse(auto$origin == "American", 1, 0)

# define predictors
x <- rk.model.matrix(~ 0 + ., data = auto[,1:7])

# define group
g <- attr(x, "assign")

# fit model (default method, response = y)
mod <- grpnet(x, y, g, family = "binomial")

# print regularization path info
mod

# plot coefficient paths
plot(mod)



######***######   family = "poisson"   ######***######

# load data
data(auto)

# fit model (formula method, response = horsepower)
mod <- grpnet(horsepower ~ ., data = auto, family = "poisson")

# print regularization path info
mod

# plot coefficient paths
plot(mod)



######***######   family = "negative.binomial"   ######***######

# load data
data(auto)

# fit model (formula method, response = horsepower)
mod <- grpnet(horsepower ~ ., data = auto, family = "negative.binomial", theta = 100)

# print regularization path info
mod

# plot coefficient paths
plot(mod)



######***######   family = "multinomial"   ######***######

# load data
data(auto)

# fit model (formula method, response = origin)
mod <- grpnet(origin ~ ., data = auto, family = "multinomial")

# print regularization path info
mod

# plot coefficient paths
plot(mod)



######***######   family = "Gamma"   ######***######

# load data
data(auto)

# fit model (formula method, response = mpg)
mod <- grpnet(mpg ~ ., data = auto, family = "Gamma")

# print regularization path info
mod

# plot coefficient paths
plot(mod)



######***######   family = "inverse.gaussian"   ######***######

# load data
data(auto)

# fit model (formula method, response = mpg)
mod <- grpnet(mpg ~ ., data = auto, family = "inverse.gaussian")

# print regularization path info
mod

# plot coefficient paths
plot(mod)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ regression }
\keyword{ smooth }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
