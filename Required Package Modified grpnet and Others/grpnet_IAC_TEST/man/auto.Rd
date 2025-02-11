\name{auto}
\alias{auto}
\docType{data}
\title{
Auto MPG Data Set
}
\description{
Miles per gallon and other characteristics of vehicles from the 1970s-1980s. A version of this dataset was used as the 1983 American Statistical Association Exposition dataset.
}
\usage{data("auto")}
\format{
  A data frame with 392 observations on the following 9 variables.
  \describe{
    \item{\code{mpg}}{miles per gallon (numeric vector)}
    \item{\code{cylinders}}{number of cylinders: 3,4,5,6,8 (ordered factor)}
    \item{\code{displacement}}{engine displacement in cubic inches (numeric vector)}
    \item{\code{horsepower}}{engine horsepower (integer vector)}
    \item{\code{weight}}{vehicle weight in of lbs. (integer vector)}
    \item{\code{acceleration}}{0-60 mph time in sec. (numeric vector)}
    \item{\code{model.year}}{ranging from 1970 to 1982 (integer vector)}
    \item{\code{origin}}{region of origin: American, European, Japanese (factor vector)}
  }
}
\details{
This is a modified version of the "Auto MPG Data Set" on the UCI Machine Learning Repository, which is a modified version of the "cars" dataset on StatLib. 

Compared to the version of the dataset in UCI's MLR, this version of the dataset has removed (i) the 6 rows with missing \code{horsepower} scores, and (ii) the last column giving the name of each vehicle (\code{car.name}).
}
\source{
The dataset was originally collected by Ernesto Ramos and David Donoho.

StatLib---Datasets Archive at Carnegie Mellon University
http://lib.stat.cmu.edu/datasets/cars.data

Machine Learning Repository at University of California Irvine
https://archive.ics.uci.edu/ml/datasets/Auto+MPG
}
\examples{
# load data
data(auto)

# display structure
str(auto)

# display header
head(auto)

# see 'cv.grpnet' for cross-validation examples
?cv.grpnet

# see 'grpnet' for fitting examples
?grpnet
}
\keyword{datasets}
