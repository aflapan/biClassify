% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CompressedQDA.R
\name{QDA}
\alias{QDA}
\title{Quadratic Discriminant Analysis (QDA)}
\usage{
QDA(TrainData, TrainCat, TestData, Method = "Full", Mode = "Automatic",
  m1 = NULL, m2 = NULL, m = NULL, s = NULL, gamma = 1e-05)
}
\arguments{
\item{TrainData}{A (n x p) numeric matrix without missing values consisting of n training samples each with p features.}

\item{TrainCat}{A vector of length n consisting of group labels of the n training samples in \code{TrainData}. Must consist of 1s and 2s.}

\item{TestData}{A (m x p) numeric matrix without missing values consisting of m training samples each with p features. The number of features must equal the number of features in \code{TrainData}.}

\item{Method}{A string of characters which determinds which version of QDA to use. Must be either "Full", "Compressed", or "Subsampled".}

\item{Mode}{A string of characters which determines how the reduced sample paramters will be inputted for each method. Must be either "Research", "Interactive", or "Automatic". Default is "Automatic".}

\item{m1}{The number of class 1 compressed samples to be generated. Must be a positive integer.}

\item{m2}{The number of class 2 compressed samples to be generated. Must be a positive integer.}

\item{m}{The number of total compressed samples to be generated. Must be a positive integer.}

\item{s}{The sparsity level used in compression. Must satify 0 < s < 1.}

\item{gamma}{A numeric value for the stabilization amount gamma * I added to the covariance matrixed used in the LDA decision rule. Default amount is 1E-5. Cannot be negative.}
}
\value{
\item{Predictions}{(m x 1) Vector of predicted class labels for the data points in \code{TestData}.}
}
\description{
A wrapper function for the various QDA implementations available in this package.

Generates class predictions for \code{TestData}.
}
\details{
Function which handles all implementations of LDA.
}
\examples{
TrainData <- QDA_Data$TrainData
TrainCat <- QDA_Data$TrainCat
TestData <- QDA_Data$TestData
plot(TrainData[,2]~TrainData[,1], col = c("blue","orange")[as.factor(TrainCat)])

#----- Full QDA -------
QDA(TrainData = TrainData,
    TrainCat = TrainCat,
    TestData = TestData,
    Method = "Full",
    gamma = 1E-5)
  
#----- Compressed QDA -------  
 m1 <- 700
 m2 <- 300
 s <- 0.01
 QDA(TrainData = TrainData,
     TrainCat = TrainCat,
     TestData = TestData,
     Method = "Compressed",
     Mode = "Research",
     m1 = m1,
     m2 = m2,
     s = s,
     gamma = 1E-5)
     
 QDA(TrainData = TrainData,
     TrainCat = TrainCat,
     TestData = TestData,
     Method = "Compressed",
     Mode = "Automatic",
     gamma = 1E-5)
 
 #----- Sub-sampled QDA ------
 m1 <- 700
 m2 <- 300
 QDA(TrainData = TrainData,
     TrainCat = TrainCat,
     TestData = TestData,
     Method = "Subsampled",
     Mode = "Research",
     m1 = m1,
     m2 = m2,
     gamma = 1E-5)
     
 QDA(TrainData = TrainData,
     TrainCat = TrainCat,
     TestData = TestData,
     Method = "Subsampled",
     Mode = "Automatic",
     gamma = 1E-5)
     
}
\references{
Lapanowski, Alexander F., and Gaynanova, Irina. ``Compressing large sample data for discriminant analysis'' arXiv preprint arXiv:2005.03858 (2020).
}
