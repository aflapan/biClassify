
## Forms Categorical Response Matrix and Diagonal group proportion matrix. Each row
## consists of 0s and a single 1.  The 1 is in the column corresponding to the group
## data point i belongs to.
## Input: length n vector of categories. Categories labeled as 1, 2, 3, ...
## Output:  (n x G) categorical response matrix Y. 
IndicatMat<- function(TrainCat){
  ### Category is a vector of length n. The ith entry lists which group the ith data
  ### point belongs to.
  E <- matrix(0 , nrow = length(TrainCat), ncol = 2)
  E[TrainCat == 1 , 1] <- 1
  E[TrainCat == 2, 2] <- 1
  Dpi <- crossprod(E, E) / nrow(E) ### This forms the group proportion matrix
  return( list(Categorical = E, Dpi = Dpi) )
}


## Forms the optimal score matrix. Works with G >= 2 
## Input: Vector of length n of training categories.
## Output: (G x G-1) matrix of score vectors. 
OptScores <- function(TrainCat) {
  Y <- IndicatMat(TrainCat)$Categorical
  nclass = colSums(Y)
  n = nrow(Y)
  theta = NULL
  for (l in 1:(ncol(Y) - 1)) {
    temp = c(rep(sqrt(n * nclass[l + 1] / (sum(nclass[1:l]) * sum(nclass[1:(l + 1)])) ),
                 l), -sqrt(n * sum(nclass[1:l]/(nclass[l + 1] * sum(nclass[ 1:(l + 1) ])) )),
             rep(0, ncol(Y) - 1 - l))
    theta = cbind(theta, temp)
  }
  return(Optimal_Scores = theta)
}

### Computes the gaussian kernel evaluation of the point x with each of the sample points.  x
### is a data point with p variables Data is a n times p matrix SigmaK is the
### bandwidth parameter of the Gaussian kernel.
Kernel <- function(x, TrainData, Sigma) {
  if(Sigma <= 0 ) stop('Gaussian kernel parameter <= 0.')
  
  DiffPart <- (t(t(TrainData) - x))^2  ## Computes the distance squared of the data and point x
  DiffPart <- rowSums(DiffPart)  # Sum of squares
  exp( - DiffPart / Sigma)  #Divide by kernel parameter and evluate exponential function
}



# Forms the kernel matrix.  (i,j) component is kernel evluated on data points i and
# j Data is n x p. SigmaKm is Gaussian kernel parameter.
KernelMat <- function(TrainData, Sigma) {
  if(Sigma <= 0 ){ stop('Gaussian kernel parameter <= 0.') }
  S = tcrossprod(as.matrix(TrainData))  ## computes XX^t
  n = nrow(S)
  D = -2 * S + matrix(diag(S), n, n, byrow = T) + matrix(diag(S), n, n)  ## computes all pairs of squared distances
  K = exp( - D / Sigma )  #evaluates kernel on distances
  return(K)
}


# Forms weighted kernel matrix. Data is (n x p) w is weight vector of length p. Each
# coordinate must lie between -1 and 1.  SigmaKw is Gaussian kernel parameter.
KwMat <- function(TrainData, w, Sigma) {
  # Check Parameters 
  if(Sigma <= 0 ) stop('Gaussian kernel parameter <= 0.')
  if(any( abs(w) > 1) ) stop('A weight coordinate is outside [-1,1].')
  w <- as.numeric(w)
  TrainData <- t(as.matrix(TrainData)) * w  #multiplies each  coordinate of the data
  # matrix by corresponding weight.
  TrainData <- t(TrainData)
  return(KernelMat(TrainData, Sigma))
}




### Function which catorgizes data points in the construction of the toy data set If
### radius of point is below 2/3-.1, the point belongs to category 1.  If radius is
### greater than or equal to 2/3, the point is in category 2.  x is a vector of length
### p.
Categorize <- function(x) {
  c <- 0
  if (sqrt(sum(x^2)) < (2/3 - 0.1)) {
    c <- 1
  } else if (sqrt(sum(x^2)) >= 2/3) {
    c <- 2
  } else {
    c <- -1
  }
  c
}


GetProjections <- function(TestData, TrainData, TrainCat, Dvec = NULL, w = rep(1, ncol(TrainData)), Kw = NULL, Sigma = NULL, Gamma = NULL){
  n <- nrow(TrainData)
  
  # Generate Parameters if none are supplied
  if(is.null(Sigma)){
    output <- SelectParams(TrainData, TrainCat)
    Gamma <- output$Gamma
    Sigma <- output$Sigma
  }
  
  #Generate ridge parameter if not provided by kernel parameter is
  else if(is.null(Gamma)){
    output <- SelectParams(TrainData, TrainCat, Sigma = Sigma)
    Gamma <- output$Gamma
  }
  
  
  #Test Parameters for inadmissible values
  if( any(w > 1) || any(w < -1)) stop("Some weight vector is outside the interval [-1,1]")
  if(Sigma <= 0 ) stop('Gaussian kernel parameter <= 0.')
  
  #form weighted kernel matrix if null
  if(is.null(Kw)){
    Kw <- KwMat(TrainData, w, Sigma)
  }
  
  #Generate One-hot encoding matrix and optimal scores
  Y <- IndicatMat(TrainCat)$Categorical
  Theta <- OptScores(TrainCat)
  YTheta <- Y %*% Theta
  
  #Generate Discriminant Vectors
  if(is.null(Dvec)){
    Dvec <- SolveKOSCPP(YTheta, Kw, Gamma)
  }
  
  #Weight the Training and Testing Data with supplied weights
  TrainData <- t( t(TrainData) * w)
  TestData <- as.matrix(TestData)
  if(ncol(TestData) != ncol(TrainData)) TestData <- t(TestData)
  TestData <- t( t(TestData) * w)
  
  #Generate Projection Values 
  PV <- apply(TestData, MARGIN = 1, FUN = function(z){
    Kx <- t(Kernel(z, TrainData, Sigma))  # Kernel evalutations of x and each vector in Data
    M1 <- colMeans(Kw)  #Create centered vector to shift Kx by
    Dvec <- scale(Dvec, center = TRUE, scale = FALSE)
    P <- ( Kx - M1 ) %*% Dvec
    P <- as.numeric(P)
    P
  })
  return(PV)
}


# Forms Q matrx and B vector in qudratic form used to update weights.  Data is n x p
# A is a n x G-1 matrix of Discriminant vectors. In the two group case, it will be a
# vector of length n.  YTheta is the n x G-1 transformed response matrix.  w is
# weight vector of length p. Each coordinate lies between -1 and 1.  Gamma is
# ridge parameter Sigma is Gaussian kernel parameter.
FormQB <- function(TrainData, Dvec, YTheta, w, Sigma, Gamma) {
  Kw <- KwMat(as.matrix(TrainData), w, Sigma)  #forms weighted kernel matrix
  #Doubly-center the weighted kernel matrix
  Kw <- scale(Kw, center = TRUE, scale = FALSE)
  Kw <- scale(t(Kw), center = TRUE, scale = FALSE)
  
  #Set parameters
  p <- length(w)
  k <- ncol(as.matrix(Dvec))
  n <- nrow(Kw)
  
  #Initialize Q matrix and B vector 
  Q <- matrix(rep(0, p^2), nrow = p)  #initialize Q matrix.
  B <- rep(0, p)  #initialize B vector
  
  for (i in 1:k) {
    Dvec_col <- Dvec[, i]
    Tmat <- TMatCPP(TrainData, Dvec_col, w, Sigma)  #Forms T matrix
    Q <- (1/n) * crossprod(Tmat) + Q  #Creates t(CT)*CT and adds it to previous terms
    New <- (1/n) * crossprod(YTheta[, i] - Kw %*%  Dvec_col + Tmat %*% w, Tmat) - (Gamma/2) * crossprod(Dvec_col, Tmat)
    # Forms Component of B vector which uses discriminant vector Dvec[,i]
    B <- B + New  #Adds new component to sum of old components
  }
  return(list(B = B, Q = Q, Tmat = Tmat))
}


SparseKernOptScore <- function(TrainData, TrainCat, w0=rep(1, ncol(TrainData)), 
                               Sigma, Gamma, Lambda, Maxniter=100,
                               Epsilon = 1e-05, Error = 1e-05) {
  Y<-IndicatMat(TrainCat)$Categorical
  # Get Everything Initialized#
  D <- (ncol(Y) - 1)  #Number of discriminant vectors, G-1
  n <- nrow(TrainData)
  p <- ncol(TrainData)
  error <- 1  #Initialize Error value
  niter <- 1
  
  #for sequence of weights and initialize first value
  Weight_Seq <- matrix(0, length(w0), Maxniter)
  Weight_Seq[, niter] <- w0

  #Form optimal scores and weighted kernel matrix
  Opt_Scores <- OptScores(TrainCat)
  YTheta <- Y %*% Opt_Scores
  Kw <- KwMat(TrainData, w0,  Sigma)
  
  #If no sparsity penalty, return KOS vector and initial weights
  if(Lambda == 0){
    Dvec <- SolveKOSCPP(YTheta, Kw, Gamma)
    return(list(Weights = w0, Dvec = Dvec))
  }
  
  #From here on, sparsity penalty Lambda > 0. 

  # Create Initial Discrimiant Vector and Quadratic Form Terms#
  Dvec_old <- SolveKOSCPP(YTheta, Kw, Gamma, Epsilon) #initialize old coefficients

  OldQB <- FormQB(TrainData, Dvec = Dvec_old, YTheta = YTheta, w = w0, Gamma = Gamma, Sigma)
  Qold <- OldQB$Q #Initialize Q
  Bold <- OldQB$B #Initialize B

  # Record Objective Function Value with Initial Terms
  OFV <- ObjectiveFuncCPP(w0, Kw, TrainData, Dvec_old, YTheta,Lambda, Gamma,  Epsilon)
  OFV_seq <- rep(0, Maxniter)
  OFV_seq[niter] <- OFV
  
  # Enter while loop for updating weights and discriminant vectors.
  # Terminates if loss function values converge or if maximum 
  # iterations reached.
  while (error > Error && niter < Maxniter) {
    if (niter > Maxniter) { 
      break #additional test for maximum iterations reached.
    }
    
    # Update the Weights
    w <- as.numeric(CoordDesCPP(w0, Qold, Bold, Lambda, 1e-06, 1e+07))
    Kw <- KwMat(TrainData, w, Sigma) #reform weighted kernel matrix

    ## Tests for decrease objective function###
    OFV_Weights_new <- ObjectiveFuncCPP(w, Kw, TrainData, Dvec_old, YTheta, Lambda, Gamma, Epsilon)
    # if(OFV_Weights_new>OFV)print('Weights Increased Obj. Func.')

    #Update Discriminant Vector
    Dvec <- SolveKOSCPP(YTheta,  Kw, Gamma, Epsilon)
    if (sum(Dvec^2) == 0) {
      print("Discriminant Vector is zero.")
    }

    OFV_dvec <- ObjectiveFuncCPP(w, Kw, TrainData, Dvec, YTheta, Lambda, Gamma, Epsilon)

    ### Update quadratic form terms ###
    Output <- FormQB(TrainData = TrainData, 
                     Dvec = Dvec, 
                     YTheta = YTheta, 
                     w = w, 
                     Gamma = Gamma, 
                     Sigma = Sigma)
    Q <- Output$Q
    B <- Output$B

    Error <- abs(OFV_dvec - OFV)

    niter <- niter + 1 #update iteration count

    OFV_seq[niter] <- OFV_dvec
    Weight_Seq[, niter] <- w

    w0 <- w
    Bold <- B
    Qold <- Q
    Dvec_old <- Dvec
    OFV <- OFV_dvec
  }
  if (niter < Maxniter) {
    Weight_Seq[, c(niter:Maxniter)] <- w
    OFV_seq[c(niter:Maxniter)] <- OFV
  }
  return(list(Weights = w0, Dvec = Dvec))
}

SelectRidge <- function(TrainData, TrainCat, Sigma, Epsilon = 1e-05) {
  n <- nrow(TrainData)
  
  YTrain <- IndicatMat(TrainCat)$Categorical
  Opt_ScoreTrain <- OptScores(TrainCat)
  YThetaTrain <- YTrain %*% Opt_ScoreTrain

  K <- KernelMat(TrainData, Sigma)
  K <- scale(K, center = TRUE, scale = FALSE)
  K <- scale(t(K), center = TRUE, scale = FALSE)
  VS <- (n / ((n - 1)^2 * (n - 2))) * (sum(diag(K)^2) - (1/n) * sum(K^2))
  denom <- (1 / (n - 1)^2) * (sum(K^2))

  t <- VS / denom
  Gamma <- t / (1 - t)
  return(Gamma)
}


RidgeGCV <- function(Data, Cat, Sigma, Epsilon = 1e-05) {
  YTrain <- IndicatMat(Cat)$Categorical
  Opt_ScoreTrain <- OptScores(Cat)
  YThetaTrain <- YTrain %*% Opt_ScoreTrain

  K <- KernelMat(Data, Sigma)
  n <- nrow(K)
  C <- diag(n) - (1/n) * matrix(rep(1, n^2), nrow = n)
  M <- (C %*% K) %*% C
  Gammaseq <- seq(from = 0, to = 0.5, by = 1e-04)
  values <- sapply(Gammaseq, FUN = function(Gamma) {
    Mat <- MASS::ginv(M %*% M + Gamma * n * (M + Epsilon * diag(n))) %*% M
    Num <- sum((YThetaTrain - Mat %*% YThetaTrain)^2)
    Denom <- (n - sum(diag(Mat)))
    return(Num/Denom)
  })
  Gammaseq[which.min(values)]
}


## Cross Validation Code
LassoCV <- function(TrainData, TrainCat, B, Gamma, Sigma,
                    Epsilon = 1e-05) {
  c <- 2 * max(abs(B))
  Lambdaseq <- seq(from = 1e-10 * c, to = c, length.out = 20)
  
  n <- nrow(TrainData)
  FoldLabels <- CreateFolds(TrainCat)
  w <- rep(1, ncol(TrainData))
  Errors <- rep(0, 20)
  for (j in 1:20) {
    totalError <- 0
    # If sparsity parameter was too large, all weights are set to 0. Set misclassification
    # Error to be maximial
    if (sum(abs(w)) == 0) {
      totalError <- length(TrainCat)
    } 
    else {
      nfold <- 5
      for (i in 1:nfold) {
        ## Make Train and Validation Folds##
        NewTrainData <- subset(TrainData, FoldLabels != i)
        NewTrainCat <- subset(TrainCat, FoldLabels != i)
        NewTestCat <- subset(TrainCat, FoldLabels == i)
        NewTestData <- subset(TrainData, FoldLabels == i)
        
        ## Scale and Center Folds ##
        output <- CenterScale(NewTrainData, NewTestData)
        NewTrainData <- output$TrainData
        NewTestData <- output$TestData
        
        YTrain <- IndicatMat(NewTrainCat)$Categorical  #Create n x G categorical response matrix
        
        # Apply kernel feature selection algorithm on training data
        output <- SparseKernOptScore(NewTrainData, 
                                     NewTrainCat, 
                                     w0 = w, 
                                     Sigma = Sigma, 
                                     Gamma = Gamma,
                                     Lambda = Lambdaseq[j], 
                                     Maxniter = 100, 
                                     Epsilon = Epsilon)
        
        w <- as.numeric(output$Weights)
        
        if (sum(abs(w)) == 0) {
          totalError <- 0
        } 
        else {
          # Scale test data by weights
          
          NewTestDataFold <- t(t(NewTestData) * w)
          NewTrainDataFold <- t(t(NewTrainData) * w)
          
          ## Need weighted kernel matrix to compute projection values
          NewKtrain <- KernelMat(NewTrainDataFold, Sigma = Sigma)
          Dvec <- output$Dvec
          
          # Create Test Projection Values
          NewTestProjections <- GetProjections(TestData = NewTestDataFold, 
                                              TrainData = NewTrainDataFold, 
                                              TrainCat = NewTrainCat, 
                                              Dvec = Dvec, 
                                              w = w, 
                                              Kw = NewKtrain, 
                                              Sigma = Sigma, 
                                              Gamma = Gamma)
          
          ### Get Train Projection Values
          TrainProjections <- GetProjections(TestData = NewTrainDataFold, 
                                            TrainData = NewTrainDataFold, 
                                            TrainCat = NewTrainCat, 
                                            Dvec = Dvec, 
                                            w = w, 
                                            Kw = NewKtrain, 
                                            Sigma = Sigma, 
                                            Gamma = Gamma)
          
          ### All of this is used to create discirminant line
          OldData <- data.frame(as.numeric(NewTrainCat), as.numeric(TrainProjections))
          colnames(OldData) <- c("Category", "Projections")
          ## fit LDA on training projections
          LDAfit <- MASS::lda(Category ~ Projections, data = OldData)
          NewData <- data.frame(as.numeric(NewTestCat), as.numeric(NewTestProjections))
          colnames(NewData) <- c("Category", "Projections")
          # Predict class membership using LDA
          predictions <- stats::predict(object = LDAfit, newdata = NewData)$class
          # Compute number of misclassified points
          FoldError <- sum(abs(as.numeric(predictions) - NewTestCat))
          totalError <- totalError + FoldError
        }
      }
    }
    Errors[j] <- totalError / n
  }
  return(list(Lambda = Lambdaseq[which.min(Errors)], Errors = Errors))
}

# Code to select kernel, ridge, and sparsity parameters.
#' @title Generates parameters.
#' @param TrainData (n x p) Matrix of training data with numeric features. Cannot have missing values.
#' @param TrainCat (n x 1) Vector of class membership. Values must be either 1 or 2.
#' @param Sigma Scalar Gaussian kernel parameter. Default set to NULL and is automatically generated if user-specified value not provided. Must be > 0. User-specified parameters must satisfy hierarchical ordering.
#' @param Gamma Scalar ridge parameter used in kernel optimal scoring. Default set to NULL and is automatically generated if user-specified value not provided. Must be > 0. User-specified parameters must satisfy hierarchical ordering.
#' @param Epsilon Numerical stability constant with default value 1e-05. Must be > 0 and is typically chosen to be small.
#' @references 
#' Lapanowski, Alexander F., and Gaynanova, Irina. ``Sparse feature selection in kernel discriminant analysis via optimal scoring'', (preprint)
#' @references Lancewicki, Tomer. "Regularization of the kernel matrix via covariance matrix shrinkage estimation." arXiv preprint arXiv:1707.06156 (2017).
#' @references Lapanowski, Alexander F., and Gaynanova, Irina. ``Sparse feature selection in kernel discriminant analysis via optimal scoring'', Artificial Intelligence and Statistics, 2019.
#' @description Generates parameters to be used in sparse kernel optimal scoring.
#' @details Generates the gaussian kernel, ridge, and sparsity parameters for use in sparse kernel optimal scoring using the methods presented in [Lapanowski and Gaynanova, preprint]. 
#' The Gaussian kernel parameter is generated using five-fold cross-validation of the misclassification error rate aross the {.05, .1, .2, .3, .5} quantiles of squared-distances between groups. 
#' The ridge parameter is generated using a stabilization technique developed in [Lapanowski and Gaynanova, preprint].
#' The sparsity parameter is generated by five-fold cross-validation over a logarithmic grid of 20 values in an automatically-generated interval.
#' @importFrom stats quantile sd
#' @return A list of 
#' \item{Sigma}{ Gaussian kernel parameter.}  
#' \item{Gamma}{ Ridge Parameter.}
#' \item{Lambda}{ Sparsity parameter.}
#' @examples 
#' Sigma <- 1.325386  #Set parameter values equal to result of SelectParam.
#' Gamma <- 0.07531579 #Speeds up example
#' 
#' TrainData <- KOS_Data$TrainData
#' TrainCat <- KOS_Data$TrainCat
#' 
#' SelectParams(TrainData = TrainData , 
#'              TrainCat = TrainCat, 
#'              Sigma = Sigma, 
#'              Gamma = Gamma)
SelectParams <- function(TrainData, TrainCat, Sigma = NULL, Gamma = NULL, Epsilon = 1e-05) {
  n <- nrow(TrainData)
  p <- ncol(TrainData)
  
  if(is.null(Sigma) & is.null(Gamma)){
    E <- matrix(0, nrow = 5, ncol = 4)
    QuantileTest <- c(0.05, 0.1, 0.2, 0.3,.5)
    Data1 <- subset(TrainData , TrainCat == 1)
    Data2 <- subset(TrainData , TrainCat == 2)
    DistanceMat <- fields::rdist(x1 = Data1, x2 = Data2)
  
    Y <- IndicatMat(TrainCat)$Categorical
    Theta <- OptScores(TrainCat)
    YTheta <- Y %*% Theta
  
    for(j in 1:5){
      Sigma <- stats::quantile(DistanceMat, QuantileTest[j])
    
      Gamma <- SelectRidge(TrainData, TrainCat, Sigma, Epsilon)
      K <- KernelMat(TrainData, Sigma)
      Dvec <- SolveKOSCPP(YTheta, K, Gamma)
      B <- FormQB(TrainData, Dvec, YTheta, w = rep(1, p), Sigma, Gamma)$B
      output <- LassoCV(TrainData = TrainData, 
                        TrainCat = TrainCat, 
                        B = B, 
                        Gamma = Gamma, 
                        Sigma = Sigma, 
                        Epsilon = Epsilon)
      
      E[j, ] <- c(min(output$Errors), Gamma, output$Lambda, Sigma)
    }
    j <- which.min(E[, 1])
    return(list(Sigma = E[j, 4], Gamma = E[j, 2], Lambda = E[j, 3]))
  }
  else if( is.null(Sigma) == FALSE & is.null(Gamma) == FALSE){
    ### Need to Create Lambda Value ###
    Y <- IndicatMat(TrainCat)$Categorical
    Theta <- OptScores(TrainCat)
    YTheta <- Y %*% Theta
    K <- KernelMat(TrainData, Sigma)
    Dvec <- SolveKOSCPP( YTheta, K, Gamma)
    B <- FormQB(TrainData, 
                Dvec, 
                YTheta, 
                w=rep(1, p), 
                Gamma, 
                Sigma)$B
    
    Lambda <- LassoCV(TrainData = TrainData , 
                      TrainCat = TrainCat, 
                      B = B, 
                      Gamma = Gamma, 
                      Sigma = Sigma, 
                      Epsilon = Epsilon)$Lambda
  }
  else if(is.null(Sigma) == FALSE & is.null(Gamma)){
    Gamma <- SelectRidge(TrainData, TrainCat, Sigma, Epsilon)
    Y <- IndicatMat(TrainCat)$Categorical
    Theta <- OptScores(TrainCat)
    YTheta <- Y %*% Theta
    K <- KernelMat(TrainData, Sigma)
    Dvec <- SolveKOSCPP(YTheta, K, Gamma)
    
    B <- FormQB(TrainData, 
                Dvec, 
                YTheta, 
                w = rep(1, p), 
                Gamma, 
                Sigma)$B
    
    Lambda <- LassoCV(TrainData = TrainData, 
                      TrainCat = TrainCat, 
                      B = B, 
                      Gamma = Gamma, 
                      Sigma = Sigma)$Lambda
  }
  else{
    stop("Hierarchical order of parameters violated.")
  }
  return(list(Sigma = Sigma, Gamma = Gamma, Lambda = Lambda))
}


### Helper Function. Computes column means and standard deviations
### of Training Data matrix. Shifts and scales Test Data matrix
### columns by those values.
CenterScale<-function(TrainData, TestData){
  ColMeans<-apply(TrainData, MARGIN=2, FUN = mean)
  
  ColSD<-apply(TrainData, MARGIN=2, FUN = stats::sd)
  
  TrainData<-scale(TrainData, scale=T)
  
  TestData<-as.matrix(TestData)
  
  TestData<-t(apply(TestData,MARGIN=1, FUN=function(x) x - ColMeans ))
  
  for(j in 1:ncol(TrainData)){
    if(ColSD[j] != 0){
      TestData[,j] <- TestData[,j] / ColSD[j]
    }
  }
  return(list(TrainData = TrainData, TestData = TestData))
}

### Helper Function. Creates fold labels which maintain class proportions
CreateFolds <- function(Cat) {
  n <- length(Cat)
  Index <- c(1:n)
  Category <- data.frame(Cat = Cat, Index = Index)
  Cat1 <- subset(Category, Cat == 1)
  n1 <- nrow(Cat1)
  Cat2 <- subset(Category, Cat == 2)
  n2 <- nrow(Cat2)
  
  nfold <- 5
  FoldLabels1 <- cut( seq(1, n1), breaks = nfold, labels = FALSE)
  FoldLabels1 <- sample( FoldLabels1, size = n1, replace = FALSE)
  Cat1 <- cbind(Cat1, FoldLabel = FoldLabels1)
  
  FoldLabels2 <- cut( seq(1, n2), breaks = nfold, labels = FALSE)
  FoldLabels2 <- sample( FoldLabels2, size = n2, replace = FALSE)
  Cat2 <- cbind( Cat2, FoldLabel = FoldLabels2)
  
  Labels <- rbind(Cat1, Cat2)
  
  return(as.numeric(Labels$FoldLabel))
}



