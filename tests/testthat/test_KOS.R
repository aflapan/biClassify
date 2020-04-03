
#--- Testing Kernel Matrix Functionality ---
test_that("KernelMat function equals Gaussian kernel matrix", {
 Data <- KOS_Data$TrainData
 n <- nrow(Data)
 sigma <- 1
 K <- KernelMat(TrainData = Data, Sigma = sigma)
 
 K_test <- matrix(0, nrow = n, ncol = n)
 for(i in 1:n){
   for(j in 1:n){
     K_test[i,j] <- exp(- sum((Data[i,] - Data[j,])^2))
   }
 }
 expect_equal(K, K_test)
})

test_that("Weighted Kernel Matrix equals kernel matrix when weights are 1", {
  Data <- KOS_Data$TrainData
  n <- nrow(Data)
  sigma <- 1
  K <- KernelMat(TrainData = Data, Sigma = sigma)
  Kw <- KwMat(TrainData = Data, w = rep(1, ncol(Data)), Sigma = sigma)
  expect_equal(K, Kw)
})

test_that("Weighted Kernel Matrix equals all 1s when weights are 0", {
  Data <- KOS_Data$TrainData
  n <- nrow(Data)
  sigma <- 1
  Kw <- KwMat(TrainData = Data, w = rep(0, ncol(Data)), Sigma = sigma)
  ones <- matrix(1, nrow = n, ncol = n)
  expect_equal(ones, Kw)
})

#--- Testing the Optimal Scoring Function ---

test_that("Optimal Scores are +- 1 for equal class sizes", {
  TrainCat <- c(rep(1, 10), rep(2,10))
  scores <- as.numeric(OptScores(TrainCat))
  test_scores <- c(1,-1)
  expect_equal(scores, test_scores)
})

# --- Testing GetProjections ---
test_that("GetProjections produces a projection for every test sample",{
  TrainData <- KOS_Data$TrainData
  TrainCat <- KOS_Data$TrainCat
  TestData <- KOS_Data$TestData
  TestCat <- KOS_Data$TestCat
  
  Projections <- GetProjections(TestData = TestData, 
                                  TrainData = TrainData, 
                                  TrainCat = TrainCat,
                                  Sigma = 1,
                                  Gamma = 0.1)
  
  expect_equal(length(Projections), nrow(TestData))
})

# --- Testing FormQB ---
test_that("FormQB returns correct dimensions",{
  TrainData <- KOS_Data$TrainData
  TrainCat <- KOS_Data$TrainCat
  TestData <- KOS_Data$TestData
  TestCat <- KOS_Data$TestCat
  
  p <- ncol(TrainData)
    
  #Generate One-hot encoding matrix and optimal scores
  Y <- IndicatMat(TrainCat)$Categorical
  Theta <- OptScores(TrainCat)
  YTheta <- Y %*% Theta
  
  Kw <- KwMat(TrainData, rep(1,p), 1)
  Dvec <- SolveKOSCPP(YTheta, Kw, 0.1)
  
  QBoutput <- FormQB(TrainData = TrainData, 
                     Dvec = Dvec, 
                     YTheta = YTheta, 
                     w = rep(1, ncol(TrainData)), 
                     Sigma = 1, 
                     Gamma = 0.1)
  testQ <- QBoutput$Q
  testB <- QBoutput$B
  testTmat <- QBoutput$Tmat
  
  expect_equal(c(nrow(testQ), length(testB), ncol(testTmat)), rep(p, 3))
})

# --- Testing SparseKernOptScore ---

test_that("SparseKernOptScore equals kernel ridge regression when Lambda = 0",{
  TrainData <- KOS_Data$TrainData
  TrainCat <- KOS_Data$TrainCat
  p <- ncol(TrainData)
  
  SparseKOS_output <- SparseKernOptScore(TrainData = TrainData,
                                         TrainCat = TrainCat,
                                         Sigma = 1,
                                         Gamma = 0.1,
                                         Lambda = 0)
  Y <- IndicatMat(TrainCat)$Categorical
  Theta <- OptScores(TrainCat)
  YTheta <- Y %*% Theta
  Kw <- KwMat(TrainData, rep(1,p), 1)
  Dvec <- SolveKOSCPP(YTheta, Kw, 0.1)
  
  expect_equal(Dvec, SparseKOS_output$Dvec)
})

test_that("SparseKernOptScore outputs list with correct dimensions",{
  TrainData <- KOS_Data$TrainData
  TrainCat <- KOS_Data$TrainCat
  p <- ncol(TrainData)
  n <- nrow(TrainData)
  
  Sigma <- 1
  Gamma <- 0.1
  Lambda <- 0.01
  
  output <- SparseKernOptScore(TrainData = TrainData,
                              TrainCat = TrainCat,
                              Sigma = Sigma,
                              Gamma = Gamma,
                              Lambda = Lambda)
  expect_equal(c(length(output$Dvec),length(output$Weights)), c(n,p))
})

# --- Test Select Ridge ---
test_that("Checking that the SelectRidge Parameter Gamma > 0 ",{
  TrainData <- KOS_Data$TrainData
  TrainCat <- KOS_Data$TrainCat
  p <- ncol(TrainData)
  n <- nrow(TrainData)
  
  Gamma <- SelectRidge(TrainData = TrainData,
                       TrainCat = TrainCat,
                       Sigma = 1)

  expect_equal(Gamma > 0, TRUE)
})



# --- Test SelectParams ---
test_that("Checking that all parameters are positive", {
  TrainData <- KOS_Data$TrainData
  TrainCat <- KOS_Data$TrainCat
  output <- SelectParams(TrainData = TrainData,
                         TrainCat = TrainCat)
  
  test_sigma <- output$Sigma
  test_ridge <- output$Gamma
  test_lambda <- output$Lambda
  
  expect_equal(all(test_sigma > 0, test_ridge > 0, test_lambda > 0), TRUE)
})