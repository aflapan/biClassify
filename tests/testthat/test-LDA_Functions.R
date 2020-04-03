

library(DAAG)
library(datasets)


test_that("Testing dimensions of compression matrix", {
  n <- 100000
  m <- 1000
  Q <- createSketchMatrix(n = n, m = m, s = 0.01)
  n_test <- ncol(Q)
  m_test <- nrow(Q)
  expect_equal(n_test, n)
  expect_equal(m_test, m)
})

test_that("Testing discrimiant vector if Null cov matrix is passed into function",{
  TrainCat <- leafshape17$arch + 1
  TrainData = data.matrix(leafshape17[, c(5, 7)])
  W_test <- formWithinGroupCov( TrainData = TrainData,
                                TrainCat = TrainCat)
  DiscrimVec1 <- formDiscrimVector(TrainData = TrainData,
                                   TrainCat = TrainCat,
                                   W = W_test)
  DiscrimVec2 <- formDiscrimVector(TrainData = TrainData,
                                   TrainCat = TrainCat)
  expect_equal(DiscrimVec1, DiscrimVec2)
})


test_that("Testing classification compared to MASS LDA on leaf data", {
  TrainCat <- leafshape17$arch + 1
  leaf17.lda <- MASS::lda(arch ~ logwid+loglen, data = leafshape17)
  MASSclass <- as.numeric(predict(leaf17.lda)$class)

  output <- Classify(TrainData = data.matrix(leafshape17[, c(5, 7)]),
                           TrainCat = TrainCat,
                           TestData = data.matrix(leafshape17[, c(5, 7)]),
                           gamma = 0)
  
  LDAclass <- output$Predictions
  expect_equal(LDAclass , MASSclass)
})


test_that("Check that Compressed Predict does not return NA",{
  TrainCat <- leafshape17$arch + 1
  TrainData = data.matrix(leafshape17[, c(5, 7)])
  output <- compressPredict(TrainData = TrainData,
                                TrainCat = TrainCat,
                                TestData = TrainData,
                                m1 = 10,
                                m2 = 10,
                                s = 1/2)
  compLabels <- output$Predictions
  expect_equal(any(is.na(compLabels)), FALSE)
})



test_that("Check that Subset Predict does not return NA",{
  TrainCat <- leafshape17$arch + 1
  TrainData = data.matrix(leafshape17[, c(5, 7)])
  output <- subsetPredict(TrainData = TrainData,
                             TrainCat = TrainCat,
                             TestData = TrainData,
                             m1 = 10,
                             m2 = 10)
  subLabels <- output$Predictions
  expect_equal(any(is.na(subLabels)), FALSE)
})



test_that("Check that subsetPredict equals Classify if subset entire data", {
  TrainCat <- leafshape17$arch + 1
  TrainData = data.matrix(leafshape17[, c(5, 7)])
  n1 <- sum(TrainCat == 1)
  n2 <- sum(TrainCat == 2)
  subLabels <- subsetPredict(TrainData = TrainData,
                             TrainCat = TrainCat,
                             TestData = TrainData,
                             m1 = n1,
                             m2 = n2)
  subLabels <- output$Predictions
  
  output <- Classify(TrainData = TrainData,
                     TrainCat = TrainCat,
                     TestData = TrainData)
  Labels <- output$Predictions
  
  expect_equal(Labels , subLabels)
})

test_that("checking that projectPredict equals LDA when compression matrix is identity",{
  TrainCat <- as.numeric(leafshape17$arch + 1)
  TrainData = data.matrix(leafshape17[, c(5, 7)])
  n1 <- sum(TrainCat == 1)
  n2 <- sum(TrainCat == 2)
  Q1 <- diag(n1)
  Q2 <- diag(n2)
  output <- projectPredict(TrainData = TrainData,
                              TrainCat = TrainCat,
                              TestData = TrainData,
                              Q1 = Q1,
                              Q2 = Q2,
                              m1 = n1,
                              m2 = n2,
                              s = 0.01,
                              gamma = 0)
  ProjClass <- output$Predictions
  
  leaf17.lda <- MASS::lda(arch ~ logwid+loglen, data = leafshape17)
  MASSclass <- as.numeric(predict(leaf17.lda)$class)
  expect_equal(MASSclass, ProjClass)
})

