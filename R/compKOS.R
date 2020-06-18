## Compressed KOS Code
formCompKMat <- function(TrainData, TrainCat, m1, m2, s, Sigma){
  TrainX1 <- TrainData[TrainCat == 1, ]
  TrainX2 <- TrainData[TrainCat == 2, ]
  
  #---- Initialize compression matrices -----
  n1 <- sum(TrainCat == 1)
  n2 <- sum(TrainCat == 2)
  
  Q1 <- createSketchMatrix(n = n1, m = m1, s = s)
  Q2 <- createSketchMatrix(n = n2, m = m2, s = s)
  
  #---- Apply Helper functions for diagonal blocks ------

  K11 <- kernelMatDiagHelper(TrainData = TrainX1, Q = Q1, Sigma = Sigma)
  K22 <- kernelMatDiagHelper(TrainData = TrainX2, Q = Q2, Sigma = Sigma)
  
  #----- Apply Helper functions for off-diagonal blocks ------
  K12 <- kernelMatOffDiagHelper(TrainData = TrainData,
                                TrainCat = TrainCat,
                                Q1 = Q1,
                                Q2 = Q2,
                                Sigma = Sigma)
  Top <- cbind(K11, K12)
  Bottom <- cbind(t(K12), K22)
  return(rbind(Top, Bottom))
} 

#A helper function for compressing the kernel matrix.
#Deals with only one class of trianing data.
kernelMatDiagHelper <- function(TrainData, Q, Sigma){
  n <- nrow(TrainData)
  m <- nrow(Q)
  M <- matrix(0, nrow = m, ncol = n)
  compK <- matrix(0, nrow = m, ncol = m)
  #--- Iterate and compress each column of kernel matrix---
  for(i in 1:n){
    x <- Kernel(x = TrainData[i, ], TrainData = TrainData, Sigma = Sigma)
    M[ ,i] <- as.numeric(Q %*% as.numeric(x))
  }
  
  #--- Compress other side ----
  M <- t(M)
  compK <- Q %*% M
  return(compK)
}

#A helper function to compute the off diagonal
#blocks of the compressed kernel matrix
kernelMatOffDiagHelper <- function(TrainData, TrainCat, Q1, Q2, Sigma){
  TrainX1 <- TrainData[TrainCat == 1, ]
  TrainX2 <- TrainData[TrainCat == 2, ]
  n1 <- nrow(TrainX1)
  n2 <- nrow(TrainX2)
  m1 <- nrow(Q1)
  m2 <- nrow(Q2)
  
  M <- matrix(0, nrow = m1, ncol = n2)
  
  for(i in 1:n2){
    x <- Kernel(x = TrainX2[i, ], TrainData = TrainX1, Sigma)
    M[ ,i] <- as.numeric(Q1 %*% x)
  }
  M <- t(M)
  compK <- t(as.matrix(Q2 %*% M))
  return(compK)
}




#-------- Compressed Kernel Optimal Scoring Implementation --------




GetProjection <- function(X, Data, Cat, Dvec, K, Sigma , Gamma) {
  
  YTheta <- Y %*% Theta
  
  n <- nrow(Kw)
  PV<-apply(X, MARGIN = 1, FUN = function(z){
    Kx <- t(Kernel(x = z, TrainData = Data, Sigma = Sigma))  # Kernel evalutations of x and each vector in Data
    M1 <- colMeans(Kw)  #Create centered vector to shift Kx by
    Dvec <- scale(Dvec, center = TRUE, scale = FALSE)
    P <- ( Kx - M1 ) %*% Dvec
    P<-as.numeric(P)
    P
  })
  return(PV)
}



compressedKOS <- function(Kcomp, m1, m2, gamma, epsilon = 1e-5){
  m <- m1 + m2
  if(m != nrow(Kcomp)) stop("Error: m1+m2 does not equal the number of rows of Kcomp")
  Kcomp <- scale(Kcomp, scale = F, center = T) 
  Kcomp <- t(scale(t(Kcomp), center = T, scale = F))

  #Generate One-hot encoding matrix and optimal scores
  Y <- IndicatMat(c(rep(1, m1), rep(2, m2)))$Categorical
  Theta <- OptScores(c(rep(1, m1), rep(2, m2)))
  YThetaComp <- Y %*% Theta
  
  Dvec <- solve(a = Kcomp %*% Kcomp+m*gamma*(Kcomp+epsilon*diag(m)), 
                b = Kcomp %*% YThetaComp)
  return(Dvec)
}


GetProjection <- function(x, Data, Cat, Q1, Q2, DVec, Kcomp , Sigma){
  m1 <- nrow(Q1)
  m2 <- nrow(Q2)
  KVec1 <- Kernel(x = x, TrainData = subset(Data, Cat == 1), Sigma = Sigma)
  KVec2 <- Kernel(x = x, TrainData = subset(Data, Cat == 2), Sigma = Sigma)
  Left1 <- as.numeric(tcrossprod(t(KVec1), Q1)+rep(mean(KVec1), m1))
  Left2 <- as.numeric(tcrossprod(t(KVec2), Q2)+rep(mean(KVec2), m2))
  LeftValue <- c(Left1, Left2) %*% DVec
  RightValue <- colMeans(Kcomp) %*% DVec
  return(as.numeric(LeftValue-RightValue))
}

#---- Nystrom Implementation -----

# Returns the n x m matrix and m x m sub-sampled matrix
NystromMat <- function(Data, m, Sigma){
  n <- nrow(Data)
  Set <- sample(x = c(1:n), size = m, replace = F)
  SubSample <- Data[Set, ]
  Km <- KernelMat(TrainData = SubSample, Sigma = Sigma)
  Km <- MASS::ginv(Km)
  Km <- expm::sqrtm(Km)
  Rectangle <- matrix(0, nrow = n, ncol =m)
  for(i in 1:m){
    Rectangle[,i] <- Kernel(x = SubSample[i,], TrainData = Data, Sigma = Sigma)
  }
  return(list(Rectangle , Km))
}


NystromKOS <- function(TrainData, TrainCat, m, gamma, Sigma){
  #Generate One-hot encoding matrix and optimal scores
  Y <- IndicatMat(TrainCat)$Categorical
  Theta <- OptScores(TrainCat)
  YTheta <- Y %*% Theta
 
  n1 <- length(TrainCat[TrainCat==1])
  n2 <- length(TrainCat[TrainCat==2])
  n <- nrow(TrainData)
  Set <- sample(x = c(1:n), size = m, replace = F)
  SubSample <- TrainData[Set, ]
  Km <- KernelMat(TrainData = SubSample, Sigma = Sigma)
  Rectangle <- matrix(0, nrow = n, ncol = m)
  for(i in 1:m){
    Rectangle[,i] <- Kernel(x = SubSample[i,], TrainData = TrainData, Sigma = Sigma)
  }
  InvMat <- gamma*Km+crossprod(Rectangle, Rectangle)
  InvMat <- MASS::ginv(InvMat)
  Vec <- crossprod(Rectangle ,YTheta)
  Dvec <- InvMat %*% Vec
  return(list(Dvec = Dvec, Rec = Rectangle, sample = SubSample))
}

NyProject <- function(x, Sample, Rectangle, Dvec, Sigma){
  CenterTerm <- colSums(Rectangle)%*%Dvec
  return(crossprod(Kernel(x = x, TrainData = Sample, Sigma = Sigma ),Dvec) - CenterTerm)
}
