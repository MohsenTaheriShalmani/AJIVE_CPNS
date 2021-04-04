library(RiemBase)

# Math functions
# function to calculate Mahalanobis distance (Hotelling Metric)
MahalanobisDistance<-function(X,Y){
  k<-dim(X)[2]
  nx<-dim(X)[1]
  ny<-dim(Y)[1]
  Sx<-cov(X)
  Sy<-cov(Y)
  meanX<-colMeans(X)
  meanY<-colMeans(Y)
  n<-nx+ny-1
  S<-((nx-1)*Sx+(ny-1)*Sy)/(nx+ny-2) #S=pooled covariance matrix
  T2<-t(meanX-meanY)%*%solve(S*(1/nx+1/ny))%*%(meanX-meanY)
  return(T2)
}

# Hotelling T2 test
HotellingT2<-function(X,Y){
  if(dim(X)[2]!=dim(Y)[2]){
    cat("Dimention Error!\n")
    break
  }
  k<-dim(X)[2]
  nx<-dim(X)[1]
  ny<-dim(Y)[1]
  Sx<-cov(X)
  Sy<-cov(Y)
  meanX<-colMeans(X)
  meanY<-colMeans(Y)
  n<-nx+ny-1
  S<-((nx-1)*Sx+(ny-1)*Sy)/(nx+ny-2) #S=pooled covariance matrix
  T2<-t(meanX-meanY)%*%solve(S*(1/nx+1/ny))%*%(meanX-meanY)
  F_value<-((n-k)/(k*(n-1)))*T2
  df1<-k
  df2<-n-k
  p_value<-1-pf(F_value,df1,df2)
  return(p_value)
}

# library(BisRNA)
# The fisher.method function is also available in BisRNA library
fisher.method<-function (pvalues)
{
  df <- 2 * length(pvalues)
  global_pValue<-pchisq(-2 * sum(log(pvalues), na.rm = TRUE), df, lower.tail = FALSE)
  # global_pValue<-1-pchisq(-2 * sum(log(pvalues), na.rm = TRUE), df, lower.tail = TRUE)
  return(global_pValue)
}

# convert vectors to unit vectors
convertVec2unitVec <- function(vec) {
  if(norm(vec,type = "2")==0){
    stop("vector is zero!")
  }
  return(vec/norm(vec,type = "2"))
}

# cross product of 2 vectors
myCrossProduct <- function(v,u) {
  return(c(v[2]*u[3]-v[3]*u[2],v[3]*u[1]-v[1]*u[3],v[1]*u[2]-v[2]*u[1]))
}

# frechet mean

frechetMean <- function(directions) {
  
  allDirTemp<-t(directions)
  data1 <- list()
  for (j in 1:dim(allDirTemp)[1]){
    data1[[j]] <-allDirTemp[j,]
  }
  data2 <- riemfactory(data1, name="sphere")
  ### Compute Fre'chet Mean
  out1<- rbase.mean(data2)
  meanFrechet<-as.vector(out1$x)
  
  return(meanFrechet)
  
}


# calculate unit normal vector of triangle mesh
unitNormalOfTriangle <- function(point1,point2,point3) {
  a<-point2-point1
  b<-point3-point1
  
  normalVec<-c((a[2]*b[3]-a[3]*b[2]),-(a[1]*b[3]-a[3]*b[1]),(a[1]*b[2]-a[2]*b[1]))
  triangleArea<-sqrt(sum(normalVec^2))/2
  unitNormal<-normalVec/sqrt(sum(normalVec^2))
  
  return(unitNormal)
  
}

