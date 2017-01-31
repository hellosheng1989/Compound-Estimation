
### Spline First Derivative

library(sfsmisc)



Lmatrix<-function(x,h)
{
  n<-length(x)
  Lmatrix<-matrix(0,n,n)
  for (k in 1:n)
  {
    tempy <- rep(0,n)
    tempy[k] <- 1
    Lmatrix[,k] <- predict(smooth.spline(x=x,y=tempy,spar=h),x,deriv=0)$y
    
  }
  return(Lmatrix)
}


L1matrix<-function(x,h)
{
  n<-length(x)
  L1matrix<-matrix(0,n,n)
  for (k in 1:n)
  {
    tempy <- rep(0,n)
    tempy[k] <- 1
    L1matrix[,k] <- predict(smooth.spline(x=x,y=tempy,spar=h),x,deriv=1)$y
  }
  return(L1matrix)
}



L2matrix<-function(x,h)
{
  n<-length(x)
  L2matrix<-matrix(0,n,n)
  for (k in 1:n)
  {
    tempy <- rep(0,n)
    tempy[k] <- 1
    L2matrix[,k] <- predict(smooth.spline(x=x,y=tempy,spar=h),x,deriv=2)$y
  }
  return(L2matrix)
}
