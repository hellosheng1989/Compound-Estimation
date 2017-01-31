##Obtain the pointwise estimator from local regression up to degree J. c_{j,a} in paper
###Obtain the matrix of pointwise estimate dim = L_n * (J+1)
library(locfit)
library(locpol)
###function to obtain pointwise estimators: Matrix have dim = L_n*(J+1) 
### we obtain c_{j,a}, and j=0,1,...J and J<=locdeg
LocalPtMat = function(Yvec,Xvec,aVec,h_par, Svec, locdeg = 2,J=2){
  In_weights = 1/Svec^2
  L_n = length(aVec)
  PtMat = matrix(rep(0,L_n*(J+1)),nrow=L_n)
  for(i in 0:J){
    locObj = locfit(Yvec~Xvec,alpha=c(0,h_par),weights = In_weights,deriv=rep(1,i),deg = locdeg, kern = "gauss")
    PtMat[,(i+1)] = predict(locObj, aVec)
  }
  return(PtMat)
}