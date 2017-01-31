#EmpricalDri is for calcualting the fisrt empirical derivative when there is no heteroskedasticity
# k is the order of the empirical derivatives, sigma is the error standard deviation,Xvec and Yvec are the covariates and responses
EmpiricalDri<-function(Xvec=Xvec,sigma=sigma,Yvec=Yvec,k=k){
     n<-length(Xvec)
     Svec = rep(sigma,n)
     EmpiY = rep(0,n)
     Cmat = matrix(rep(0,n*n),ncol=n)
     #i start from 3 to avoid boundary issue
     for(i in 2:(n-1)){
       k1 = min(k,i-1,n-i)
       Xvec1 = Xvec[(i+1):(i+k1)]; Xvec2 = Xvec[(i-1):(i-k1)]
       Svec1 = Svec[(i+1):(i+k1)]; Svec2 = Svec[(i-1):(i-k1)]
       Yvec1 = Yvec[(i+1):(i+k1)]; Yvec2 = Yvec[(i-1):(i-k1)]
       wVec = (c(1:k1))^2/(k1*(k1+1)*(2*k1+1)/6)
       if(i>=k+1&i<=n-k){
         Cmat[i,(i-1):(i-k)] = -wVec/(Xvec1-Xvec2)
         Cmat[i,(i+1):(i+k)] = wVec/(Xvec1-Xvec2)
       }
       if(i<=k){
         Cmat[i,(i-1):(i-k1)] = -wVec/(Xvec1-Xvec2)
         Cmat[i,(i+1):(i+k1)] = wVec/(Xvec1-Xvec2)
       }
       if(i>n-k){
         Cmat[i,(i-1):(i-k1)] = -wVec/(Xvec1-Xvec2)
         Cmat[i,(i+1):(i+k1)] = wVec/(Xvec1-Xvec2)
       }
       EmpiY[i] = sum(wVec* ((Yvec1-Yvec2)/(Xvec1-Xvec2)) )
     }
     EmpiY[1] = EmpiY[2]; EmpiY[n] = EmpiY[n-1]
     return(list(EmpiY=EmpiY,Cmat=Cmat))
}



#EmpricalDri2 is for calcualting the fisrt empirical derivative when there is heteroskedasticity
# k is the order of the empirical derivatives, with heteroskedasticity
# Xvec and Yvec are the values of covariate and responses. Svec is a vector of error standard deviations.
EmpiricalDri2<-function(Xvec=Xvec,Svec=Svec,Yvec=Yvec,k=k){
  n<-length(Xvec)
  EmpiY = rep(0,n)
  Cmat = matrix(rep(0,n*n),ncol=n)
  #i start from 3 to avoid boundary issue
  for(i in 2:(n-1)){
    k1 = min(k,i-1,n-i)
    Xvec1 = Xvec[(i+1):(i+k1)]; Xvec2 = Xvec[(i-1):(i-k1)]
    Svec1 = Svec[(i+1):(i+k1)]; Svec2 = Svec[(i-1):(i-k1)]
    Yvec1 = Yvec[(i+1):(i+k1)]; Yvec2 = Yvec[(i-1):(i-k1)]
    Denom = sum( (Xvec1-Xvec2)^2 / (Svec1^2+Svec2^2) )
    wVec = ( (Xvec1-Xvec2)^2 / (Svec1^2+Svec2^2) )/Denom
    if(i>=k+1&i<=n-k){
      Cmat[i,(i-1):(i-k)] = -wVec/(Xvec1-Xvec2)
      Cmat[i,(i+1):(i+k)] = wVec/(Xvec1-Xvec2)
    }
    if(i<=k){
      Cmat[i,(i-1):(i-k1)] = -wVec/(Xvec1-Xvec2)
      Cmat[i,(i+1):(i+k1)] = wVec/(Xvec1-Xvec2)
    }
    if(i>n-k){
      Cmat[i,(i-1):(i-k1)] = -wVec/(Xvec1-Xvec2)
      Cmat[i,(i+1):(i+k1)] = wVec/(Xvec1-Xvec2)
    }
    EmpiY[i] = sum(wVec* ((Yvec1-Yvec2)/(Xvec1-Xvec2)) )
  }
  EmpiY[1] = EmpiY[2]; EmpiY[n] = EmpiY[n-1]
  return(list(EmpiY=EmpiY,Cmat=Cmat))
}



