#Get_Wvec is for calculating weight function in paper
#x is a number, aVec and beta_nVec should be the same length
Get_Wvec = function(aVec,x,beta_nVec){
  Denominate = sum( exp(-beta_nVec*(x-aVec)^2) )
  Wvec = exp(-beta_nVec*(x-aVec)^2)/Denominate
  return(Wvec)
}

#derivative of Weight W
Get_Wderi = function(aVec,x,beta_nVec){
  A = sum( exp(-beta_nVec*(x-aVec)^2))
  B = sum( beta_nVec*(x-aVec)*exp(-beta_nVec*(x-aVec)^2))
  Wvec_Deri = (-2*beta_nVec*(x-aVec)*exp(-beta_nVec*(x-aVec)^2)*A 
               + 2*exp(-beta_nVec*(x-aVec)^2)*B)/A^2
  return(Wvec_Deri)
}

###calculate Compound estimator with any input pointwise estimators, PtMat is dim(L_n*(J+1)), and x is a value of prediction point.
CPEcompute = function(PtMat=PtMat,beta_nVec=beta_nVec,aVec=aVec,Yvec=Yvec,polyOd=polyOd,x=x,
                      Svec=Svec){
  L_n=length(aVec); J=polyOd-1
  polyMat = matrix(rep(0,(J+1)*L_n),ncol=L_n)
  for(j in 0:J){
    polyMat[j+1,] = (x-aVec)^j 
  }
  Wvec = Get_Wvec(aVec=aVec,x=x,beta_nVec=beta_nVec)
  polynoVec = diag(PtMat%*%polyMat)   # \sum_{j=0}^J c_{j,a}(x-a)^j
  CPE = sum(Wvec*polynoVec)
  return(CPE)
}


##first derivative - compound estimation with any input pointwise estimators
CPEcompute_D1 = function(PtMat=PtMat,beta_nVec=beta_nVec,aVec=aVec,Yvec=Yvec,polyOd=polyOd,x=x,
                         Svec=Svec){
  L_n=length(aVec); J=polyOd-1
  polyMat = matrix(rep(0,(J+1)*L_n),ncol=L_n)
  for(j in 0:J){
    polyMat[j+1,] = (x-aVec)^j 
  }
  PtMat1 = PtMat[,-1]
  polyMat1 = matrix(rep(0,J*L_n),ncol=L_n)
  for(j in 1:J){
    polyMat1[j,] = j*(x-aVec)^(j-1)
  }
  polynoVec = diag(PtMat%*%polyMat)   # \sum_{j=0}^J c_{j,a}(x-a)^j
  polynoVecD = diag(PtMat1%*%polyMat1)   # \sum_{j=1}^J c_{j,a} j(x-a)^(j-1)
  Wvec = Get_Wvec(aVec=aVec,x=x,beta_nVec=beta_nVec)
  WvecD = Get_Wderi(aVec=aVec,x=x,beta_nVec=beta_nVec)
  CPE_D1 = sum(Wvec*polynoVecD + WvecD*polynoVec)
  return(CPE_D1)
}


###############calculate compound estimation for a vector Xvec
###PtMat is the matrix of pointwise estimate dim = L_n * (J+1), Xvec is a vector of covariates, xvec is a vector of prediction points. 
CPEcomputeVecG = function(PtMat=PtMat,beta_nVec=beta_nVec,aVec=aVec,Yvec=Yvec,polyOd=polyOd,xvec=Xvec,
                         Svec=Svec){
  L_n=length(aVec); J=polyOd-1 ; n=length(xvec)
  MatXdA = apply(as.array(xvec),1,function(x){x-aVec})   
  BigMatC = matrix(rep(0,n*L_n),ncol=n)
  for(j in 0:J){
    TempBigMatC = PtMat[,j+1]*(MatXdA^j)
    BigMatC = BigMatC + TempBigMatC
  }
  BigWmatVec = apply(as.array(xvec),1,Get_Wvec,aVec=aVec,beta_nVec=beta_nVec)
  BigWmat = matrix(BigWmatVec,ncol=L_n,byrow = TRUE)
  CPEvec = rowSums(BigWmat*t(BigMatC))
  return(CPEvec)
}

###PtMat is the matrix of pointwise estimate dim = L_n * (J+1)
##################### J should >=2
###first derivative - compound estimation with same bandwidth on a vector xvec = Xvec
CPEcompute_D1VecG = function(PtMat=PtMat,beta_nVec=beta_nVec,aVec=aVec,Yvec=Yvec,polyOd=polyOd,xvec=Xvec,
                            Svec=Svec){
  L_n=length(aVec); J=polyOd-1;n=length(xvec)
  MatXdA = apply(as.array(xvec),1,function(x){x-aVec})   
  BigMatC = matrix(rep(0,n*L_n),ncol=n)
  for(j in 0:J){
    TempBigMatC = PtMat[,j+1]*(MatXdA^j)
    BigMatC = BigMatC + TempBigMatC
  }
  PtMat1 = as.matrix(PtMat[,-1])
  BigMatC1 = matrix(rep(0,n*L_n),ncol=n)
  for(j in 1:J){
    TempBigMatC1 = PtMat1[,j]*j*(MatXdA^(j-1))
    BigMatC1 = BigMatC1 + TempBigMatC1
  }
  BigWmatVec = apply(as.array(xvec),1,Get_Wvec,aVec=aVec,beta_nVec=beta_nVec)
  BigWmat = matrix(BigWmatVec,ncol=L_n,byrow = TRUE)
  #########################
  BigWmatVecD = apply(as.array(xvec),1,Get_Wderi,aVec=aVec,beta_nVec=beta_nVec)
  BigWmatD = matrix(BigWmatVecD,ncol=L_n,byrow = TRUE)
  CPE_D1vec = rowSums(BigWmatD*t(BigMatC) + BigWmat*t(BigMatC1))
  return(CPE_D1vec)
}