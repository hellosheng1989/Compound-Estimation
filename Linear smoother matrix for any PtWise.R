source("CPE with any PtWise.R")
###### G_D1Matrix calculate the linear smoother for compound estimation of 1st derivative
#####Yvec ~ Xvec
#####polyOd is the order J we want for pointwise estimators.
####xvec is 
G_D1Matrix = function(PtMat=PtMat,beta_nVec=beta_nVec,aVec=aVec,Yvec=Yvec,Xvec=Xvec,polyOd=polyOd,Svec=Svec){
  n=length(Yvec)
  TempYmat = diag(n)
  GMatD1=apply(TempYmat,1,CPEcompute_D1VecG,PtMat=PtMat,beta_nVec=beta_nVec,aVec=aVec,Yvec=Yvec,polyOd=polyOd,xvec=Xvec,
               Svec=Svec)
  return(GMatD1)
}

#######G_Matrix calculate the linear smoother for compound estimation 
G_Matrix = function(PtMat=PtMat,beta_nVec=beta_nVec,aVec=aVec,Yvec=Yvec,Xvec=Xvec,polyOd=polyOd,Svec=Svec){
  n=length(Yvec)
  TempYmat = diag(n)
  GMat=apply(TempYmat,1,CPEcomputeVecG,PtMat=PtMat,beta_nVec=beta_nVec,aVec=aVec,Yvec=Yvec,polyOd=polyOd,xvec=Xvec,
             Svec=Svec)
  return(GMat)
}