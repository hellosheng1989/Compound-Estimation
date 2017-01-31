source("CPE with any PtWise.R")
###### G_D1Matrix calculate the linear smoother for compound estimation of 1st derivative
#####Yvec ~ Xvec 
G_D1Matrix = function(h=0.5,beta_nVec=beta_nVec,aVec=aVec,Yvec=Yvec,Xvec=Xvec,polyOd=polyOd,xvec=Xvec,Svec=Svec){
  n=length(Yvec)
  TempYmat = diag(n)
  GMatD1=apply(TempYmat,1,CPEcompute_D1VecG,cMat=cMat,beta_nVec=beta_nVec,
               aVec=aVec,Xvec=Xvec,polyOd=polyOd,xvec=Xvec,Svec=Svec)
  return(GMatD1)
}

#######G_Matrix calculate the linear smoother for compound estimation 
G_Matrix = function(h=0.5,beta_nVec=beta_nVec,aVec=aVec,Yvec=Yvec,Xvec=Xvec,polyOd=polyOd,xvec=Xvec,Svec=Svec){
  n=length(Yvec)
  TempYmat = diag(n)
  GMatD1=apply(TempYmat,1,CPEcomputeVecG,cMat=cMat,beta_nVec=beta_nVec,
               aVec=aVec,Xvec=Xvec,polyOd=polyOd,xvec=Xvec,Svec=Svec)
  return(GMatD1)
}