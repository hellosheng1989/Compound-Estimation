source("EmpiFunc.R") 
# Cmat from empirical derivative 
# DeriEst are the estimated derivatives
#EmpiY are the empirical derivatives
#n is sample size
#lmat is the linear smoother for fisrt derivatives

DCpVal = function(ind,lmat,Cmat,Svec,n,DeriEst,EmpiY){
     Tvec = rowSums(Cmat*(2*lmat-Cmat)*((matrix(rep(Svec,n),ncol=n,byrow = TRUE))^2))
     DCp = sum((EmpiY[ind] - DeriEst[ind])^2) + sum(Tvec[ind])
     print(DCp)
     return(DCp)
}


# G_D1Matrix = function(PtMat=PtMat,beta_nVec=beta_nVec,aVec=aVec,Yvec=Yvec,Xvec=Xvec,polyOd=polyOd,xvec=Xvec,Svec=Svec){
#   n=length(Yvec)
#   TempYmat = diag(n)
#   GMatD1=apply(TempYmat,1,CPEcompute_D1VecG,PtMat=PtMat,beta_nVec=beta_nVec,
#                aVec=aVec,Xvec=Xvec,polyOd=polyOd,xvec=Xvec,Svec=Svec)
#   return(GMatD1)
# }

lmat = LmatrixD1(x=Xvec,h=h,weights=1/Svec^2,polyOd=polyOd)
EmpiYCal = EmpiricalDri2(Xvec=Xvec,Svec=Svec,Yvec=Yvec,k=25)
EmpiY = EmpiYCal$EmpiY
Cmat = EmpiYCal$Cmat
DeriEst = predict(locfit.raw(x=lp(Xvec,h=h),y=Yvec,weights=1/Svec^2,deriv=1,kern="gauss",deg=polyOd),Xvec)









