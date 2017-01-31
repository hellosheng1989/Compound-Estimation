#calulate hat matrix L, order=7
Lmatrix<-function(x,h,weights=weight,polyOd=5)
{
  n<-length(x)
  Lmatrix<-matrix(0,n,n)
  for (k in 1:n)
  {
    tempy <- rep(0,n)
    tempy[k] <- 1
    Lmatrix[,k] <- predict(locfit.raw(x,tempy,alpha=c(0,h),deg=polyOd,weights=weights,kern="gauss"),x)
  }
  return(Lmatrix)
}

#calulate hat matrix L for first derivative estimation
LmatrixD1<-function(x,h,weights=weight,polyOd=5)
{
  n<-length(x)
  Lmatrix<-matrix(0,n,n)
  for (k in 1:n)
  {
    tempy <- rep(0,n)
    tempy[k] <- 1
    Lmatrix[,k] <- predict(locfit.raw(x,tempy,alpha=c(0,h),deriv=1,deg=polyOd,weights=weights,kern="gauss"),x)
  }
  return(Lmatrix)
}

#get LP Estimator and weight diagram only at each design points, corresponding to G(X_i,h)
Local_Gvec = function(Yvec,Xvec,h_par, Svec=(abs(1-Xvec)),polyOd=5){
  In_weights=1/Svec^2
  locObj = locfit.raw(lp(Xvec,h=h_par),y=Yvec,weights=In_weights,deg=polyOd,kern="gauss")
  Evec = predict(locObj,Xvec)
  W_Diagram = Lmatrix(x=Xvec,h=h_par,weights=In_weights,polyOd=polyOd)
  Gvec = diag(W_Diagram)
  return(list(Evec=Evec,Gvec=Gvec))
}


#get LP first derivative Estimator and weight diagram W(Xi,Xj)
Local_WdigD1 = function(Yvec,Xvec,h_par, Svec=(abs(1-Xvec)),polyOd=5){
  In_weights=1/Svec^2
  locObj = locfit.raw(lp(Xvec,h=h_par),y=Yvec,weights=In_weights,deriv=1,deg=polyOd,kern="gauss")
  ED1vec = predict(locObj,Xvec)
  WD1_Diagram = LmatrixD1(x=Xvec,h=h_par,weights=In_weights,polyOd=polyOd)
  return(list(ED1vec=ED1vec,WD1_Diagram=WD1_Diagram))
}

####Cp for local regression
LpCpCompute = function(ADh_par=0,Yvec=Yvec,Xvec=Xvec,polyOd=polyOd,Svec=Svec){
  n=length(Xvec); h_par = exp(ADh_par)/(1+exp(ADh_par))
  TempResults = Local_Gvec(Yvec=Yvec,Xvec=Xvec,h_par=h_par,Svec=Svec,polyOd=polyOd)
  Evec = TempResults$Evec
  Gvec = TempResults$Gvec
  Cp = 1/n*sum((Yvec-Evec)^2) + 2/n*sum((Svec^2)*Gvec)
  print(Cp)
  return(Cp)
}




#Gcp for local regression of first derivative

