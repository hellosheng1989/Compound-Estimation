CPEcpCompute = function(CPEvec,Gvec,Svec=Svec){
  Cp = 1/n*sum((Yvec-CPEvec)^2) + 2/n*sum((Svec^2)*Gvec)
  print(Cp)
  return(Cp)
}
