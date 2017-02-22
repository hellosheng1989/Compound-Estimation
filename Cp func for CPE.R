# x is a number, aVec and beta_nVec should be the same length
Get_Wvec = function(aVec, x, beta_nVec) {
    Denominate = sum(exp(-beta_nVec * (x - aVec)^2))
    Wvec = exp(-beta_nVec * (x - aVec)^2)/Denominate
    return(Wvec)
}

# derivative of Weight W

Get_Wderi = function(aVec, x, beta_nVec) {
    A = sum(exp(-beta_nVec * (x - aVec)^2))
    B = sum(beta_nVec * (x - aVec) * exp(-beta_nVec * (x - aVec)^2))
    Wvec_Deri = (-2 * beta_nVec * (x - aVec) * exp(-beta_nVec * (x - aVec)^2) * A + 2 * exp(-beta_nVec * 
        (x - aVec)^2) * B)/A^2
    return(Wvec_Deri)
}


## get linear smoother for X_i at point X_i : G_ni(X_i,x=X_i)
Get_G = function(i = 1, aVec = aVec, Xvec = Xvec, Yvec = Yvec, h = 0.5, polyOd = 5, beta_nVec = beta_nVec, 
    Svec = Svec) {
    L_n = length(aVec)
    n = length(Xvec)
    J = polyOd - 1
    lMat = matrix(rep(0, L_n * (J + 1)), ncol = J + 1)
    polyMat = matrix(rep(0, (J + 1) * L_n), ncol = L_n)
    tempy = rep(0, n)
    tempy[i] = 1
    for (j in 0:J) {
        lMat[, j + 1] = predict(locfit.raw(Xvec, tempy, alpha = c(0, h), deriv = rep(1, j), weights = 1/Svec^2, 
            deg = polyOd, kern = "gauss"), aVec)
        polyMat[j + 1, ] = (Xvec[i] - aVec)^j
    }
    polynoVec = diag(lMat %*% polyMat)
    Wvec = Get_Wvec(aVec, x = Xvec[i], beta_nVec)
    tempG = sum(Wvec * polynoVec)
    return(tempG)
}

## get linear smoother for X_i at point X_k,: G_n_i_k(X_k,x=x_i)
Get_G2 = function(k = 1, i = 1, aVec = aVec, Xvec = Xvec, Yvec = Yvec, h = 0.5, polyOd = 5, beta_nVec = beta_nVec, 
    Svec = Svec) {
    L_n = length(aVec)
    n = length(Xvec)
    J = polyOd - 1
    lMat = matrix(rep(0, L_n * (J + 1)), ncol = J + 1)
    polyMat = matrix(rep(0, (J + 1) * L_n), ncol = L_n)
    tempy = rep(0, n)
    tempy[k] = 1
    for (j in 0:J) {
        lMat[, j + 1] = predict(locfit.raw(Xvec, tempy, alpha = c(0, h), deriv = rep(1, j), weights = 1/Svec^2, 
            deg = polyOd, kern = "gauss"), aVec)
        polyMat[j + 1, ] = (Xvec[i] - aVec)^j
    }
    polynoVec = diag(lMat %*% polyMat)
    Wvec = Get_Wvec(aVec, x = Xvec[i], beta_nVec)
    tempG = sum(Wvec * polynoVec)
    return(tempG)
}

## Get linear smoother for X_i at point X_i, for the derivative estimation
Get_GD = function(i = 1, aVec = aVec, Xvec = Xvec, Yvec = Yvec, h = 0.5, polyOd = 5, beta_nVec = beta_nVec, 
    Svec = Svec) {
    L_n = length(aVec)
    n = length(Xvec)
    J = polyOd - 1
    lMat = matrix(rep(0, L_n * (J + 1)), ncol = J + 1)
    polyMat = matrix(rep(0, (J + 1) * L_n), ncol = L_n)
    polyMat2 = matrix(rep(0, (J + 1) * L_n), ncol = L_n)
    tempy = rep(0, n)
    tempy[i] = 1
    for (j in 0:J) {
        lMat[, j + 1] = predict(locfit.raw(Xvec, tempy, alpha = c(0, h), deriv = rep(1, j), weights = 1/Svec^2, 
            deg = polyOd, kern = "gauss"), aVec)
        polyMat[j + 1, ] = (Xvec[i] - aVec)^j
        polyMat2[j + 1, ] = j * (Xvec[i] - aVec)^(j - 1)
    }
    polynoVec = diag(lMat %*% polyMat)
    polynoVec2 = diag(lMat %*% polyMat2)
    Wvec = Get_Wvec(aVec, x = Xvec[i], beta_nVec)
    WvecD = Get_Wderi(aVec, x = Xvec[i], beta_nVec)
    tempG = sum(WvecD * polynoVec + Wvec * polynoVec2)
    return(tempG)
}




# calculate compound estimation with same bandwidth
CPEcompute = function(h = 0.5, beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, 
    x = x, Svec = Svec) {
    L_n = length(aVec)
    J = polyOd - 1
    cMat = matrix(rep(0, L_n * (J + 1)), ncol = J + 1)
    polyMat = matrix(rep(0, (J + 1) * L_n), ncol = L_n)
    for (j in 0:J) {
        cMat[, j + 1] = predict(locfit.raw(Xvec, Yvec, alpha = c(0, h), deriv = rep(1, j), weights = 1/Svec^2, 
            deg = polyOd, kern = "gauss"), aVec)
        polyMat[j + 1, ] = (x - aVec)^j
    }
    Wvec = Get_Wvec(aVec = aVec, x = x, beta_nVec = beta_nVec)
    polynoVec = diag(cMat %*% polyMat)  # \sum_{j=0}^J c_{j,a}(x-a)^j
    CPE = sum(Wvec * polynoVec)
    return(CPE)
}


## first derivative - compound estimation with same bandwidth
CPEcompute_D1 = function(h = 0.5, beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, 
    x = x, Svec = Svec) {
    L_n = length(aVec)
    J = polyOd - 1
    cMat = matrix(rep(0, L_n * (J + 1)), ncol = J + 1)
    polyMat = matrix(rep(0, (J + 1) * L_n), ncol = L_n)
    for (j in 0:J) {
        cMat[, j + 1] = predict(locfit.raw(Xvec, Yvec, alpha = c(0, h), deriv = rep(1, j), weights = 1/Svec^2, 
            deg = polyOd, kern = "gauss"), aVec)
        polyMat[j + 1, ] = (x - aVec)^j
    }
    cMat1 = matrix(rep(0, L_n * J), ncol = J)
    polyMat1 = matrix(rep(0, J * L_n), ncol = L_n)
    for (j in 1:J) {
        cMat1[, j] = predict(locfit.raw(Xvec, Yvec, alpha = c(0, h), deriv = rep(1, j), weights = 1/Svec^2, 
            deg = polyOd, kern = "gauss"), aVec)
        polyMat1[j, ] = j * (x - aVec)^(j - 1)
    }
    polynoVec = diag(cMat %*% polyMat)  # \sum_{j=0}^J c_{j,a}(x-a)^j
    polynoVecD = diag(cMat1 %*% polyMat1)  # \sum_{j=1}^J c_{j,a} j(x-a)^(j-1)
    Wvec = Get_Wvec(aVec = aVec, x = x, beta_nVec = beta_nVec)
    WvecD = Get_Wderi(aVec = aVec, x = x, beta_nVec = beta_nVec)
    CPE_D1 = sum(Wvec * polynoVecD + WvecD * polynoVec)
    return(CPE_D1)
}



# calculate compound estimation with different bandwidth



#### choose common h and common beta
CPEcpCompute = function(par = c(ADh_par, ADbeta_n), aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = 3, Svec = Svec) {
    ADh_par = par[1]
    ADbeta_n = par[2]
    h_par = exp(ADh_par)/(1 + exp(ADh_par))
    beta_n = exp(ADbeta_n)
    beta_nVec = rep(beta_n, length(aVec))
    n = length(Xvec)
    CPEvec = rep(0, n)
    Gvec = rep(0, n)
    for (i in 1:n) {
        CPEvec[i] = CPEcompute(h = h_par, beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, 
            x = Xvec[i], Svec = Svec)
        Gvec[i] = Get_G(i = i, aVec = aVec, Xvec = Xvec, Yvec = Yvec, h = h_par, polyOd = polyOd, beta_nVec = beta_nVec, 
            Svec = Svec)
    }
    Cp = 1/n * sum((Yvec - CPEvec)^2) + 2/n * sum((Svec^2) * Gvec)
    print(Cp)
    return(Cp)
}


# choose common h and different beta! par = c(ad_h,ad_betaVec)
CPEcpCompute1 = function(par = c(ADh_par, ADbeta_nVec), aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = 3, 
    Svec = Svec) {
    n = length(Xvec)
    CPEvec = rep(0, n)
    Gvec = rep(0, n)
    ADh_par = par[1]
    h_par = exp(ADh_par)/(1 + exp(ADh_par))
    ADbeta_nVec = par[-1]
    beta_nVec = exp(ADbeta_nVec)
    for (i in 1:n) {
        CPEvec[i] = CPEcompute(h = h_par, beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, 
            x = Xvec[i], Svec = Svec)
        Gvec[i] = Get_G(i = i, aVec = aVec, Xvec = Xvec, Yvec = Yvec, h = h_par, polyOd = polyOd, beta_nVec = beta_nVec, 
            Svec = Svec)
    }
    Cp = 1/n * sum((Yvec - CPEvec)^2) + 2/n * sum((Svec^2) * Gvec)
    print(Cp)
    return(Cp)
}

## choose beta when beta is reparameterized by log(beta) = r0 + r1*a + r2*a^2...
CPEcpComputeRP = function(par = c(ADh_par, g0, g1, g2), aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = 3, 
    Svec = Svec) {
    n = length(Xvec)
    CPEvec = rep(0, n)
    Gvec = rep(0, n)
    ADh_par = par[1]
    h_par = exp(ADh_par)/(1 + exp(ADh_par))
    r0 = par[2]
    r1 = par[3]
    r2 = par[4]
    ADbeta_nVec = r0 + r1 * aVec + r2 * aVec^2
    beta_nVec = exp(ADbeta_nVec)
    for (i in 1:n) {
        CPEvec[i] = CPEcompute(h = h_par, beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, 
            x = Xvec[i], Svec = Svec)
        Gvec[i] = Get_G(i = i, aVec = aVec, Xvec = Xvec, Yvec = Yvec, h = h_par, polyOd = polyOd, beta_nVec = beta_nVec, 
            Svec = Svec)
    }
    Cp = 1/n * sum((Yvec - CPEvec)^2) + 2/n * sum((Svec^2) * Gvec)
    print(Cp)
    return(Cp)
}


## choose beta when beta is reparameterized by log(beta) = C - 2*log(Svec(a))
CPEcpComputeRPV = function(par = c(ADh_par, C), aVec = aVec, Svec.a = Svec.a, Yvec = Yvec, Xvec = Xvec, polyOd = 3, 
    Svec = Svec) {
    n = length(Xvec)
    CPEvec = rep(0, n)
    Gvec = rep(0, n)
    ADh_par = par[1]
    h_par = exp(ADh_par)/(1 + exp(ADh_par))
    C = par[2]
    ADbeta_nVec = C + 2 * log(Svec.a)
    beta_nVec = exp(ADbeta_nVec)
    for (i in 1:n) {
        CPEvec[i] = CPEcompute(h = h_par, beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, 
            x = Xvec[i], Svec = Svec)
        Gvec[i] = Get_G(i = i, aVec = aVec, Xvec = Xvec, Yvec = Yvec, h = h_par, polyOd = polyOd, beta_nVec = beta_nVec, 
            Svec = Svec)
    }
    Cp = 1/n * sum((Yvec - CPEvec)^2) + 2/n * sum((Svec^2) * Gvec)
    print(Cp)
    return(Cp)
}







## choose beta when fixed h
CPEcpComputeF1 = function(par = ADbeta_n, h_par = h_par, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = 3, 
    Svec = Svec) {
    ADbeta_n = par
    beta_n = exp(ADbeta_n)
    beta_nVec = rep(beta_n, length(aVec))
    n = length(Xvec)
    CPEvec = rep(0, n)
    Gvec = rep(0, n)
    for (i in 1:n) {
        CPEvec[i] = CPEcompute(h = h_par, beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, 
            x = Xvec[i], Svec = Svec)
        Gvec[i] = Get_G(i = i, aVec = aVec, Xvec = Xvec, Yvec = Yvec, h = h_par, polyOd = polyOd, beta_nVec = beta_nVec, 
            Svec = Svec)
    }
    Cp = 1/n * sum((Yvec - CPEvec)^2) + 2/n * sum((Svec^2) * Gvec)
    print(Cp)
    return(Cp)
}



## choose different beta when fixed h
CPEcpComputeF2 = function(par = ADbeta_nVec, h_par = h_par, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = 3, 
    Svec = Svec) {
    n = length(Xvec)
    CPEvec = rep(0, n)
    Gvec = rep(0, n)
    h_par = h_par
    ADbeta_nVec = par
    beta_nVec = exp(ADbeta_nVec)
    for (i in 1:n) {
        CPEvec[i] = CPEcompute(h = h_par, beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, 
            x = Xvec[i], Svec = Svec)
        Gvec[i] = Get_G(i = i, aVec = aVec, Xvec = Xvec, Yvec = Yvec, h = h_par, polyOd = polyOd, beta_nVec = beta_nVec, 
            Svec = Svec)
    }
    Cp = 1/n * sum((Yvec - CPEvec)^2) + 2/n * sum((Svec^2) * Gvec)
    print(Cp)
    return(Cp)
}

## choose beta when reparameterized, h is fixed
CPEcpComputeRP1 = function(par = c(r0, r1, r2), h_par = h_par, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = 3, 
    Svec = Svec) {
    n = length(Xvec)
    CPEvec = rep(0, n)
    Gvec = rep(0, n)
    h_par = h_par
    r0 = par[1]
    r1 = par[2]
    r2 = par[3]
    ADbeta_nVec = r0 + r1 * aVec + r2 * aVec^2
    beta_nVec = exp(ADbeta_nVec)
    for (i in 1:n) {
        CPEvec[i] = CPEcompute(h = h_par, beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, 
            x = Xvec[i], Svec = Svec)
        Gvec[i] = Get_G(i = i, aVec = aVec, Xvec = Xvec, Yvec = Yvec, h = h_par, polyOd = polyOd, beta_nVec = beta_nVec, 
            Svec = Svec)
    }
    Cp = 1/n * sum((Yvec - CPEvec)^2) + 2/n * sum((Svec^2) * Gvec)
    print(Cp)
    return(Cp)
}


## choose beta when beta is reparameterized by log(beta) = C - 2*log(Svec(a)),h is fixed
CPEcpComputeRPV1 = function(par = C, aVec = aVec, h_par = h_par, Svec.a = Svec.a, Yvec = Yvec, Xvec = Xvec, 
    polyOd = 3, Svec = Svec) {
    n = length(Xvec)
    CPEvec = rep(0, n)
    Gvec = rep(0, n)
    ADbeta_nVec = C + 2 * log(Svec.a)
    beta_nVec = exp(ADbeta_nVec)
    for (i in 1:n) {
        CPEvec[i] = CPEcompute(h = h_par, beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, 
            x = Xvec[i], Svec = Svec)
        Gvec[i] = Get_G(i = i, aVec = aVec, Xvec = Xvec, Yvec = Yvec, h = h_par, polyOd = polyOd, beta_nVec = beta_nVec, 
            Svec = Svec)
    }
    Cp = 1/n * sum((Yvec - CPEvec)^2) + 2/n * sum((Svec^2) * Gvec)
    print(Cp)
    return(Cp)
}

######################################################################################################## 






######################################################################################################## get linear smoother for X_i at point X_i , when with different bandwidth h
Get_G1 = function(i = 1, aVec = aVec, Xvec = Xvec, Yvec = Yvec, polyOd = 5, beta_nVec = beta_nVec, Svec = Svec) {
    L_n = length(aVec)
    n = length(Xvec)
    J = polyOd - 1
    lMat = matrix(rep(0, L_n * (J + 1)), ncol = J + 1)
    polyMat = matrix(rep(0, (J + 1) * L_n), ncol = L_n)
    tempy = rep(0, n)
    tempy[i] = 1
    for (j in 0:J) {
        lMat[, j + 1] = predict(locfit(tempy ~ Xvec, alpha = c(0, 0, 2), deriv = rep(1, j), acri = "cp", 
            weights = 1/Svec^2, deg = polyOd, kern = "gauss"), aVec)
        polyMat[j + 1, ] = (Xvec[i] - aVec)^j
    }
    polynoVec = diag(lMat %*% polyMat)
    Wvec = Get_Wvec(aVec, x = Xvec[i], beta_nVec)
    tempG = sum(Wvec * polynoVec)
    return(tempG)
}

### Compute compound estimation with different h from local Cp in the first step
CPEcompute2 = function(beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, x = x, 
    Svec = Svec) {
    L_n = length(aVec)
    J = polyOd - 1
    cMat = matrix(rep(0, L_n * (J + 1)), ncol = J + 1)
    polyMat = matrix(rep(0, (J + 1) * L_n), ncol = L_n)
    for (j in 0:J) {
        cMat[, j + 1] = predict(locfit(Yvec ~ Xvec, alpha = c(0, 0, 2), deriv = rep(1, j), acri = "cp", weights = 1/Svec^2, 
            deg = polyOd, kern = "gauss"), aVec)
        polyMat[j + 1, ] = (x - aVec)^j
    }
    Wvec = Get_Wvec(aVec = aVec, x = x, beta_nVec = beta_nVec)
    polynoVec = diag(cMat %*% polyMat)  # \sum_{j=0}^J c_{j,a}(x-a)^j
    CPE = sum(Wvec * polynoVec)
    return(CPE)
}

### Calculate the first derivative for above






### choose common beta,input beta a scale here
CPEcpComputeDHCB = function(par = ADbeta_n, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = 3, Svec = Svec) {
    ADbeta_n = par
    beta_n = exp(ADbeta_n)
    beta_nVec = rep(beta_n, length(aVec))
    n = length(Xvec)
    CPEvec = rep(0, n)
    Gvec = rep(0, n)
    for (i in 1:n) {
        CPEvec[i] = CPEcompute2(beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, 
            x = Xvec[i], Svec = Svec)
        Gvec[i] = Get_G1(i = i, aVec = aVec, Xvec = Xvec, Yvec = Yvec, polyOd = polyOd, beta_nVec = beta_nVec, 
            Svec = Svec)
    }
    Cp = 1/n * sum((Yvec - CPEvec)^2) + 2/n * sum((Svec^2) * Gvec)
    print(Cp)
    return(Cp)
}


# choose different beta! input beta is a vector here
CPEcpComputeDHDB = function(par = ADbeta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = 3, Svec = Svec) {
    n = length(Xvec)
    CPEvec = rep(0, n)
    Gvec = rep(0, n)
    ADbeta_nVec = par
    beta_nVec = exp(ADbeta_nVec)
    for (i in 1:n) {
        CPEvec[i] = CPEcompute2(beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, 
            x = Xvec[i], Svec = Svec)
        Gvec[i] = Get_G1(i = i, aVec = aVec, Xvec = Xvec, Yvec = Yvec, polyOd = polyOd, beta_nVec = beta_nVec, 
            Svec = Svec)
    }
    Cp = 1/n * sum((Yvec - CPEvec)^2) + 2/n * sum((Svec^2) * Gvec)
    print(Cp)
    return(Cp)
}




######################## NOT U

#### choose different h and common beta # not finished
CPEcpCompute2 = function(par = c(ADh_parVec, ADbeta_n), aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = 3, 
    Svec = Svec) {
    Ln = length(aVec)
    ADh_par = par[-(Ln + 1)]
    ADbeta_n = par[Ln + 1]
    h_parVec = exp(ADh_par)/(1 + exp(ADh_par))
    beta_n = exp(ADbeta_n)
    beta_nVec = rep(beta_n, length(aVec))
    n = length(Xvec)
    CPEvec = rep(0, n)
    Gvec = rep(0, n)
    for (i in 1:n) {
        CPEvec[i] = CPEcompute(h = h_par, beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, 
            x = Xvec[i], Svec = Svec)
        Gvec[i] = Get_G(i = i, aVec = aVec, Xvec = Xvec, Yvec = Yvec, h = h_par, polyOd = polyOd, beta_nVec = beta_nVec, 
            Svec = Svec)
    }
    Cp = 1/n * sum((Yvec - CPEvec)^2) + 2/n * sum((Svec^2) * Gvec)
    print(Cp)
    return(Cp)
}

# choose different h and differen beta! par = c(ad_h,ad_betaVec) ##not finished
CPEcpCompute3 = function(par = c(ADh_parVec, ADbeta_nVec), aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = 3, 
    Svec = Svec) {
    n = length(Xvec)
    CPEvec = rep(0, n)
    Gvec = rep(0, n)
    ADh_par = par[1]
    h_par = exp(ADh_par)/(1 + exp(ADh_par))
    ADbeta_nVec = par[-1]
    beta_nVec = exp(ADbeta_nVec)
    for (i in 1:n) {
        CPEvec[i] = CPEcompute(h = h_par, beta_nVec = beta_nVec, aVec = aVec, Yvec = Yvec, Xvec = Xvec, polyOd = polyOd, 
            x = Xvec[i], Svec = Svec)
        Gvec[i] = Get_G(i = i, aVec = aVec, Xvec = Xvec, Yvec = Yvec, h = h_par, polyOd = polyOd, beta_nVec = beta_nVec, 
            Svec = Svec)
    }
    Cp = 1/n * sum((Yvec - CPEvec)^2) + 2/n * sum((Svec^2) * Gvec)
    print(Cp)
    return(Cp)
}




