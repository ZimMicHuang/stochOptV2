stochLiabMatch = function(
  l_df, # simulated liability 
  simCF, # simulated cash flows
  markDat, # market data
  EPS = 0.005, # relative to the mean of liabilities' NPV
  plotit=TRUE,
  n_iter=100 # max iteration
) {
  discount = markDat$RF[1]
  N=length(simCF)
  S=nrow(simCF[[1]])
  TF=ncol(simCF[[1]])
  NPV = numeric(TF) #NPV: vector of cash flow --> net present value
  for (i in 1:TF) {
    NPV[i] = 1/(1+discount/12)^i
  }
  NPV = as.matrix(NPV,ncol=1) 
  
  tau = matrix(nrow=TF,ncol=1) #tau: vector of cash flow --> numerator of duration
  for (i in 1:TF) {
    tau[i,1] = i/(1+discount/12)^i
  }
  
  l = sort(as.numeric(as.matrix(l_df)%*%NPV)) # scenarios--NPV of Liability
  l_t = as.numeric((apply(l_df,2,sum)/S)%*%tau) # Expected duration of liability
  # block matrix algebra
  # block matrix for stoch. dom. constraint
  M = matrix(nrow=S,ncol=N) # M: u -> R^S, map to realizations of cash flow
  # block matrix for immunization constraint
  W = matrix(nrow=S,ncol=N)
  
  for (i in 1:S) {
    Cf = matrix(nrow=N,ncol=TF)
    for (j in 1:N) { 
      Cf[j,] = simCF[[j]][i,] 
    }
    M[i,]=t(Cf%*%NPV)
    W[i,]=t(Cf%*%tau)
  }
  e=matrix(rep(1/S,times=S),ncol=1)
  W = t(e)%*%W
  # we have that M%*%u=empirical distribution of NPV; W%*%u=sample mean duration
  
  
  EPS = 0.01
  findCut = function(u,liab=l) {
    # given the decision variable, 
    # find the quantile which NPV cash flow >(2) NPV liability is violated, 
    # and find auxillary vector e_lam such that the index set lambda comtains all the scenarios 
    # which violates the stoch, dom. constraint
    # input: u: decision var
    # output: e_lam: index row vector, 1*S 
    #         lq: sum of l at level q
    
    # u is the PRINCIPE VALUE/FACE VALUE in each bond
    
    dist_xnpv = M%*%u
    max_diff = 0
    res=list()
    for (i in 2:S) {
      # i/S is the probability level
      q = liab[i] 
      tmp = (dist_xnpv<=q) # index of cut
      if (sum(tmp)==0) {
        FX = -1 # no need to use -Inf because l,x always non-negative
      } else {
        FX = mean(q-dist_xnpv[which(tmp)]) 
      }
      FL = mean(q-liab[which(liab<=q)])
      if (FX-FL>=max_diff) {
        max_diff=FX-FL
        res = list(e_lam=tmp/sum(tmp),lq=mean(liab[which(liab<=q)]))
      }
      # if (FX>FL+EPS) {
      #   #print(i)
      #   #break
      #   return(list(e_lam=tmp/sum(tmp),lq=mean(l[which(l<=q)])))
      # }
    }
    if (max_diff<=EPS) {
      print("Stochastic Constraint Satified! Returning -1.")
      return(list(e_lam=-1))
    } else {
      return(res)
    }
  }  
  constr = function(I,B) {
    # given e_lam's stored in matrix I, where each row is an e_lam,
    # the index vector for which scenarios violated the stoch. dom. constr.,
    # get the set of constraints to be satisfied in standard form.
    # Ax>=b
    # input: I: matrix storing sets of cuts that the stochastic dominance constraint is violated
    #         B: vector storing the expected shortfalls of the liabilities
    # output: A: matrix, #constr*N
    #         b: vector, #constr
    A=matrix(nrow=nrow(I),ncol=N)
    b=numeric(nrow(I))
    for (i in 1:nrow(I)) {
      e_lam = I[i,]
      A[i,] = e_lam%*%M
      b[i] = B[i]
    }
    return(list(A=A,b=b))
  }
  
  # iterative cutting method with simplex
  require(boot)
  I = matrix(rep(1/S,times=S),nrow=1)
  B = c(mean(l))
  col=terrain.colors(n=n_iter)
  for (i in 1:n_iter) {
    #print(i)
    inequ = constr(I,B)
    sol = boot::simplex(a=markDat$mp,A1=-1*markDat$mp,b1=0,
                        A2=inequ$A,
                        b2=inequ$b,
                        A3=W,
                        b3=l_t)
    tmp = findCut(sol$soln)
    e_lam = tmp$e_lam
    lq = tmp$lq
    if (plotit) {
      if (i ==1) {
        plot(ecdf(M %*% sol$soln),xlim=c(0,4000),ylim=c(0,1),lwd=1,
             main="NPV of Portfolio Cash-Flow vs Liability",col="pink",do.points=FALSE)
        plot(ecdf(l),add=TRUE,lwd=3,do.points=FALSE)
      } else if (i %% 1 == 0) {
        plot(ecdf(M %*% sol$soln),add=TRUE,col="pink",do.points=FALSE,lwd=1)
      }
    }
    if (e_lam==-1) {
      break
    } else {
      I = rbind(I,t(e_lam))
      B = c(B,lq)
    }
  }
  if (plotit) {
    legend("right",col=c("red","black"),legend = c("Portfolio", "Liability"),pch=15 )
    plot(ecdf(M %*% sol$soln),add=TRUE, col="red",do.points=FALSE,lwd=3)
  }
  
  stoch_mod = sol
  names(stoch_mod$soln)=names(simCF)
  stoch_mod$b2 = inequ$b
  return(stoch_mod)
}  


regLiabMatch = function(
  l_df, # simulated liability 
  simCF, # simulated cash flows
  markDat # market data
) {
  require(boot)
  discount = markDat$RF[1]
  N=length(simCF)
  S=nrow(simCF[[1]])
  TF=ncol(simCF[[1]])
  NPV = numeric(TF) #NPV: vector of cash flow --> net present value
  for (i in 1:TF) {
    NPV[i] = 1/(1+discount/12)^i
  }
  NPV = as.matrix(NPV,ncol=1) 
  tau = matrix(nrow=TF,ncol=1) #tau: vector of cash flow --> numerator of duration
  for (i in 1:TF) {
    tau[i,1] = i/(1+discount/12)^i
  }
  l = sort(as.numeric(as.matrix(l_df)%*%NPV)) # scenarios--NPV of Liability
  l_t = as.numeric((apply(l_df,2,sum)/S)%*%tau) # Expected duration of liability
  
  # block matrix for stoch. dom. constraint
  M = matrix(nrow=1,ncol=N) # M: u -> R^S, map to realizations of cash flow
  # block matrix for immunization constraint
  W = matrix(nrow=1,ncol=N)
  Cf = matrix(nrow=N,ncol=TF)
  for (j in 1:N) { # j is the index of asset
    Cf[j,] = c(
      rep(markDat$coupon[j],times=markDat$nper[j]),
      rep(0,TF-markDat$nper[j])
    )
    Cf[j,markDat$nper[j]] = 1+Cf[j,markDat$nper[j]]
  }
  M=t(Cf%*%NPV)
  W=t(Cf%*%tau)
  
  sol = boot::simplex(a=markDat$mp,A1=rep(0,times=N),b1=1,
                      A2=M,
                      b2=mean(l),
                      A3=W,
                      b3=l_t)
  base_mod=sol
  names(base_mod$soln)=names(simCF)
  base_mod$b2 = mean(l)
  return(base_mod)
}


# # plot stochastic dominance relation
# 
# cf_npv = sort(M%*%stoch_mod$soln)
# l_npv = sort(l)
# z=1500
# s = 1000
# trapezoid = function(f,a,b,nInt=500) {
#   # composite trapezoid for numerical integration
#   h = (b-a)/nInt
#   res = 0 
#   for (i in 1:nInt) {
#     res = res+1/2*(f(a)+f(a+h))*h
#     a=a+h
#   }
#   return(res)
# }
# FX = function(s,X=cf_npv) {
#   return(mean(X<=s))
# }
# FL = function(s,L=l_npv) {
#   return(mean(L<=s))
# }
# 
# DZ = function(z) {
#   return(trapezoid(f=function(s) {FL(s)-FX(s)},a=0,b=z))
# }
# 
# plot(seq(1,1000,by=10),sapply(seq(1,1000,by=10),DZ),type="l")
# 
