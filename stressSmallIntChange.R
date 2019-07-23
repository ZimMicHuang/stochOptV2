runTest_smallInt = function(
  simCF, # simulated cash flows
  l_df, # simulated liabilities
  stoch_mod, # stochastic portfolio
  base_mod, # regulat portfolio
  int_mod, # model for interest rates. Must have simulate(.,n) defined for integer n
  RF, # risk-free rate used in optimization
  plotit=TRUE
) {
  require(xts)
  S = nrow(simCF[[1]])
  TF=ncol(simCF[[1]])
  N=length(simCF)
  CF_stoch = matrix(nrow=S,ncol=TF) 
  CF_base = matrix(nrow=S,ncol=TF) 
  l_mat = matrix(nrow=S,ncol=TF) 
  for (i in 1:S) {
    Cf = matrix(nrow=N,ncol=TF)
    for (j in 1:N) { # j is the index of asset
      Cf[j,] = simCF[[j]][i,] # Cf_i matrix's j-th row is i-th realizations of cash flows of j-th asset
    }
    simRates = simulate(int_mod,n=TF+1)/100 # percentage conversion!!!!!!
    simRates = simRates+(RF-simRates[1]) # make sure simulation started at discount rate used in model
    simRates = simRates[-1]/12 # converts to monthly
    CF_stoch[i,] = t(stoch_mod$soln)%*%Cf/cumprod(1+simRates) # discount -> PV's
    CF_base[i,] = t(base_mod$soln)%*%Cf/cumprod(1+simRates) 
    l_mat[i,] = l_df[i,]/cumprod(1+simRates)
  }
  CF_stoch = xts(t(CF_stoch),order.by = seq.Date(from=Sys.Date(),by="month",length=TF))
  CF_base = xts(t(CF_base),order.by = seq.Date(from=Sys.Date(),by="month",length=TF))
  l_mat = xts(t(l_mat),order.by = seq.Date(from=Sys.Date(),by="month",length=TF))
  stoch_NPV = apply(CF_stoch,2,sum)
  base_NPV = apply(CF_base,2,sum)
  l_NPV = apply(l_mat,2,sum)
  if (plotit) {
    h1 = ecdf(base_NPV)
    h2 = ecdf(stoch_NPV)
    h3 = ecdf(l_NPV)
    plot(h1,main="Distribution of NPV of Cash Flows",
         xlim=c(min(base_NPV,stoch_NPV,l_NPV),1.5*max(base_NPV,stoch_NPV)),
         lwd=2,col="green",do.points=FALSE)
    plot(h2,col="blue",add=TRUE,do.points=FALSE,lwd=2)
    plot(h3,col="red",add=TRUE,lwd=2,do.points=FALSE)
    legend("topleft",legend=c("random","non-random","liability"),col=c("blue","green","red"),pch=16)
  }
  return(list(stoch_NPV=stoch_NPV,base_NPV=base_NPV,l_NPV=l_NPV))
}