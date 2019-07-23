getDefDat = function(type=c("read","gen"), # read from local CSV or generate data
                       ANNUAL_DEF_UB=0.4, # highest annual default prob
                       ANNUAL_DEF_LB=0.03, # lowest default rate among basket
                       N_COMP = 100,  # how many companies to consider
                       f="defProb.csv") {
  if (type=="gen") {
    df0 = numeric(N_COMP)
    df0[1:N_COMP] = runif(n=N_COMP,ANNUAL_DEF_LB,ANNUAL_DEF_UB)
    names(df0)=sapply(1:N_COMP,FUN=function(x) {paste0("CORP",x)})
    return(df0)
  } else if (type=="read") {
    df1 = read.csv(f)
    df0 = df1[,-1]
    names(df0) = df1[,1]
    return(df0)
  }
}

runSimDef = function(df0,N_SIM=5000,save=FALSE) {
  # df: numeric vector containing all relevant info (Time horizon and default prob)
  # returns: simDef: list of simulated times of default
  N_COMP = length(df0)
  defProb = df0[1:N_COMP]
  simDef = vector("list",length=N_COMP)
  for (i in 1:N_COMP) {
    p = 1-(1-defProb[i])^(1/12) # convert to monthly default prob
    simDef[[i]] = rgeom(n=N_SIM,prob=p)+1
  }
  names(simDef)=names(df0)
  if (save) {
    save(simDef,"simulated default times.Rdata")
  }
  return(simDef)
}



