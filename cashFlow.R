runSimCF = function(
  markDat, # bond specs
  defProb=c(), # list of default probabilities
  defTimes=c(), # simulated default times 
  doneComp=c(), # list of already default and already paid bonds
  N_SIM=5000,
  TF=36
) {
  if (!is.null(defProb)) {
    # simulate using default probabilities
    N = length(defProb)
    S = N_SIM
    simDef = runSimDef(defProb,save=FALSE,N_SIM=S)
    simCF = vector("list",length=N)
    for (i in 1:N) {
      mat = matrix(nrow=S,ncol=TF)
      for (j in 1:S) { 
        timeOfDefault = simDef[[i]][j]
        nper = markDat$nper[i]
        if (timeOfDefault>nper) { # case company successful repays
          mat[j,] = c(rep(markDat$coupon[i],times=nper),rep(0,times=TF-nper))
          mat[j,nper] = mat[j,nper]+1 # principle payment 
        } else {
          # at some point the bond defaults
          mat[j,] = c(rep(markDat$coupon[i],times=timeOfDefault),rep(0,times=TF-timeOfDefault))
          mat[j,timeOfDefault] = mat[j,timeOfDefault]+markDat$rec[i]
        }
      }
      simCF[[i]] = mat
    }
    names(simCF) = names(simDef)
  } else if (!is.null(defTimes)) {
    # simulate using simulated default times
    N = length(defTimes)
    S = length(defTimes[[1]])
    simDef = defTimes
    simCF = vector("list",length=N)
    for (i in 1:N) {
      mat = matrix(nrow=S,ncol=TF)
      for (j in 1:S) { 
        timeOfDefault = simDef[[i]][j]
        nper = markDat$nper[i]
        if (timeOfDefault>nper) { # case company successful repays
          mat[j,] = c(rep(markDat$coupon[i],times=nper),rep(0,times=TF-nper))
          mat[j,nper] = mat[j,nper]+1 # principle payment 
        } else {
          # at some point the bond defaults
          mat[j,] = c(rep(markDat$coupon[i],times=timeOfDefault),rep(0,times=TF-timeOfDefault))
          mat[j,timeOfDefault] = mat[j,timeOfDefault]+markDat$rec[i]
        }
      }
      simCF[[i]] = mat
    }
    names(simCF) = names(simDef)
  } else {
    print("Insufficient data provided, returning -1.")
    return(-1)
  }
  
  if (!is.null(doneComp)) {
    # if a company is done, then all cashflows remain at 0
    for (i in 1:N) {
      if (names(simCF)[i] %in% doneComp) {
        simCF[[i]] = matrix(0,nrow=S,ncol=TF)
      }
    }
  }
  return(simCF)
}


