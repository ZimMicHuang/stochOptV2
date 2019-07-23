runSimLiab = function(
  type=c("read","gen"), 
  f = "liab.csv",
  MU_liab = 2, # parameters for the liability streams
  SD_liab = 0.7,
  N_SIM=5000, # number of scenarios
  TF=36, # time periods to simulate
  distribution = c("invLog","norm") # inverse log normal or normal
) {
  S = N_SIM
  if (type=="gen") {
    l_df = matrix(nrow=S,ncol=TF)
    if (distribution=="invLog") {
      for (i in 1:S) {
        l_s = exp(rnorm(TF+1,mean = MU_liab,sd=SD_liab))
        l_s = (max(l_s)-l_s)
        l_s = l_s[-which(l_s==0)]
        #l_s = rnorm(TF,mean=MU_liab,sd=SD_liab)
        l_df[i,] = l_s
      }
    } else {
      for (i in 1:S) {
        l_s = rnorm(TF,mean=MU_liab,sd=SD_liab)
        l_s[which(l_s<0)] = 0
        l_df[i,] = l_s
      }
    }
    return(l_df)
  } else {
    l_df1 = read.csv(f)
    l_df = l_df1[,-1]
    return(l_df)
  }
}


