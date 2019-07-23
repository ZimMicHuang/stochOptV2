source("bond.R")
source("def.R")
source("liab.R")
source("cashFlow.R")
source("optimize.R")
source("impliedRisk.R")
N=5000
TF=36
df0 = getDefDat(type="read")
markDat = getBondDat(type="read")
n_iter = 100
p_lev = matrix(nrow=n_iter,ncol=100)
kappa = numeric(n_iter)
lambda = matrix(nrow=n_iter,ncol=100)
for (i in 1:n_iter) {
  l_df = runSimLiab(type="gen",N_SIM=N)
  simCF = runSimCF(markDat=markDat,defProb=df0,N_SIM = N,TF=TF)
  stoch_mod = stochLiabMatch(l_df = l_df,simCF=simCF,markDat = markDat)
  res = getImpliedRisk(stoch_mod = stoch_mod,l_df=l_df,markDat = markDat)
  lambda[i,] = c(res$lamda,rep(0,times=100-length(res$lamda)))
  p_lev[i,] = c(res$p_lev,rep(0,times=100-length(res$p_lev)))
  kappa[i] = res$kappa
}