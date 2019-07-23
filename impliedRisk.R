source("bond.R")
source("def.R")
source("liab.R")
source("cashFlow.R")
source("optimize.R")
N=5000
TF=36
df0 = getDefDat(type="read")
markDat = getBondDat(type="read")
l_df = runSimLiab(type="gen",N_SIM=N)
simCF = runSimCF(markDat=markDat,defProb=df0,N_SIM = N,TF=TF)
stoch_mod = stochLiabMatch(l_df = l_df,simCF=simCF,markDat = markDat)

q_lev = as.numeric(stoch_mod$b2)
# we need to reverse engineer the probability levels
p_lev = numeric()
l = sort(l_df %*% (1/(1+markDat$RF[1]/12)^(1:TF)))
l_loren = sapply(l,function(s) {mean(l[which(l<=s)])})
for (q in q_lev) {
  which(l_loren>q)[1]
  if (mean(l[which(l<=q)])%in%q_lev) {
    p_lev = c(p_lev,i/length(l))
  }
} # take also the optimal langrage
