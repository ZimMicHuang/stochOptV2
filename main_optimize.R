source("bond.R")
source("def.R")
source("liab.R")
source("cashFlow.R")
source("optimize.R")

N_SIM=5000
TF=36
df0 = getDefDat(type="read")
markDat = getBondDat(type="read")
l_df = runSimLiab(type="read")
# defTimes = runSimDef(df0,save=FALSE,N_SIM=N_SIM)
# simDef = runSimDef(df0,save=FALSE,N_SIM=N_SIM)
simCF = runSimCF(markDat=markDat,defProb=df0,N_SIM = N_SIM,TF=TF)

stoch_mod = stochLiabMatch(l_df = l_df,simCF=simCF,markDat = markDat)
base_mod = regLiabMatch(l_df = l_df,simCF=simCF,markDat = markDat)
