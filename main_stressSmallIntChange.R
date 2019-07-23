source("bond.R")
source("def.R")
source("liab.R")
source("cashFlow.R")
source("optimize.R")
source("stressSmallIntChange.R")
N_TRAIN=500
TF=36
N_TEST=10000
df0 = getDefDat(type="read")
markDat = getBondDat(type="read")
l_df = runSimLiab(type="gen",N_SIM=N_TRAIN+N_TEST)
l_train = l_df[1:N_TRAIN,]
simCF_train = runSimCF(markDat=markDat,defProb=df0,N_SIM = N_TRAIN,TF=TF)
stoch_mod = stochLiabMatch(l_df = l_train,simCF=simCF_train,markDat = markDat)
base_mod = regLiabMatch(l_df = l_train,simCF=simCF_train,markDat = markDat)


simCF_test = runSimCF(markDat=markDat,defProb=df0,N_SIM = N_TEST,TF=TF)
l_test = l_df[(N_TRAIN+1):(N_TRAIN+N_TEST),]

############# Interest rate models #######
{
require(forecast)
require(xts)
require(lubridate)
df=read.csv("DGS1MO.csv",colClasses = "character")
df=df[-which(df$DGS1MO=="."),]
r = xts(as.numeric(df$DGS1MO),order.by = ymd(df$DATE))
r = apply.monthly(r,FUN=function(x) {mean(na.omit(x))})
plot(r)

rising1 =auto.arima(r["2015/2019",],d=1,max.p=1,max.q = 0,ic="aic") # d=1 to fit vasicek
rising2 = auto.arima(r["2005/2007-03-30",],d=1,max.p=1,max.q = 0,ic="aic")
falling1=auto.arima(r["2000/2004-05",],d=1,max.p=1,max.q = 0,ic="aic")
falling2=auto.arima(r["2007-01/2008-12",],d=1,max.p=1,max.q = 0,ic="aic")
normal=auto.arima(r["2008-12/2014-12",],d=1,max.p=1,max.q = 0,ic="aic")
############################################
}

res_normal = runTest_smallInt(simCF=simCF_test,l_df=l_test,
                 stoch_mod = stoch_mod,base_mod=base_mod,
                 int_mod=normal,RF=markDat$RF[1])
res_rising = runTest_smallInt(simCF=simCF_test,l_df=l_test,
                              stoch_mod = stoch_mod,base_mod=base_mod,
                              int_mod=rising1,RF=markDat$RF[1])
res_falling = runTest_smallInt(simCF=simCF_test,l_df=l_test,
                                            stoch_mod = stoch_mod,base_mod=base_mod,
                                            int_mod=falling1,RF=markDat$RF[1])

