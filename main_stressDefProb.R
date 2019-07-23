source("bond.R")
source("def.R")
source("liab.R")
source("cashFlow.R")
source("optimize.R")
source("stressDefProb.R")
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

# systematically overestimate default probability
df1 = (df0+rnorm(n=length(df0),mean=0,sd=0.3*min(df0)))*0.4
simCF_test = runSimCF(markDat=markDat,defProb=df1,N_SIM = N_TEST,TF=TF)
l_test = l_df[sample.int(n=N_TRAIN+N_TEST,size=N_TEST),]

res_smallDef = runTest_defProb(simCF=simCF_test,l_df=l_test,
                 stoch_mod = stoch_mod,base_mod=base_mod,RF=markDat$RF[1])


# good estimate on average, noise added
df2 = df0+rnorm(n=length(df0),mean=0,sd=0.3*min(df0))
simCF_test = runSimCF(markDat=markDat,defProb=df2,N_SIM = N_TEST,TF=TF)
l_test = l_df[sample.int(n=N_TRAIN+N_TEST,size=N_TEST),]

res_normDef = runTest_defProb(simCF=simCF_test,l_df=l_test,
                               stoch_mod = stoch_mod,base_mod=base_mod,RF=markDat$RF[1])



# under estimate 
df3 = (df0+rnorm(n=length(df0),mean=0,sd=0.3*min(df0)))*1.5
simCF_test = runSimCF(markDat=markDat,defProb=df3,N_SIM = N_TEST,TF=TF)
l_test = l_df[sample.int(n=N_TRAIN+N_TEST,size=N_TEST),]

res_largeDef = runTest_defProb(simCF=simCF_test,l_df=l_test,
                              stoch_mod = stoch_mod,base_mod=base_mod,RF=markDat$RF[1])
