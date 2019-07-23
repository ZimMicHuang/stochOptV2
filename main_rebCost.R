source("bond.R")
source("def.R")
source("liab.R")
source("cashFlow.R")
source("optimize.R")
source("rebCost.R")
N_SIM=400
TF=36
df0 = getDefDat(type="read")
markDat = getBondDat(type="read")
l_df_all = runSimLiab(type="gen",N_SIM = 10000)



time_change_vec = seq(1,22,by=2)
mag_change_vec = seq(-0.05,0.05,by=0.01)
reb_cost_grid = expand.grid(time_change=time_change_vec,
                            mag_change=mag_change_vec)
reb_cost_grid = cbind(reb_cost_grid,reb_cost_stoch=numeric(length(time_change_vec)*length(mag_change_vec)),reb_cost_reg = numeric(length(time_change_vec)*length(mag_change_vec)))
rc=0
for (i in 1:nrow(reb_cost_grid)) {
  l_df = l_df_all[sample(x=1:10000,size=N_SIM,replace=FALSE),]
  tc = reb_cost_grid[i,1]
  mc = reb_cost_grid[i,2]
  rc = 0
  res = getRebCost(time_change=tc,mag_change=mc,N_SIM = N_SIM,TF=TF,
                   df0=df0,markDat = markDat,l_df=l_df,
                   parallel=TRUE,cores=6)
  av_reb_stoch = mean(res$stoch$reb_cost[,1]-res$stoch$reb_cost[,2]) # negative means cash outflow
  av_reb_reg = mean(res$reg$reb_cost[,1]-res$reg$reb_cost[,2])
  reb_cost_grid[i,3] = av_reb_stoch
  reb_cost_grid[i,4] = av_reb_reg
  print(reb_cost_grid[i,])
  print(i)
}