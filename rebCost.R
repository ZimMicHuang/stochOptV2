getDoneComp = function(
  markDat, # bond data
  simCF, # simulated defaults
  s, # index of scenario
  time_change # time of interest rate change
) {
  defaulted_companies = c()
  paid_companies = c()
  N = length(simCF)
  tk = names(simCF)
  for (i in 1:N) {
    maturity=markDat$nper[i]
    # horizon closer by time_change
    if (maturity<=time_change) {
      # company_i bond matured already: cash flow remains 0
      paid_companies = c(paid_companies,tk[i])
    } else if (simCF[[i]][s,time_change+1]==0) {
      # company_i defaulted already: always default, cash flow remains 0
      # mark this by changing cp rate and recovery rate
      defaulted_companies = c(defaulted_companies,tk[i])
    } 
  }
  return(c(defaulted_companies,paid_companies))
}

updateMarkDat = function(
  markDat, # bond
  time_change, # timing of change
  mag_change # magnitude of change
) {
  require(optiRum)
  markDat_new = markDat
  markDat_new$RF = markDat$RF+mag_change
  markDat_new$ytm = markDat$ytm+mag_change
  for (i in 1:nrow(markDat)) {
    
    if (markDat$nper[i]>time_change) {
      # bond still active
      markDat_new$nper[i] = markDat$nper[i]-time_change
      markDat_new$mp[i] = optiRum::PV(rate=markDat_new$ytm[i]/12,
                                      nper=markDat_new$nper[i],
                                      pmt=-1*markDat_new$coupon[i],
                                      fv=-1)
    } else {
      # bond already paid off
      markDat_new$nper[i] = 0
      markDat_new$mp[i] = 999
    }
  }
  return(markDat_new)
}

updateLiab = function(
  l_df, # simulated liab
  time_change,
  mag_change
) {
  l_df_new = l_df[,(time_change+1):TF]
  return(l_df_new)
}

getPortVal = function(
  simCF, # simulated cash flows, stage 1
  markDat, # stage 1 bond data
  markDat_new, # stage 2 bond data
  stoch_mod, # stage 1 solution
  doneComp, # inactive companies
  s,
  time_change
) {
  # get end of stage 1 (before rebalancing) portfolio value
  N=length(simCF)
  cash = 0
  val = 0
  cf = numeric(time_change)
  for (i in 1:N) {
    cf = simCF[[i]][s,1:time_change]*stoch_mod$soln[i]
    cash=cash+sum( cf*(1+markDat$RF[1]/12)^(seq(from=time_change-1,to=0,by=-1)) )
    if (!(names(simCF)[i] %in% doneComp)) {
      val = val+stoch_mod$soln[i] * markDat_new$mp[i]
    }
  }
  return(as.numeric(cash+val))
}
getRebCost = function(
  time_change, # timing of interest rate change
  mag_change, # magnitude of interest rate change
  N_SIM, # number of scenarios to simulate
  TF, # investment horizon
  df0, # default probability estimation
  markDat, # stage 1 market data
  l_df, # simulated liabilities
  parallel = FALSE, # use parallel computing
  cores = 8 # number of cores to use
) {
  time_change_ = time_change
  mag_change_ = mag_change
  N_SIM_ = N_SIM
  TF_ = TF
  markDat1 = markDat
  l_df1 = l_df
  # stage 1 #############
  simCF1 = runSimCF(markDat=markDat1,defProb=df0,N_SIM = N_SIM_,TF=TF_)
  N_COMP = length(df0)
  stoch_mod1 = stochLiabMatch(l_df = l_df1,simCF=simCF1,markDat = markDat1)
  base_mod1 = regLiabMatch(l_df = l_df1,simCF=simCF1,markDat = markDat1)
  # Stage 2 #############
  markDat2 = updateMarkDat(markDat1,time_change=time_change_,
                              mag_change = mag_change_)
  l_df2 = updateLiab(l_df=l_df1,time_change=time_change_,
                        mag_change = mag_change_)
  
  new_port_stoch = matrix(nrow=N_SIM_,ncol=N_COMP)
  colnames(new_port_stoch)=names(df0)
  reb_cost_stoch=matrix(nrow=N_SIM_,ncol=2)
  colnames(reb_cost_stoch)=c("preReb","postReb")
  
  new_port_reg = matrix(nrow=N_SIM_,ncol=N_COMP)
  colnames(new_port_reg)=names(df0)
  reb_cost_reg=matrix(nrow=N_SIM_,ncol=2)
  colnames(reb_cost_reg)=c("preReb","postReb")
  
  
  
  if (parallel) {
    require(doParallel)
    cl = makeCluster(cores)
    registerDoParallel(cl)
    
    tmp = foreach (scenario = 1:N_SIM_,
                   .export = c("getDoneComp","getPortVal","runSimCF",
                               "stochLiabMatch","regLiabMatch","runSimDef"),
                   .packages = c("boot")) %dopar% {
      doneComp1 = getDoneComp(markDat=markDat1,simCF = simCF1,s=scenario,time_change = time_change_)
      
      markDatTmp = markDat2
      if (!is.null(doneComp1)) {
        for (tk in doneComp1) {
          index = which(rownames(markDatTmp)==tk)
          markDatTmp$coupon[index]=0
          markDatTmp$rec[index]=0
        }
      }
      
      befReb_stoch = getPortVal(simCF=simCF1,markDat=markDat1,markDat_new = markDatTmp,
                                stoch_mod=stoch_mod1,doneComp=doneComp1,s=scenario,time_change=time_change_)
      befReb_reg = getPortVal(simCF=simCF1,markDat=markDat1,markDat_new = markDatTmp,
                              stoch_mod=base_mod1,doneComp=doneComp1,s=scenario,time_change=time_change_)
      simCF2 = runSimCF(markDat=markDatTmp,
                        defProb=df0,
                        N_SIM = N_SIM_,
                        TF=TF_-time_change_,
                        doneComp = doneComp1)
      
      stoch_mod2 = stochLiabMatch(l_df=l_df2,simCF=simCF2,markDat=markDatTmp,EPS=0.005,plotit=FALSE)
      base_mod2 = regLiabMatch(l_df=l_df2,simCF=simCF2,markDat=markDatTmp)
      list(
        new_port_stoch = stoch_mod2$soln,
        reb_cost_stoch = c(befReb_stoch,stoch_mod2$value),
        
        new_port_reg = base_mod2$soln,
        reb_cost_reg = c(befReb_reg,base_mod2$value)
      )
    }
    print("parallell computing successful!")
    for (s in 1:N_SIM_) {
      new_port_stoch[s,]=tmp[[s]]$new_port_stoch
      reb_cost_stoch[s,]=tmp[[s]]$reb_cost_stoch
      
      new_port_reg[s,]=tmp[[s]]$new_port_reg
      reb_cost_reg[s,]=tmp[[s]]$reb_cost_reg
    }
    return(list(stoch=list(new_port=new_port_stoch,reb_cost=reb_cost_stoch),
                reg=list(new_port=new_port_reg,reb_cost=reb_cost_reg)))
  } else {
    for (scenario in 1:N_SIM_) {
      doneComp1 = getDoneComp(markDat=markDat1,simCF = simCF1,s=scenario,time_change = time_change_)
      
      markDatTmp = markDat2
      if (!is.null(doneComp1)) {
        for (tk in doneComp1) {
          index = which(rownames(markDatTmp)==tk)
          markDatTmp$coupon[index]=0
          markDatTmp$rec[index]=0
        }
      }
      
      befReb_stoch = getPortVal(simCF=simCF1,markDat=markDat1,markDat_new = markDatTmp,
                          stoch_mod=stoch_mod1,doneComp=doneComp1,s=scenario,time_change=time_change_)
      befReb_reg = getPortVal(simCF=simCF1,markDat=markDat1,markDat_new = markDatTmp,
                              stoch_mod=base_mod1,doneComp=doneComp1,s=scenario,time_change=time_change_)
      simCF2 = runSimCF(markDat=markDatTmp,
                        defProb=df0,
                        N_SIM = N_SIM_,
                        TF=TF_-time_change_,
                        doneComp = doneComp1)
      
      stoch_mod2 = stochLiabMatch(l_df=l_df2,simCF=simCF2,markDat=markDatTmp,EPS=0.005,plotit=FALSE)
      base_mod2 = regLiabMatch(l_df=l_df2,simCF=simCF2,markDat=markDatTmp)
      new_port_stoch[scenario,]=stoch_mod2$soln
      reb_cost_stoch[scenario,1]=befReb_stoch
      reb_cost_stoch[scenario,2]=stoch_mod2$value
      new_port_reg[scenario,]=base_mod2$soln
      reb_cost_reg[scenario,1]=befReb_reg
      reb_cost_reg[scenario,2]=base_mod2$value
      print(scenario)
    }
    return(list(stoch=list(new_port=new_port_stoch,reb_cost=reb_cost_stoch),
           reg=list(new_port=new_port_reg,reb_cost=reb_cost_reg)))
  }
}
