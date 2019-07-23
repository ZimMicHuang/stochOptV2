getBondDat = function(
  type=c("read","gen"), 
  f="markDat.csv",
  defProb, # simulated default probabilities for each company
  RF = 0.03, # risk free rate
  CREDIT_RISK_PREMIUM=0.025, # parameters for generating bond prices (see paper by John Hull)
  TF_UB = 36, # longest maturity
  TF_LB = 6, #shortest maturity
  cp_LB = 0, #lowest coupon
  cp_UB = 0.15, #highest coupon
  rec_LB=0.2, # lowest recovery rate
  rec_UB=0.6 # largest recovery rate
) {
  require(optiRum)
  if (type=="gen") {
    N=length(defProb)
    # randomize: coupon, recovery, periods
    cp=runif(N,cp_LB,cp_UB)/12 # coupon rates
    rec = runif(N,rec_LB,rec_UB) # recovery rates
    nper = sample(x=TF_LB:TF_UB,size=N,replace=TRUE) # maturity (months)
    
    # fix: recovery, default prob, premium --> ytm
    ytm = numeric(N)
    for (i in 1:N) {
      p = defProb[i]
      ytm[i] = p*(1-rec[i])+RF #risk neutral yield
    }
    ytm = ytm+rnorm(N,mean=CREDIT_RISK_PREMIUM,sd=0.2*CREDIT_RISK_PREMIUM)
    # solve: ytm, coupon, periods, fv (normalized to 1) --> price
    # in reality ytm is implied, and market price is known directly
    markPrice = optiRum::PV(
      rate=ytm/12,
      nper=nper,
      pmt=-1*cp,
      fv=-1
    )
    markDat = data.frame(coupon=cp, 
                         mp = markPrice, 
                         rec = rec, 
                         ytm=ytm,
                         nper=nper,
                         RF=RF
    )
    rownames(markDat) = names(defProb)
    return(markDat)
  } else {
    markDat = read.csv(f,row.names = 1)
    # ytm MUST BE CORRECT
    return(markDat)
  }
  
}