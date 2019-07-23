getImpliedRisk = function(
  stoch_mod, # result of stoch opt
  l_df, # simulated liability scenarios
  markDat # actually only risk-free rate is needed
) {
  # q_lev: active q = avar_p constraints
  active = which(stoch_mod$surplus==0)  # surplus == 0 means active >= constraint
  q_lev = as.numeric(stoch_mod$b2[active]) 
  shadow_prices = stoch_mod$a[(length(stoch_mod$soln)+length(stoch_mod$slack)+1):
                                length(stoch_mod$a)][active] # note the formatting of simplex.object
  # we need to reverse engineer the probability levels
  p_lev = numeric()
  l = sort(l_df %*% (1/(1+markDat$RF[1]/12)^(1:TF)))
  l_loren = sapply(l,function(s) {mean(l[which(l<=s)])})
  for (q in q_lev) {
    ind = which(l_loren>=q)[1]
    p_lev = c(p_lev,ind/length(l))
  } # take also the optimal langrage
  risk_aversion = sum(shadow_prices) # kappa: risk-aversion coefficient
  norm_shadow = shadow_prices/risk_aversion # normalized shawdow price
  res = list(
    lamda = norm_shadow,
    p_level = p_lev,
    kappa = risk_aversion
  )
  return(res)
}
