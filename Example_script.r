source('Poisson_rate_estimation_functions_whole_biota_data_cf_Bacon_et_al_2015_PNAS.r')

# we simulate a vector with 20 timings of events arising uniformly from -5 to 0 time units
times=runif(n=20,min=-5,max=0)


############################################################
# 1) Constant Poisson rate
fun0=neg_lnL_poiss_constant(t=times)
ML1=optim(par=c(15),fn=fun0)
-ML1 # careful: ML1 is the opposite of the log-likelihood

############################################################
# 2) Exponential Poisson rate
fun_lin_poiss=neg_lnL_linear_poiss(t=times)
ML2=optim(par=c(0.1,1),fn=fun_lin_poiss)
-ML2

############################################################
# 3) Poisson rate shift at Pleistocene
fun_pleisto_poiss=neg_lnL_step_poiss(t=times)
ML3=optim(par=c(0.1,1),fn=fun_pleisto_poiss)
-ML3

############################################################
# 4) Poisson rate shift at -400,000 yrs

fun_400kyrs_poiss=neg_lnL_step_poiss_400kyrs(t=times)
ML4=optim(par=c(0.1,1),fn=fun_400kyrs_poiss)
-ML4


########################
# AIC comparison 4 models, corrected for small sample size
AICc1=2*(1+ML1$value)+2*1*(1+1)/(length(times)-1-1)
AICc2=2*(2+ML2$value)+2*2*(2+1)/(length(times)-2-1)
AICc3=2*(3+ML3$value)+2*2*(2+1)/(length(times)-2-1)
AICc4=2*(3+ML4$value)+2*2*(2+1)/(length(times)-2-1)

AICc1 ; AICc2 ; AICc3 ; AICc4 
ML1$par # the parameters of model 1 (probably the one with the lowest AICc) --> there is only one: the constant Poisson rate

