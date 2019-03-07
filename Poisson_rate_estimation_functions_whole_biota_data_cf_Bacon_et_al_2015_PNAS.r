############################################################
# The R functions below implement different forms of the likelihood functions used by Bacon et al. 2015 Proc. Nat. Acad. Sci.
# They model a series of dated events (e.g. divergence times, dated migration events) as arising from a Poisson process and mesure the most likely (Poisson) rate at which these events happened
# This model is typically useful for studying the rate at which some events (e.g. divergence, dispersal) happened across an entire biota (i.e. in widely different clades) rather than within a single clade for which we would have a phylogeny
# The derivation of this likelihood function, as well as some important discussion of its limitations, is described in the Supplementary Information of the paper by Bacon et al. 

# In all of these R functions, we have the following parameters:
# -'times' is a vector of times in the past (they must all be negative!) that need not be ordered 
# - 'lambda0' is a numerical value that gives the constant Poisson rate across time
# - 'lambda' is a function that gives the shape of the Poisson rate through time

############################################################
############################################################
############# GENERAL LIKELIHOOD FUNCTIONS #################
############################################################
############################################################

############################################################
# 1) The most general function: can accomodate any shape of the Poisson rate through time (lambda)
loglikelihood_any_poiss=function(times,lambda){
  if (max(times)>0) {stop('Divergence times must all be negative!')}
  else{
    timessorted=c(sort(times),0)
    lnL=0
    for (i in 1:length(times)){
      temp=-integrate(f=lambda,lower=timessorted[i],upper=timessorted[i+1])$value+log(lambda(timessorted[i]))
      lnL=lnL+temp
    }
    return(lnL)
  }
}

############################################################
# 2) A very special case of the general function: constant Poisson rate through time (lambda(t)=lambda0), will often be used as a null model
loglikelihood_poiss_constant=function(times,lambda0){
  if (max(times)>0) {stop('Divergence times must all be negative!')}
  else{
    timessorted=c(sort(times),0)
    lnL=0
    for (i in 1:length(times)){
      temp=log(exp(-lambda0*(timessorted[i+1]-timessorted[i])))+log(lambda0)
      lnL=lnL+temp
    }
    return(lnL)
  }
}

###########################################################################
###########################################################################
############# FUNCTIONS FOR MAXIMUM LIKELIHOOD ESTIMATION #################
###########################################################################
###########################################################################

# All four functions below build on the two likelihood functions above: given the observed vector of times 't' they build a case-specific likelihood function that only depends on the Poisson rate (i.e. times are fixed since they were observed) 
# Be carefull: they all return the opposite of the log-likelihood!!! This is done because you will probably optimize the parameters if these functions using the R base function 'optim', which actually does minimization --> by searching the minimum of the opposite of the likelihood you're searching the maximumn of the likelihood function!

############################################################
# 1) Constant Poisson rate through time
neg_lnL_poiss_constant=function(t){
  FUN=function(x){return(-loglikelihood_poiss_constant(times=t,lambda0 = x))}
  return(FUN)
}
  
############################################################
# 2) Exponential Poisson rate through time
neg_lnL_expo_poiss=function(t){
  FUN=function(X){return(-loglikelihood_any_poiss(times=t,lambda = function(t){return(X[1]+X[2]*exp(-X[2]*t))}))}
  return(FUN)
}

############################################################
# 3) Poisson rate shift at Pleistocene (-2.1 Ma): there is a constant Poisson rate before -2.1 and another constant Poisson rate after -2.1 (until the present)
# Be careful with the units in which your times are measured. 
# You can of course change the value of '-2.1' to any value of interest to test for a shift in the Poisson rate, as done below in function 4)
neg_lnL_step_poiss=function(t){
  FUN=function(X){return(-loglikelihood_any_poiss(times=t,lambda=stepfun(x=-2.1,y=X)))}
  return(FUN)
}

############################################################
# 4) Poisson rate shift at -400,000 yrs
# Same spirit as in 3)
neg_lnL_step_poiss_400kyrs=function(t){
  FUN=function(X){return(-loglikelihood_any_poiss(times=t,lambda=stepfun(x=-0.4,y=X)))}
  return(FUN)
}
