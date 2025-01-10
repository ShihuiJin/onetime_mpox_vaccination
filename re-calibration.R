#model calibration using the DRC data

#data preparation----------
#vaccination coverage for ppl born in 1901-2000
vacc_c=coverage_birth_mean[['COD']]
n.sim=1000

#DRC surveillance data & contact matrices
#from the manuscript by Murayama et al.
load('data/DRC_data.RData')

#functions needed
source('functions.R')
library(parallel)
library(stringr)
library(dplyr)
library(purrr)
library(hhh4contacts)
library(parallel)
numCores <- detectCores()

cms=list(drc2015_home, drc2024_home) #synthetic, home contact
cases=list(tshuapa_h2hag$cases, drc_endemic_ag$cases)
times=c(2015,2024)

#calibration-----------
#step 1: susceptibility
paras2=mclapply(1:n.sim, function(i){
  resampling=function(cases,seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    rmultinom(1,sum(cases),cases/sum(cases))
  }
  fit_ngm_pair=function(x){
    do.call('sum',lapply(1:length(cms), function(k){
      cm=cms[[k]]
      cm$parameters$s_infant=x[1]
      cm$parameters$s_vax=x[2]
      suscept=immu_fun(cm, times[k],'COD')
      fit=eigen(ngm_cal(suscept,cm$matrix))
      vec=abs(fit$vectors[,1])
      vec=vec/sum(vec)
      -dmultinom(resampling(cases[[k]], i), prob=vec, log=T)
    }))
  }
  opt=optim(c(1,1), fit_ngm_pair, method='L-BFGS-B', 
            lower=rep(1e-6,2), hessian = T)
  opt$par
}, mc.cores = numCores)%>%do.call('rbind',.)%>%as.matrix

#step 2: sexual-contact-related
cm=cms[[2]]
cm=addsexualcontact(cm, c(15,seq(20,40,10)),countryname = 'COD', year = 2024)
cm=propmix(cm)
cases_sex=unlist(lapply(1:2, function(i) drc_kivu$cases[i,]))

sex_paras2=mclapply(1:n.sim, function(i){
  resampling=function(cases,seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    rmultinom(1,sum(cases),cases/sum(cases))
  }
  fit_ngm_sex=function(y){
    par_names=names(cm$parameters)
    par_names=par_names[str_detect(par_names,'_addmat')]
    cm$parameters[par_names]=y[1:length(par_names)+2]
    cm$parameters$addmat_v=y[1]
    cm$parameters$addmat_w=y[2]
    cm=propmix(cm)
    suscept=age_suscept(paras2[i,], cm, 2024,'COD')
    suscept=rep(suscept,4)
    fit=eigen(ngm_cal(suscept,cm$sexmat))
    vec=abs(fit$vectors[,1])
    vec=vec/sum(vec)
    n.age=length(cm$ageinterval)
    vec1=vec[c(1:n.age, 1:n.age+n.age*2)]+vec[n.age+c(1:n.age, 1:n.age+n.age*2)]
    -dmultinom(resampling(cases_sex, i), prob=vec1, log=T)-propmix_ll(cm)
  }
  y=c(30,30,rep(1.1,8))
  opt_sex=optim(y, fit_ngm_sex, 
                method='L-BFGS-B', lower=rep(1e-3,length(y)),
                control = list(maxit=400),hessian = T)
  opt_sex$par
}, mc.cores = numCores)%>%do.call('rbind',.)%>%as.matrix

#output----------
paras=cbind(paras2,sex_paras2)%>%as.matrix() #for 2-step approach
#find the R0_scales for each parameter set
R0_scales=lapply(1:n.sim, function(i){
  cm0=cm_setup(2015, 'COD', paras[i,],sex=F)
  ngm_vacc(cm0, sex=F)$val/0.82
})%>%unlist()
save(paras, R0_scales, file="data/fitted_DRC_1k_boot.RData")

