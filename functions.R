#functions
#data loading---------
future_pop=function(country='COD',year=2025, age_end=75,gender=F){
  year=min(year, 2050)
  pop_dat=proj_pop%>%filter(ISO==country)%>%
    dplyr::select(all_of(c('ISO','Age','Sex',paste0('X',year))))%>%
    rename(n=paste0('X',year))
  age_end=age_end/5+1
  {
    pop_raw=as.matrix(do.call('cbind',lapply(c('MA','FE'), function(i){
      pop_dat%>%filter(Sex==i)%>%pull(n)%>%as.numeric()
    })))
    pop_raw[age_end,]=colSums(pop_raw[-c(1:age_end-1),])
    pop_raw=pop_raw[1:age_end,]/1e3
  }
  if(!gender) pop_raw=rowSums(pop_raw)
  pop_raw
}

merge_future_pop=function(ageinterval, country='COD', year=2023, gender=F){
  pop_raw=future_pop(country, year, 75,gender)
  ind=match(ageinterval, 0:15*5)
  ind=c(ind, 17)
  if(!gender){
    unlist(lapply(2:length(ind), function(i){
      sum(pop_raw[ind[i-1]:(ind[i]-1)]) }))
  }else{
    as.matrix(do.call('cbind',lapply(1:2, function(j){
      unlist(lapply(2:length(ind), function(i){
        sum(pop_raw[ind[i-1]:(ind[i]-1),j]) }))
    })))
  }
}

#summary function-------------
estim=function(dat){
  if(is.vector(dat)){
    c(mean(dat),quantile(dat, c(0.025,0.975,0.5,0.25,0.75)))
  }else{
    #matrix
    lapply(1:ncol(dat), function(i){
      c(mean(dat[,i]),quantile(dat[,i], c(0.025,0.975,0.5,0.25,0.75)))
    })%>%do.call('rbind',.)%>%as.matrix
  }
}

#model fitting functions-------------
#calculate vaccine coverage rate for specific age groups
coverage=function(ageint,year, country='COD',vacc_p=vacc_c){
  #first convert to age groups of 5 yrs
  vacc_p=c(vacc_p,rep(0,100))
  full_cov=unlist(lapply(0:15, function(t){
    if(t<15){
      mean(vacc_p[year-t*5-1900+(-4):0])
    }else{
      mean(vacc_p[1:(year-t*5-1900)])
    }
  }))
  pop_full=future_pop(country, year,age_end = 75)
  pop_full=pop_full/sum(pop_full)
  #merge age groups
  ind=match(ageint, 0:15*5)
  ind=c(ind, 17)
  unlist(lapply(2:length(ind)-1, function(i){
    idx=ind[i]:(ind[i+1]-1)
    sum(pop_full[idx]*full_cov[idx])/sum(pop_full[idx]) #weighted average
  }))
}

#next generation matrix
ngm_cal=function(suscept, matrix){
  diag(suscept)%*%t(matrix)
}

#age(interval)-specific immunity based on vaccine coverage rates (by birth year)
immu_fun=function(cm, year=2024, country='COD',vacc_p=vacc_c){
  vacc=coverage(cm$ageinterval,year,country, vacc_p)
  out0=1-cm$parameters$s_vax*vacc
  out0[1]=out0[1]*cm$parameters$s_infant
  out0
}

#incorporate sexual contact parameters
construct_sex=function(na,pa, coef, Q_vec=NULL){
  vec=coef*na*pa/sum(na*pa)
  if(is.null(Q_vec)) Q_vec=rep(1, length(na))
  mt=as.matrix(do.call('cbind',lapply(1:length(pa), function(i){
    if(pa[i]==0){
      rep(0,length(pa))
    }else{
      vec*Q_vec[i]
    }
  })))
  mt
}

addsexualcontact <- function(cmt, activeage = seq(10, 40, by = 10), countryname = "COD", year = 2024) {
  ind <- which(cmt$ageinterval %in% activeage)
  parkeys <- paste0(rep(c("1_", "2_"), each = length(ind)), "addmat", ind)
  cmt$parameters[parkeys] <- 1.1
  cmt$parameters$addmat_v <- 30
  cmt$parameters$addmat_w <- 30
  cmt$misc$pop <- merge_future_pop(cmt$ageinterval, country = countryname, year = year, gender = T)
  cmt$misc$p_sex <- sex_data%>%filter(ISO==countryname)%>%pull(p_sex)
  return(cmt)
}

#high risk proportion
cal_sex=function(cm) {
  addmatkeys <- grep("_addmat", names(cm$parameters), value = TRUE)
  ind <- strsplit(addmatkeys, "_addmat")
  ind <- lapply(ind, function(x) as.integer(x[2]))
  paramat <- matrix(sapply(addmatkeys, function(key) 
    cm$parameters[[key]]), nrow = length(unique(sapply(ind, `[`, 1))))
  paramat <- (1 - paramat)^2 / 2 #some transformation-->proportion of active population in each age group
  if(is.null(cm$parameters$scaling)){
    scaling=rep(1,2)
  }else{
    scaling=cm$parameters$scaling
  }
  #proportion list
  paqa <- lapply(1:2, function(i) rep(0, length(cm$ageinterval))) 
  paqa_substitute <- split(paramat, col(paramat))
  for (i in 1:2) {
    paqa[[i]][unique(sapply(ind, `[`, 1))] <- paqa_substitute[[i]]*scaling[i]
  }
  paqa
}

#derive sexual contact matrix
#add scaling for different countries
propmix <- function(cm) {
  paqa=cal_sex(cm)
  
  na <- cm$misc$pop #by gender
  sum_napaqa <- sapply(seq_along(paqa), function(x) sum(na[,x] * paqa[[x]]))
  
  vw <- c(cm$parameters$addmat_v, cm$parameters$addmat_w) #w_F, w_M
  v0 <- vw[1] / ((sqrt(12^2 + 22^2) / (17.61 + 5.51))^2 + 1) #coeffient for Sigma_MF
  w0 <- v0 * (sum_napaqa[2] / sum_napaqa[1]) #coeffient for Sigma_FM
  
  cm$misc$high_p=sum_napaqa/colSums(na[4:7,])
  
  s_cmt=list(
    construct_sex(na[,1],paqa[[1]],vw[1]), #S_MF
    construct_sex(na[,2],paqa[[2]],vw[2]), #S_FM
    construct_sex(na[,1],paqa[[1]],v0, paqa[[2]]), #Sigma_MF*Q
    construct_sex(na[,2],paqa[[2]],w0, paqa[[1]]) #Sigma_FM*P
  )
  n.col=length(paqa[[1]])
  null_mat=matrix(0,n.col,n.col)
  p_male=cm$misc$pop[,1]/rowSums(cm$misc$pop)
  contact_mats=as.matrix(do.call('cbind',lapply(1:4, function(i){
    if(i<=2){
      t(cm$matrix)%*%diag(p_male)
    }else{
      t(cm$matrix)%*%diag(1-p_male)
    }
  })))
  cm$sexmat <- rbind(
    cbind(null_mat,null_mat, s_cmt[[1]],s_cmt[[3]]),
    contact_mats,
    cbind(s_cmt[[2]],s_cmt[[4]],null_mat,null_mat),
    contact_mats
  )%>%as.matrix%>%t()
  
  cm
}

#likelihood for boundary conditions
propmix_ll=function(cm, bcond = c(0.1, 0.02)){
  if(!is.null(cm$misc$p_sex)) bcond[2]=cm$misc$p_sex
  -sum(((log(cm$misc$high_p) - log(bcond)) / log(bcond))^2) * 10000
}

#age group-specific susceptbility
age_suscept=function(x,cm, year,country='COD',vacc_p=vacc_c){
  cm$parameters$s_infant=x[1]
  cm$parameters$s_vax=x[2]
  immu_fun(cm,year,country,vacc_p)
}

#fit the optimized values
ngm_fit=function(xy, cm, time, sex=T,country='COD'){
  cm$parameters$s_infant=xy[1]
  cm$parameters$s_vax=xy[2]
  suscept=immu_fun(cm, time,country,vacc_c)
  if(sex==F){ #community contact model
    ngm=ngm_cal(suscept, cm$matrix)
    fit=eigen(ngm)
    vec=abs(fit$vectors[,1])
    vec=vec/sum(vec)
    val=max(abs(fit$values))
    list(val=val,vec=vec,ngm=ngm)
  }else{
    cm=addsexualcontact(cm, c(15,seq(20,40,10)),countryname = country, year = time)
    par_names=names(cm$parameters)
    par_names=par_names[str_detect(par_names,'_addmat')]
    cm$parameters[par_names]=xy[1:length(par_names)+4]
    cm$parameters$addmat_v=xy[3]
    cm$parameters$addmat_w=xy[4]
    cm=propmix(cm)
    suscept=rep(suscept,4)
    ngm=ngm_cal(suscept,cm$sexmat)
    fit=eigen(ngm)
    vec0=abs(fit$vectors[,1]) #unscaled
    vec=vec0/sum(vec0)
    n.age=length(cm$ageinterval)
    vec1=vec[c(1:n.age, 1:n.age+n.age*2)]+vec[n.age+c(1:n.age, 1:n.age+n.age*2)]
    p_sex=vec[c(1:n.age, 1:n.age+n.age*2)]%>%sum
    val=max(abs(fit$values))
    list(val=val, vec_combined=vec1, p_sex=p_sex,vec=vec0, ngm=ngm)
  }
}

#projection, for all Sub-Saharan countries----------
#sex: whether sexual transmission is modeled
cm_setup=function(year=2025, country='COD', xy=full_model_fits, 
                  sex=T, activeage=c(15,seq(20,40,10))){
  ageint=drc_endemic_ag$ageinterval
  vacc_c=coverage_birth_mean[[country]]
  cmt=contact_home[[country]]
  rownames(cmt)=0:15*5
  pop=future_pop(country, year, age_end=75) #combined population, also allow for future predictions
  temp=match(ageint,0:15*5)
  groupmap=c(temp[-1], 17)-temp
  cmt_new <- aggregateC(cmt, groupmap, pop)
  new_order=rownames(cmt_new)%>%as.numeric%>%order
  cmt_new=cmt_new[new_order, new_order]
  cm=ContactMatrix(cmt_new, ageint, rep(1,length(ageint)),list(), 
                   list(issynthetic = T))
  cm$parameters$s_infant=xy[1]
  cm$parameters$s_vax=xy[2]
  cm$susceptibility=immu_fun(cm, year,country,vacc_c) #baseline susceptiblity, without new vaccination
  pop <- merge_future_pop(cm$ageinterval,country,year, T)
  cm$misc$pop=pop
  if(sex>0){
    p_sex=sex_data%>%filter(ISO==country)%>%pull(p_sex)
    cm$misc$p_sex=p_sex
    idx <- which(cm$ageinterval %in% activeage)
    parkeys <- paste0(rep(c("1_", "2_"), each = length(idx)), "addmat", idx)
    cm$parameters[parkeys] <- xy[1:length(parkeys)+4]
    cm$parameters$addmat_v <- xy[3]
    cm$parameters$addmat_w <- xy[4]*sex #scaling
    paqa=cal_sex(cm)
    DRC_pop=merge_future_pop(cm$ageinterval,'COD',2024,gender=T)
    #scaling
    k2=p_sex/sum(pop[,2]*paqa[[2]])*sum(pop[idx,2]) #female
    k1=sum(pop[,2]*paqa[[2]])*k2/sum(DRC_pop[,2]*paqa[[2]])*sum(DRC_pop[,1]*paqa[[1]])/sum(pop[,1]*paqa[[1]])
    cm$parameters$scaling <- c(k1,k2)
    cm$misc$paqa=cal_sex(cm)
  }
  cm
}

#calculating actual vaccination rate based on target coverage rates each year
#vaccinating high-sexual-activity female (vacc_sex), general population (all), in percentage from 0 to 100
vacc_rate=function(cm=cm, age=c(20,30),v_sex=NULL, v_all=NULL,vacc_plan=NULL,exclude=T){
  n.age=length(cm$ageinterval)
  #initialization
  if(is.null(vacc_plan)){
    vacc_plan=list( #vaccine uptake rate by age group
      ageinterval=cm$ageinterval,
      vacc_sex=rep(0,n.age), #high-sexual-acitivty females
      vacc_allm=rep(0,n.age),
      vacc_allf=rep(0,n.age) 
    ) 
  }
  if(is.null(cm$misc$paqa)) cm$misc$paqa=cal_sex(cm)
  idx=match(age,cm$ageinterval) #sexually active index
  if(is.null(v_sex)&is.null(v_all)) return(vacc_plan) #no vaccination
  p_sexf=cm$misc$paqa[[2]][idx]
  n.age=length(idx)
  if(is.null(v_all)) v_all=rep(0, n.age)
  if(is.null(v_sex)) v_sex=rep(0, n.age)
  v_sex=unlist(lapply(1:n.age, function(k){
    if(p_sexf[k]==0) return(0)
    max(v_sex[k],v_all[k],vacc_plan$vacc_sex[idx][k])
  }))
  v_allm=unlist(lapply(1:n.age, function(k){
    max(v_all[k],vacc_plan$vacc_allm[idx][k])
  }))
  former_allf=(vacc_plan$vacc_allf[idx]-vacc_plan$vacc_sex[idx]*p_sexf)/(1-p_sexf)
  if(exclude){
    #target overall uptake rate exclude female sex workers
    v_allf=unlist(lapply(1:n.age, function(k){
      v_sex[k]*p_sexf[k]+max(v_all[k],former_allf[k])*(1-p_sexf[k])
    }))
  }else{
    #target overall uptake rate exclude female sex workers
    v_allf=unlist(lapply(1:n.age, function(k){
      max(v_sex[k]*p_sexf[k],v_all[k],vacc_plan$vacc_allf[idx][k])
    }))
  }
  vacc_plan$vacc_sex[idx]=v_sex
  vacc_plan$vacc_allf[idx]=v_allf
  vacc_plan$vacc_allm[idx]=v_allm
  return(vacc_plan)
}

# % protected due to vaccination, under different assumed VE
vaccine_protect=function(suscept0,uptake, ve=0.8){
  #uptake in percetage (i.e., 0-100)
  p.prot=uptake/100*ve
  out=c(suscept0[1]*(1-p.prot[1]), (suscept0*(1-p.prot))[-1])
  out[out<0]=0
  out
}

#next generation matrix for vaccination
ngm_vacc=function(cm, vacc_plan=NULL, sex=T, ve=0.85, R0_scale=1){
  suscept=cm$susceptibility
  if(is.null(vacc_plan)){
    n.age=length(cm$ageinterval)
    vacc_plan=list( #vaccine uptake rate by age group
      ageinterval=cm$ageinterval,
      vacc_sex=rep(0,n.age), 
      vacc_allm=rep(0,n.age),
      vacc_allf=rep(0,n.age) 
    ) 
  }
  if(sex==F){ #community contact model
    vaccined=(vacc_plan$vacc_allf+vacc_plan$vacc_allm)/2
    suscept1=vaccine_protect(suscept,vaccined,ve)
    ngm=ngm_cal(suscept1,cm$matrix)
    fit=eigen(ngm)
    vec=abs(fit$vectors[,1])
    vec=vec/sum(vec)
    val=max(abs(fit$values))
    list(val=val/R0_scale,vec=vec,ngm=ngm)
  }else{ #with sexual transmission
    cm=propmix(cm)
    suscept1=c(vaccine_protect(suscept,vacc_plan$vacc_allm,ve),
               vaccine_protect(suscept,vacc_plan$vacc_allm,ve),
               vaccine_protect(suscept,vacc_plan$vacc_sex,ve),
               vaccine_protect(suscept,vacc_plan$vacc_allf,ve))
    ngm=ngm_cal(suscept1,cm$sexmat)
    fit=eigen(ngm)
    vec0=abs(fit$vectors[,1]) #unscaled
    vec=vec0/sum(vec0)
    n.age=length(cm$ageinterval)
    vec1=vec[c(1:n.age, 1:n.age+n.age*2)]+vec[n.age+c(1:n.age, 1:n.age+n.age*2)]
    p_sex=vec[c(1:n.age, 1:n.age+n.age*2)]%>%sum
    val=max(abs(fit$values))
    list(val=val/R0_scale, vec_combined=vec1, p_sex=p_sex,vec=vec0, ngm=ngm)
  }
}

#vaccine allocation functions-----------
#determine the best plan through gradient descent-------------
#thre: max vaccine uptake rate; m: step size (in thousand)
find_optim_plan=function(cm,thre=80,m=10, sex=T, n.vac=8,scale=R0_scale){
  init=rep(0,n.vac)
  R0_doses=function(v_rates){
    vacc_plan=vacc_rate(cm, age=cm$ageinterval[1:n.vac], v_sex=rep(0,n.vac), v_all=v_rates)
    res=ngm_vacc(cm, vacc_plan,sex=sex, ve=ve, R0_scale = scale)
    res$val
  }
  if(R0_doses(rep(thre, n.vac))>=1) return(rep(NA,n.vac))
  pop=rowSums(cm$misc$pop)[1:n.vac]
  full_list=c()
  m0=m
  while(R0_doses(init)>=1)
  {
    grads=grad(R0_doses, init)
    if(m0<=0) m0=m
    idx=order(grads/pop)
    i=idx[!idx%in%full_list][1]
    if((init[i]+m0/pop[i]*100)>=thre){
      m0=m0-(thre-init[i])*pop[i]/100
      init[i]=thre
      full_list=c(full_list,i)
    }else{
      init[i]=init[i]+m0/pop[i]*100
      m0=m
    }
  }
  init
}

#optimization while prioritizing high-sexual-active females
find_optim_plan_sex=function(cm,sex_age=c(20,30),thre=80,m=10, n.vac=8, scale=R0_scale){
  vacc_plan=vacc_rate(cm, age=sex_age, v_sex=rep(thre,length(sex_age)), v_all=rep(0,length(sex_age)))
  base_r=ngm_vacc(cm, vacc_plan,sex=T, ve=ve, R0_scale = scale)
  sex_idx=match(sex_age,cm$ageinterval)
  pop_all=rowSums(cm$misc$pop[1:n.vac,])
  pop_sex=(cm$misc$pop[,2]*cm$misc$paqa[[2]])[1:n.vac]
  pop_nonsex=pop_all-pop_sex*(c(1:n.vac)%in%sex_idx)
  v_sex=rep(0,n.vac)
  if(base_r$val<1){
    #vaccinating sex worker is sufficient
    R0_doses=function(v_rates){
      vacc_plan=vacc_rate(cm, age=sex_age, v_sex=v_rates, v_all=rep(0,length(sex_age)))
      res=ngm_vacc(cm, vacc_plan,sex=T, ve=ve, R0_scale = scale)
      res$val
    }
    init=rep(0,length(sex_age))
    pop=pop_sex[sex_idx]
    full_sex=F
    m=m/20
  }else{
    v_sex[sex_idx]=thre
    R0_doses=function(v_rates){
      vacc_plan=vacc_rate(cm, age=cm$ageinterval[1:n.vac], v_sex=v_sex, v_all=v_rates)
      res=ngm_vacc(cm, vacc_plan,sex=T, ve=ve, R0_scale = scale)
      res$val
    }
    if(R0_doses(rep(thre, n.vac))>=1) return(rep(NA,n.vac*2))
    init=rep(0,n.vac)
    pop=pop_nonsex
    full_sex=T
  }
  
  full_list=c()
  m0=m
  while(R0_doses(init)>=1)
  {
    #print(init)
    grads=grad(R0_doses, init)
    if(m0<=0) m0=m
    idx=order(grads/pop)
    i=idx[!idx%in%full_list][1]
    if((init[i]+m0/pop[i]*100)>=thre){
      m0=m0-(thre-init[i])*pop[i]/100
      init[i]=thre
      full_list=c(full_list,i)
    }else{
      init[i]=init[i]+m0/pop[i]*100
      m0=m
    }
  }
  if(full_sex){
    return(c(v_sex, init))
  }else{
    v_sex[sex_idx]=init
    return(c(v_sex, rep(0, n.vac)))
  }
}

#optimization without children
find_optim_plan_child=function(cm,thre=80,m=10, sex=T, n.vac=8, adult_age=20,scale=R0_scale){
  vac.start=match(adult_age, cm$ageinterval)
  init=rep(0,n.vac-vac.start+1)
  R0_doses=function(v_rates){
    vacc_plan=vacc_rate(cm, age=cm$ageinterval[vac.start:n.vac], v_sex=rep(0,n.vac-vac.start+1), v_all=v_rates)
    res=ngm_vacc(cm, vacc_plan,sex=sex, ve=ve, R0_scale = scale)
    #max(res$val,1)
    res$val
  }
  #grads=grad(R0_doses, init)
  if(R0_doses(rep(thre, n.vac))>=1) return(rep(NA,n.vac))
  pop=rowSums(cm$misc$pop)[vac.start:n.vac]
  full_list=c()
  m0=m
  while(R0_doses(init)>=1)
  {
    grads=grad(R0_doses, init)
    if(m0<=0) m0=m
    idx=order(grads/pop)
    i=idx[!idx%in%full_list][1]
    if((init[i]+m0/pop[i]*100)>=thre){
      m0=m0-(thre-init[i])*pop[i]/100
      init[i]=thre
      full_list=c(full_list,i)
    }else{
      init[i]=init[i]+m0/pop[i]*100
      m0=m
    }
  }
  c(rep(0,vac.start-1),init)
}

#fixed vaccination plans (mass)-----------
#calculate vaccination rates
rate_cal=function(age, all_p, cm, max_p=80, sequential=T, convert_to_plan=T){
  pop=rowSums(cm$misc$pop)
  doses=all_p/100*sum(pop)
  idx=match(age, cm$ageinterval)
  rates=rep(0, length(cm$ageinterval))
  if(sequential){
    i=1
    while(doses>0 & i<=length(idx)){
      rate1=min(max_p,doses/pop[idx[i]]*100)
      rates[idx[i]]=rate1
      doses=doses-rate1*pop[idx[i]]/100
      i=i+1
    }
  }else{
    rates[idx]=min(max_p, doses/sum(pop[idx])*100)
  }
  if(convert_to_plan){
    vacc_plan=list( #vaccine uptake rate by age group
      ageinterval=cm$ageinterval,
      vacc_sex=rates, 
      vacc_allm=rates,
      vacc_allf=rates 
    ) 
    return(vacc_plan)
  }else{
    return(rates)
  }
}

#calculate Reff for different simulations, fixed year, country, scenario
fix_vacc=function(year, country, vacc_age, p_vacc, thre=80, n_sim=n.sim, sequen=T, sex_p=T){
  cm=cm_setup(year, country, paras[1,],sex=max(sex_p,0.01))
  vacc_plan=rate_cal(vacc_age,p_vacc,cm,max_p=thre,sequential = sequen)
  mclapply(1:n_sim, function(i){
    par=paras[i,]
    cm=cm_setup(year, country, paras[i,],sex=max(sex_p,0.01))
    res=ngm_vacc(cm, vacc_plan,sex=(sex_p>0), ve=ve, R0_scale = R0_scales[i])
    res$val
  },mc.cores = numCores,mc.cleanup = TRUE)%>%unlist()
}





