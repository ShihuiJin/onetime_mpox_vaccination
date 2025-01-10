#vaccination effect estimations
ve=0.8; p_max=80 #can vary

#baseline scenario: no vaccination-------------
base_R0=list()
for(k in seq_along(scaling)) #levels of sexual transmission
{
  k0=scaling[k]
  print(paste0(k,': ',Sys.time()))
  base_R0[[k]]=lapply(countries, function(c){
    #vacc_c=coverage_birth[[c]]%>%colMeans()
    lapply(trange, function(t){
      mclapply(1:n.sim, function(i){
        cm=cm_setup(t,c,xy=paras[i,],sex=k0)
        res=ngm_vacc(cm,sex=(k0>0), ve=ve, R0_scale = R0_scales[i])
        c(res$val,res$p_sex)
      },mc.cores = numCores)%>%do.call('rbind',.)%>%as.matrix
    })
  })
}


#mass vaccination optimization------------
ndose_opt=list()
for(k in seq_along(scaling))
{
  k0=scaling[k]
  print(paste0(k,': ',Sys.time()))
  ndose_opt[[k]]=lapply(range(trange), function(t){
    lapply(countries, function(c){
      out=mclapply(1:n.sim, function(i){
        cm=cm_setup(t, c, paras[i,],sex=max(k0,0.01))
        m=10^(ceiling(log(sum(cm$misc$pop))/log(10))-3) #for a finer resolution can change 3 to 4
        find_optim_plan(cm,thre=p_max,m=m,sex=(k0>0),n.vac=8,scale=R0_scales[i])
      },mc.cores = numCores)%>%do.call('rbind',.)%>%as.matrix()
      if(any(is.na(out))) print(paste0('not reached: scale',k0, ';  time: ', t, ':  country: ', c))
      out
    })
  })
}

#prioritising sex workers--------
ndose_opt_sex=list()
for(k in 2:length(scaling)-1)
{
  k0=scaling[k+1]
  print(paste0(k,': ',Sys.time()))
  ndose_opt_sex[[k]]=lapply(c(2025,2050), function(t){
    lapply(countries, function(c){
      out=mclapply(1:n.sim, function(i){
        cm=cm_setup(t, c, paras[i,],sex=max(k0,0.01))
        m=10^(ceiling(log(sum(cm$misc$pop))/log(10))-3) #for a finer resolution can change 3 to 4
        find_optim_plan_sex(cm,sex_age = c(15,2:4*10),thre=p_max,m=m,n.vac=8, scale=R0_scales[i])
      },mc.cores = numCores)%>%do.call('rbind',.)%>%as.matrix()
      if(any(is.na(out))) print(paste0('not reached: scale',k0, ';  time: ', t, ':  country: ', c))
      out
    })
  })
}

#excluding young children-----------
ndose_opt_child=list()
for(k in seq_along(scaling))
{
  k0=scaling[k]
  print(paste0(k,': ',Sys.time()))
  ndose_opt_child[[k]]=lapply(range(trange), function(t){
    lapply(countries, function(c){
      out=mclapply(1:n.sim, function(i){
        cm=cm_setup(t, c, paras[i,],sex=max(k0,0.01))
        m=10^(ceiling(log(sum(cm$misc$pop))/log(10))-3) #for a finer resolution can change 3 to 4
        find_optim_plan_child(cm,thre=p_max,m=m,sex=(k0>0),n.vac=8,adult_age = 20,scale=R0_scales[i])
      },mc.cores = numCores)%>%do.call('rbind',.)%>%as.matrix()
      if(any(is.na(out))) print(paste0('not reached: scale',k0, ';  time: ', t, ':  country: ', c))
      out
    })
  })
}


#fixed mass vaccination strategies------------
fixed_res=list()
for(k in seq_along(scaling)) #levels of sexual transmission
{
  k0=scaling[k]
  print(paste0(k,': ',Sys.time()))
  fixed_res[[k]]=lapply(1:n.fixed, function(m){ 
    lapply(countries, function(c){
      lapply(mass_rate, function(p){
        lapply(trange, function(t){
          fix_vacc(t,c,fixed_list[[m]],p,p_max,n.sim,sequen = (m<=(n.fixed/2)), sex_p = k0)
        })
      })
    })
  })
}
