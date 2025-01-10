#data preparation
library(dplyr)
library(stringr)
library(purrr)
library(hhh4contacts)
library(parallel)
library(numDeriv) #numerical derivative
numCores <- detectCores()

#contact matrix
load('data/contact_home.RData')

#population
proj_pop=read.csv('data/pop_proj_wb.csv')

#countries
load('data/country_list.RData')

#high-sexual-activity population
sex_data=read.csv('data/African_sex_prop.csv')

#vaccine coverage
load('data/coverage_birth_mean.RData')

#parameter estimates (1000 draws)
load('data/fitted_DRC_1k_boot.RData')

#countries of interest
countries=country_list%>%
  filter(region!='Northern Africa'|country=='Sudan',ISO%in%sex_data$ISO)%>%
  pull(ISO)

#vaccination-related settings (invariant)
n.sim=1000
trange=seq(2025,2050,5)
scaling=c(0,0.25,0.5,1)
mass_rate=c(5,10,20,30)
n.fixed=10
fixed_list=list(
  #sequential
  c(0,5,10),
  c(0,20,30),
  c(20,0,30),
  c(20,15,30),
  c(20,30,40),
  #average
  c(0,5,10),
  c(0,20),
  c(20,30),
  c(20,30,40), 
  c(0,5,10,15,20,30,40)
)


