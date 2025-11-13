rm(list=ls())
source('01_functions_S5.R')#converting mutualism matrix to ones and zeros
require(deSolve) ## for integrating ordinary differential equations
require(tidyverse) ## for efficient data manipulation & plotting
require(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(ggplot2)
library(viridis)
library(igraph)
library(bipartite)


theme_set(theme_bw()) 


#dat <- paste('network.csv',sep='')
adj.mat<-function(data){
  d <- read.csv(file=data,header=FALSE )
  dat<-as.matrix(d)
  dat[dat > 0] = 1
  dat[dat < 0] = 
    dat<-apply(dat,2,as.numeric)
  return(dat)}

mydir = 'datasets_1'

#-------------------------------------------------------------------------------
webfiles = list.files(path = mydir, pattern = "*.csv", full.names = TRUE)



#print(r)
g<-adj.mat(webfiles[which(webfiles =="datasets_1/M_PL_061_09.csv")]) #ork web names
Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.animals<-degree.plants<-numeric()

#degree of plants and anichmals
for(i in 1:Plantspecies){
  degree.plants[i]<-sum(g[i,])} # degree of plants
for(j in 1:Aspecies){
  degree.animals[j]<-sum(g[,j]) # degree of animals
}

degree<-c(degree.animals, degree.plants)

#parameters for modelling: intialisation of the model

## vector of species trait standard deviations

N <- runif( (Aspecies+Plantspecies) , 1,1)  ## initial species densities
#muA<- #runif(Aspecies, 12, 30)  #initial mean phenotypic optimum trait values
#muP<- #runif(Plantspecies, 15, 30)   #intial mean phenotypic optimum trait values
na <-  runif( (Aspecies) , 1,1) # (init_data %>% filter(webname == as.character(fact$web[r])))$density[1:Aspecies]
np<- runif( (Plantspecies) , 1,1)# (init_data %>% filter(webname == fact$web[r]))$density[ (Aspecies+1):(Aspecies+Plantspecies)]

muA<- rnorm(Aspecies, mean = 22, sd= 2) #initial mean phenotypic optimum trait values
muP<- rnorm(Plantspecies, mean = 22, sd= 2)


#incorporating NOISE in the dynamics
tmax<-500
time_range<-c(0,tmax)
duration <-500 #fact$forcing_duration[r]
d<- c(rep(1,duration),rep(0,(deltat-duration)))
deltat<- (time_range[2]-time_range[1])/1 + 1
d<- c(rep(1,duration),rep(0,(deltat-duration)))
temp<-22 #Temperature
sd_temp<-2
noise_c_d<-replicate(1,rnorm(d,temp,sd_temp))
times<-seq(0, tmax, 1)
noise_T<-as.data.frame(list(times=times,import=rep(0,length(times))))
noise_T$import <- noise_c_d

temp<-22 #fact$Temperature[r] #20 degrees
nestedness<-nestedness_NODF(g)
C<-Connectance(g)
web.name<-"datasets_1/M_PL_061_09.csv"
dganimals<-degree.animals
dgplants<-degree.plants
sig <-runif((Aspecies+Plantspecies),0.009,0.03)
envA<-runif(Aspecies,0.009, 0.03)
envP<- runif(Plantspecies, 0.009, 0.03)
sa<- sig[1:Aspecies]
sp<-sig[(Aspecies+1): (Aspecies+Plantspecies)]
pa<-envA+sa
pp<-envP+sp
trait_evolve<-1 #no evoluation of mean trait
var_evolve<-1 #no evolution of variance 

t1<-plot_snapshot(Na=na, Np=np, m = c(muA,muP), 
                  sigma = c(pa,pp), Temp=22, limits = c(20,25),res = 1001)+
  annotate("text",x=22, y =5, label="time == 2000",parse=T)

t1


bw  <-4
aw  <- 0
gi <- 1.5
ki <-0.4#mortality rate
w<- 1
ic <-c(na, np, muA, muP, sa , sp )
mut.strength<-1
species_index <- which(degree==max(degree))
params <- list(time=time,matrix=g,Temp=22,tmax=tmax,sig=sig,bw=bw,gi=gi,ki=ki,variation=fact$variation[1],
               web.name=web.name,w=w,trait_evolve=trait_evolve,var_evolve=1,
               mut.strength=mut.strength,C=C,nestedness=nestedness,h2="evo_trait_variance",
               web.name=web.name,dganimals=dganimals,degree=degree,w_c=0.02,
               dgplants=dgplants,envA=envA,envP=envP,noise_T=noise_T)

species_index <- which(degree==max(degree))

OUT<-ode(func=eqs_type2_new, y=ic, parms=params,
         times=seq(0, tmax, by=20)) %>% 
  organize_results(pars = params) 

OUT %>% plot_density()

sol_1<-  OUT %>% filter(time == tmax)
Na1<- (sol_1 %>% filter(type=="N"))$v
Np1<-(sol_1 %>% filter(type =="P"))$v
ma1<-(sol_1 %>% filter(type == "muA"))$v
mp1<-(sol_1 %>% filter(type == "muP"))$v
sa1<-(sol_1 %>% filter(type == "sA"))$v
sp1<-(sol_1 %>% filter(type == "sP"))$v



#@@@@@@@@############################ Environmental change of 2 degrees: ########################

tmax<-500
time_range<-c(0,tmax)
duration <-500 #fact$forcing_duration[r]
d<- c(rep(1,duration),rep(0,(deltat-duration)))
deltat<- (time_range[2]-time_range[1])/1 + 1
d<- c(rep(1,duration),rep(0,(deltat-duration)))
temp<-25 #fact$Temperature[r] #20 degrees
sd_temp<-2
noise_c_d<-replicate(1,rnorm(d,temp,sd_temp))
times<-seq(0, tmax, 1)
noise_T<-as.data.frame(list(times=times,import=rep(0,length(times))))
noise_T$import <- noise_c_d


bw  <-4
aw  <- 0
gi <- 1.25
ki <-0.4#mortality rate
w<- 1
ic <-c(Na1, Np1, ma1, mp1, sa1 , sp1 )
mut.strength<-1
species_index <- which(degree==max(degree))
params <- list(time=time,matrix=g,Temp=24,tmax=tmax,sig=sig,bw=bw,gi=gi,ki=ki,variation=fact$variation[1],
               web.name=web.name,w=w,trait_evolve=trait_evolve,var_evolve=1,
               mut.strength=mut.strength,C=C,nestedness=nestedness,h2="evo_trait_variance",
               web.name=web.name,dganimals=dganimals,degree=degree,w_c=0.02,
               dgplants=dgplants,envA=envA,envP=envP,noise_T=noise_T)

species_index <- which(degree==max(degree))

OUT2<-ode(func=eqs_type2_new, y=ic, parms=params,
         times=seq(0, tmax, by=10)) %>% 
  organize_results(pars = params) 

OUT2 %>% plot_density()

sol_2<-  OUT2 %>% filter(time == tmax)
Na2<- (sol_2 %>% filter(type=="N"))$v
Np2<-(sol_2 %>% filter(type =="P"))$v
ma2<-(sol_2 %>% filter(type == "muA"))$v
mp2<-(sol_2 %>% filter(type == "muP"))$v
sa2<-(sol_2 %>% filter(type == "sA"))$v
sp2<-(sol_2 %>% filter(type == "sP"))$v



