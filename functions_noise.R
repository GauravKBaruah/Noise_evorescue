#rm(list=ls())

require(statmod)
cutoff <- function(x) ifelse(x<1, (1*(x>0))*(x*x*x*(10+x*(-15+6*x))), 1)
adj.mat<-function(data){
  #dat <- paste('network.csv',sep='')
  d <- read.csv(file=data,header=FALSE )
  dat<-as.matrix(d)
  dat[dat > 0] = 1
  dat<-apply(dat,2,as.numeric)
  return(dat)}


min_no_path_to_target_node<- function(graph_matrix, start_node, target_node ){
  
  
  g1 <- graph_from_biadjacency_matrix(graph_matrix, weighted = T)
  #all_paths <- all_simple_paths(g1, from = start_node, to = target_node, mode = "out")
  shortest_path<- length((shortest_paths(g1, from = start_node, to = target_node))$vpath[[1]])-1
  
  
  #paths<-numeric()
  #if(length(all_paths) == 0){
  # paths <- NA
  #}else{ 
  # for(i in 1:length(all_paths)){
  #paths[i]<- length(all_paths[[i]])-1
  
  #}
  return(shortest_path)
}

gausquad.animals<-function(m,sigma,w,h,np,na,mut.strength,
                           points,mat,degree.animal,P){
  
  
  temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
  if(mat == 0){
    return(list(G= 0, B = 0,V=0,Ja=0 ))
  }
  else if(mat == 1 ){
    #nodes oir points in the abscissa where the integral will be evaluated numerically
    
    z1<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sqrt(sigma$sa))$nodes #z'
    z2<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sqrt(sigma$sp))$nodes #z''
    
    #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
    w1<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$ma,sigma =sqrt(sigma$sa))$weights #pi(z')
    w2<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$mp,sigma =sqrt(sigma$sp))$weights #pj(z'')
    
    
    #for the pairwise model however there are only two species interacting and hence i and j
    #or in other words the integral goes over z and z'
    for (i in 1: points){
      
      temp2[i]<-sum(np*exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]/(1+h*exp(-(z1[i]- z2)^2/w^2)*np))
      f <-  exp(-(z1[i]- z2)^2/w^2) 
      j1[i]<- sum(na*(mut.strength)*f/(1+h*np*(mut.strength/degree.animal)*f)^2*w2*w1[i])
      gvar[i]<-1/P^2*sum(np*(mut.strength/degree.animal)*((z1[i]-m$ma)^2-P)*exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]/(1+h*exp(-(z1[i]- z2)^2/w^2)*np))
      dat2[i]<- sum(np*(mut.strength/degree.animal)*(z1[i]-m$ma)*(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i])/(1+h*exp(-(z1[i]- z2)^2/w^2)*np))
      
      
    }
    G = sum(temp2)
    B = sum(dat2) 
    V = sum(gvar)
    
    
    
    return(list(G= G, B = B,V=V, Ja=sum(j1)))
  }else{
    
    z1<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sqrt(sigma$sa))$nodes #z'
    z2<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sqrt(sigma$sp))$nodes #z''
    
    #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
    w1<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$ma,sigma =sigma$sa)$weights #pi(z')
    w2<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$mp,sigma =sigma$sp)$weights #pj(z'')
    
    
    #for the pairwise model however there are only two species interacting and hence i and j
    #or in other words the integral goes over z and z'
    for (i in 1: points){
      
      temp2[i]<-sum(np*exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]/(1+h*exp(-(z1[i]- z2)^2/w^2)*np))
      f <-  exp(-(z1[i]- z2)^2/w^2) 
      j1[i]<- sum(na*(mut.strength)*f/(1+h*np*(mut.strength/degree.animal)*f)^2*w2*w1[i])
      gvar[i]<-1/P^2*sum(np*(mut.strength/degree.animal)*((z1[i]-m$ma)^2-P)*exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]/(1+h*exp(-(z1[i]- z2)^2/w^2)*np))
      dat2[i]<- sum(np*(mut.strength/degree.animal)*(z1[i]-m$ma)*(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i])/(1+h*exp(-(z1[i]- z2)^2/w^2)*np))
      
      
    }
    G = sum(temp2)
    B = sum(dat2) 
    V = sum(gvar)
    
    
    
    return(list(G= G, B = B,V=V, Ja=sum(j1)))
    
  }
}


gausquad.plants<-function(m,sigma,w,h,np,na,mut.strength,points,mat,degree.plant,P){
  
  temp2<-dat2<-x3<-x4<-gvar<-j1<-array(dim=c(points))
  
  if(mat == 0){
    
    return(list(G= 0, 
                B = 0,V=0,
                Jp=0))
    
  }
  else if (mat==1){
    
    #nodes oir points in the abscissa where the integral will be evaluated numerically
    z1<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sqrt(sigma$sp))$nodes #z'
    z2<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sqrt(sigma$sa))$nodes #z''
    
    #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
    w1<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$mp,sigma =sqrt(sigma$sp))$weights #pi(z')
    w2<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$ma,sigma =sqrt(sigma$sa))$weights #pj(z'')
    
    
    #for the pairwise model however there are only two species interacting and hence i and j
    #or in other words the integral goes over z and z'
    for (i in 1: points){
      f <-  exp(-(z1[i]- z2)^2/w^2) 
      temp2[i]<- sum(na*exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]/(1+h*exp(-(z1[i]- z2)^2/w^2)*na))
      dat2[i]<- sum(na*(mut.strength/degree.plant)*(z1[i]-m$mp)*(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i])/(1+h*exp(-(z1[i]- z2)^2/w^2)*na))
      gvar[i]<-1/P^2*sum(na*(mut.strength/degree.plant)*((z1[i]-m$mp)^2-P)*exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]/(1+h*exp(-(z1[i]- z2)^2/w^2)*na))
      j1[i]<-sum(np*(mut.strength/degree.plant)*f/(1+h*na*(mut.strength/degree.plant)*f)^2*w2*w1[i])
      
    }
    
    G = sum(temp2)#sum(temp2*na*(mut.strength/degree.plant)/(1+x3*(mut.strength/degree.plant)*na))http://127.0.0.1:8283/graphics/a77cf26d-25d3-4aa0-858d-ec22c28b9498.png
    B = sum(dat2)#na*(mut.strength/degree.plant)/(1+x3*(mut.strength/degree.plant)*na))http://127.0.0.1:8283/graphics/a77cf26d-25d3-4aa0-858d-ec22c28b9498.png
    V = sum(gvar) #sum(gvar*na*(mut.strength/degree.plant)/(1+x3*(mut.strength/degree.plant)*na))
    
    
    return(list(G= G, 
                B = B,
                V=V,
                Jp=sum(j1)))
  }
  
}


eqs <- function(time, state, pars) {
  A <- dim(pars$matrix)[2]  ## number of animal species
  P <-dim(pars$matrix)[1]
  a <- state[1:A] ## species densities of animals
  p <- state[(A+1):(A+P)] ## species densities of plants
  s <- pars$sig ## species' trait standard deviations
  ma<-state[(A+P+1):(A+P+A)]
  mp<-state[(A+P+A+1):(A+P+A+P)]
  ## define g, where g[i] is the selection pressure on species i from growth
  alphaA<-pars$Amatrix ## alpha matrix
  alphaP<-pars$Pmatrix
  dt<-0.01
  
  aij<-bij<-vij<-matrix(0, nrow=A,ncol=P) 
  aji<-bji<-vji<-matrix(0, nrow=P,ncol=A) 
  Na<-muA<- Pa <- vA<-matrix(0, nrow = time, ncol = A )
  Np<-muP<-Pp<-vP<-matrix(0, nrow =time, ncol = P )
  Np[1,]<-state[(A+1):(A+P)]
  Na[1,]<-state[1:A]
  muA[1,]<-ma
  muP[1,]<-mp
  vA[1,]<-s[1:A]
  vP[1,]<-s[(A+1):(A+P)]
  Pa[1,]<-vA[1,]+pars$envA[1:A]
  Pp[1,]<-vP[1,]+ pars$envP
  aj<-bj<-ai<-bi<-vi<-vj<-xar_bp<-xar_ba<-ba<-bp<-bar_bp<-bar_ba<-numeric()
  Temp<- rnorm(time, 20, 4)
  for (t in 1:(time-1)){
    
    
    for(r in 1:A){
      ba[r]<- pars$gi/(pars$bw)*(pars$bw)/(sqrt((pars$bw)^2+Pa[t,r]))*exp(-(Temp[t]- muA[t,r])^2/(2*(pars$bw)^2+Pa[t,r])) - pars$ki
      bar_ba[r]<- pars$gi/(pars$bw)*exp(-(Temp[t]- muA[t,r])^2/(2*(pars$bw)^2+Pa[t,r]))*sqrt(Pa[t,r])*(pars$bw)*(Temp[t] -muA[t,r])/((pars$bw)^2 + Pa[t,r])^1.5
      xar_ba[r]<- (2*pars$gi*(2*muA[t,r]^2-4*Temp[t]*muA[t,r]+2*Temp[t]^2-2*pars$bw^2-Pa[t,r]^2)*exp(-(Temp[t]-muA[t,r])^2/(2*pars$bw^2+Pa[t,r]^2)))/(sqrt(pars$bw^2+Pa[t,r]^2)*(2*pars$bw^2+Pa[t,r]^2)^2)
      
    }
    for(j in 1:P){
      bp[j]<- pars$gi/(pars$bw)*(pars$bw)/(sqrt((pars$bw)^2+ Pp[t,j]))*exp(-(Temp[t]- muP[t,j])^2/(2*(pars$bw)^2+Pp[t,j])) - pars$ki
      bar_bp[j]<-pars$gi/(pars$bw)*exp(-(Temp[t]- muP[t,j])^2/(2*(pars$bw)^2+Pp[t,j]))*sqrt(Pp[t,j])*(pars$bw)*(Temp[t] - muP[t,j])/((pars$bw)^2 + Pp[t,j])^1.5 
      xar_bp[j]<- (2*pars$gi*(2*muP[t,j]^2-4*Temp[t]*muP[t,j]+2*Temp[t]^2-2*pars$bw^2-Pp[t,j]^2)*exp(-(Temp[t]-muP[t,j])^2/(2*pars$bw^2+Pp[t,j]^2)))/(sqrt(pars$bw^2+Pp[t,j]^2)*(2*pars$bw^2+Pp[t,j]^2)^2)
      
    }
    
    for(r in 1:A){
      for(l in 1:P){
        #
        m.temp<-list(ma=muA[t,r],mp=muP[t,l])
        sigma1<-list(sa=Pa[t,r],sp=Pp[t,l])
        temp1<-gausquad.animals(m=m.temp,sigma=sigma1,w=pars$w,h=0.25,np=Np[t,l],na=Na[t,r],
                                mut.strength=pars$mut.strength, points=6
                                ,mat=pars$matrix[l,r],
                                degree.animal = pars$dganimals[r],
                                P=Pa[t,r])
        aij[r,l]<-temp1$G
        bij[r,l]<-temp1$B
        vij[r,l]<-temp1$V
        
      }
      ai[r]<-sum(aij[r,])
      bi[r]<-sum(bij[r,])
      vi[r]<-sum(vij[r,])
    }
    for(k in 1:P){
      for(m in 1:A){
        m2.temp<-list(ma=muA[t,m],mp=muP[t,k])
        sigma2<-list(sa=Pa[t,m],sp=Pp[t,k])
        temp2<-gausquad.plants(m=m2.temp,sigma=sigma2,w=pars$w,h=0.25,np=Np[t,k],na=Na[t,m],
                               mut.strength=pars$mut.strength,
                               points=6,mat=pars$matrix[k,m], 
                               degree.plant =pars$dgplants[k],
                               P=Pp[t,k])
        aji[k,m] <-temp2$G
        bji[k,m]<-temp2$B
        vji[k,m]<-temp2$V
      }
      aj[k]<-sum(aji[k,])
      bj[k]<-sum(bji[k,])
      vj[k]<-sum(vji[k,])
    }
    #print(t)
    
    
    
    Na[t+1, ]<-Na[t,] + Na[t,]*(ba-alphaA%*%Na[t,]+ai)*dt #+ rnorm(A, 0,sd=0.05)*Na[t,]*dt## density eqs
    Np[t+1, ]<-Np[t,] + Np[t,]*(bp-alphaP%*%Np[t,]+aj)*dt #+ rnorm(P, 0,sd=0.05)*Np[t,]*dt ## trait mean eqs
    muA[t+1, ]<-muA[t,] +(vA[t,]/Pa[t,])*(bar_ba+ bi)*dt #+ rnorm(A, 0,sd=0)*dt ## trait mean eqs
    muP[t+1, ]<- muP[t,]+(vP[t,]/Pp[t,])*(bar_bp+ bj)*dt #+ rnorm(P, 0,sd=0)*dt ## trait mean eqs
    vA[t+1,]<-vA[t,] + 0.5*(vA[t,]/Pa[t,])^2*(xar_ba + vi)*dt*cutoff(vA[t,]/(1e-10))
    vP[t+1,]<-vP[t,] + 0.5*(vP[t,]/Pp[t,])^2*(xar_bp + vj)*dt*cutoff(vP[t,]/(1e-10))
    
    #+ 
    #
    Pa[t+1, ] <- vA[t+1,] + pars$envA
    Pp[t+1, ] <- vP[t+1,] + pars$envP
    
    Na[t+1,which(Na[t+1,] < 0)]<-0
    Np[t+1,which(Np[t+1,] < 0)]<-0
    
    
    
  } 
  ## return equations by first flattening them back into a single vector
  output= list(Plants = Np[1:time,],Animals=Na[1:time,], Plant.trait = muP[1:time,], 
               Animal.trait=muA[1:time,],Pa=Pa[1:time,], Pp=Pp[1:time,], aij=aij,ba=pars$ba)
  return(output)
}

Mcommunity = function(iter, time, ...){
  set.seed(rnorm(1,as.numeric(Sys.time())-Sys.getpid(),10000)) 
  init = time
  replicate =try(eqs(time=init, ...))
  replicate$start = init
  replicate$iter = iter
  return(replicate)
}




eqs_type2 <- function(time, state, pars) {
  A <- dim(pars$matrix)[2]  ## number of animal species
  P <-dim(pars$matrix)[1]
  Np<-state[(A+1):(A+P)]
  Na<-state[1:A]
  s <- pars$sig ## species' trait standard deviations
  ma<-state[(A+P+1):(A+P+A)]
  mp<-state[(A+P+A+1):(A+P+A+P)]
  gvarA<-state[(A+P+A+P+1):(A+P+A+P+A)]
  gvarP<-state[(A+P+A+P+A+1):(A+P+A+P+A+P)]
  
  s <- pars$sig ## species' trait standard deviations
  alpha.a<-pars$Amatrix ## alpha matrix
  alpha.p<-pars$Pmatrix
  
  aij<-bij<-vij<-matrix(0, nrow=A,ncol=P) 
  aji<-bji<-vji<-matrix(0, nrow=P,ncol=A) 
  muA<-ma
  muP<-mp
  aj<-bj<-ai<-bi<-vj<-vi<-numeric()
  #*cutoff(gvarA/(1e-7))
  varA <- ( (gvarA + pars$envA))
  varP <- ( (gvarP + pars$envP))
  
  # noise_f_P<-approxfun(pars$noise_T$import[,1],method="linear", rule=2)
  Temp<- pars$Temp#replicate(1,noise_f_P(time))
  #pars$Temp #
  ba<- pars$gi/(pars$bw)*(pars$bw)/(sqrt((pars$bw)^2+varA))*exp(-(Temp- muA)^2/(2*(pars$bw)^2+varA)) - pars$ki
  bar_ba<- pars$gi/(pars$bw)*exp(-(Temp- muA)^2/(2*(pars$bw)^2+varA))*(varA)*(pars$bw)*(Temp -muA)/((pars$bw)^2 + varA)^1.5
  xar_ba<- (2*pars$gi*(2*muA^2-4*Temp*muA+2*Temp^2-2*pars$bw^2-varA)*exp(-(Temp-muA)^2/(2*pars$bw^2+varA)))/(sqrt(pars$bw^2+varA)*(2*pars$bw^2+varA)^2)
  
  bp<- pars$gi/(pars$bw)*(pars$bw)/(sqrt((pars$bw)^2+ varP))*exp(-(Temp- muP)^2/(2*(pars$bw)^2+varP)) - pars$ki
  bar_bp<-pars$gi/(pars$bw)*exp(-(Temp- muP)^2/(2*(pars$bw)^2+varP))*(varP)*(pars$bw)*(Temp - muP)/((pars$bw)^2 + varP)^1.5 
  xar_bp<- (2*pars$gi*(2*muP^2-4*Temp*muP+2*Temp^2-2*pars$bw^2-varP)*exp(-(Temp-muP)^2/(2*pars$bw^2+varP)))/(sqrt(pars$bw^2+varP)*(2*pars$bw^2+varP)^2)
  
  for(r in 1:A){
    for(l in 1:P){
      #
      m.temp<-list(ma=muA[r],mp=muP[l])
      sigma1<-list(sa=varA[r],sp=varP[l])
      temp1<-gausquad.animals(m=m.temp,sigma=sigma1,w=pars$w,h=0.2,np=Np[l],na=Na[r],
                              mut.strength=pars$mut.strength, points=5
                              ,mat=pars$matrix[l,r],
                              degree.animal = pars$dganimals[r],
                              P=varA[r])
      aij[r,l]<-temp1$G
      bij[r,l]<-temp1$B
      vij[r,l]<-temp1$V
      
    }
    ai[r]<-sum(aij[r,])
    bi[r]<-sum(bij[r,])
    vi[r]<-sum(vij[r,])
  }
  for(k in 1:P){
    for(m in 1:A){
      m2.temp<-list(ma=muA[m],mp=muP[k])
      sigma2<-list(sa=varA[m],sp=varP[k])
      temp2<-gausquad.plants(m=m2.temp,sigma=sigma2,w=pars$w,h=0.2,np=Np[k],na=Na[m],
                             mut.strength=pars$mut.strength,
                             points=5,mat=pars$matrix[k,m], 
                             degree.plant =pars$dgplants[k],
                             P=varP[k])
      aji[k,m] <-temp2$G
      bji[k,m]<-temp2$B
      vji[k,m]<-temp2$V
    }
    aj[k]<-sum(aji[k,])
    bj[k]<-sum(bji[k,])
    vj[k]<-sum(vji[k,])
  }
  
  dndt_a<- Na*(ba-pars$Amatrix%*%Na+ai)*cutoff(Na/(1e-6)) #  na*(ba-alpha.a%*%na+ai)*cutoff(na/(1e-8)) #population dynamics
  dndt_p<- Np*(bp-pars$Pmatrix%*%Np+aj)*cutoff(Np/(1e-6))  #population dynamics
  dudt_A<- pars$trait_evolve*gvarA/varA*(bar_ba+ bi) #mean trait dynamics
  dudt_P<- pars$trait_evolve*gvarP/varP*(bar_bp+ bj) #mean trait dynamics
  dudt_vA<- pars$var_evolve*1/2*(gvarA/varA)^2*(xar_ba + vi)*cutoff(gvarA/(1e-7))
  dudt_vP<- pars$var_evolve*1/2*(gvarP/varP)^2*(xar_bp + vj)*cutoff(gvarP/(1e-7))
  
  ## return equations by first flattening them back into a single vector
  return(list(c(dndt_a, dndt_p,dudt_A,dudt_P,dudt_vA,dudt_vP)))
}




eqs_type2_noise <- function(time, state, pars) {
  A <- dim(pars$matrix)[2]  ## number of animal species
  P <-dim(pars$matrix)[1]
  Np<-state[(A+1):(A+P)]
  Na<-state[1:A]
  s <- pars$sig ## species' trait standard deviations
  ma<-state[(A+P+1):(A+P+A)]
  mp<-state[(A+P+A+1):(A+P+A+P)]
  gvarA<-state[(A+P+A+P+1):(A+P+A+P+A)]
  gvarP<-state[(A+P+A+P+A+1):(A+P+A+P+A+P)]
  
  s <- pars$sig ## species' trait standard deviations
  alpha.a<-pars$Amatrix ## alpha matrix
  alpha.p<-pars$Pmatrix
  
  aij<-bij<-vij<-matrix(0, nrow=A,ncol=P) 
  aji<-bji<-vji<-matrix(0, nrow=P,ncol=A) 
  muA<-ma
  muP<-mp
  aj<-bj<-ai<-bi<-vj<-vi<-numeric()
  #*cutoff(gvarA/(1e-7))
  varA <- ( (gvarA + pars$envA))
  varP <- ( (gvarP + pars$envP))
  
  noise_f_P<-approxfun(pars$noise_T$import[,1],method="linear", rule=2)
  Temp<- replicate(1,noise_f_P(time))
  #pars$Temp #
  ba<- pars$gi/(pars$bw)*(pars$bw)/(sqrt((pars$bw)^2+varA))*exp(-(Temp- muA)^2/(2*(pars$bw)^2+varA)) - pars$ki
  bar_ba<- pars$gi/(pars$bw)*exp(-(Temp- muA)^2/(2*(pars$bw)^2+varA))*(varA)*(pars$bw)*(Temp -muA)/((pars$bw)^2 + varA)^1.5
  xar_ba<- (2*pars$gi*(2*muA^2-4*Temp*muA+2*Temp^2-2*pars$bw^2-varA)*exp(-(Temp-muA)^2/(2*pars$bw^2+varA)))/(sqrt(pars$bw^2+varA)*(2*pars$bw^2+varA)^2)
  
  bp<- pars$gi/(pars$bw)*(pars$bw)/(sqrt((pars$bw)^2+ varP))*exp(-(Temp- muP)^2/(2*(pars$bw)^2+varP)) - pars$ki
  bar_bp<-pars$gi/(pars$bw)*exp(-(Temp- muP)^2/(2*(pars$bw)^2+varP))*(varP)*(pars$bw)*(Temp - muP)/((pars$bw)^2 + varP)^1.5 
  xar_bp<- (2*pars$gi*(2*muP^2-4*Temp*muP+2*Temp^2-2*pars$bw^2-varP)*exp(-(Temp-muP)^2/(2*pars$bw^2+varP)))/(sqrt(pars$bw^2+varP)*(2*pars$bw^2+varP)^2)
  
  for(r in 1:A){
    for(l in 1:P){
      #
      m.temp<-list(ma=muA[r],mp=muP[l])
      sigma1<-list(sa=varA[r],sp=varP[l])
      temp1<-gausquad.animals(m=m.temp,sigma=sigma1,w=pars$w,h=0.1,np=Np[l],na=Na[r],
                              mut.strength=pars$mut.strength, points=5
                              ,mat=pars$matrix[l,r],
                              degree.animal = pars$dganimals[r],
                              P=varA[r])
      aij[r,l]<-temp1$G
      bij[r,l]<-temp1$B
      vij[r,l]<-temp1$V
      
    }
    ai[r]<-sum(aij[r,])
    bi[r]<-sum(bij[r,])
    vi[r]<-sum(vij[r,])
  }
  for(k in 1:P){
    for(m in 1:A){
      m2.temp<-list(ma=muA[m],mp=muP[k])
      sigma2<-list(sa=varA[m],sp=varP[k])
      temp2<-gausquad.plants(m=m2.temp,sigma=sigma2,w=pars$w,h=0.1,np=Np[k],na=Na[m],
                             mut.strength=pars$mut.strength,
                             points=5,mat=pars$matrix[k,m], 
                             degree.plant =pars$dgplants[k],
                             P=varP[k])
      aji[k,m] <-temp2$G
      bji[k,m]<-temp2$B
      vji[k,m]<-temp2$V
    }
    aj[k]<-sum(aji[k,])
    bj[k]<-sum(bji[k,])
    vj[k]<-sum(vji[k,])
  }
  
  dndt_a<- Na*(ba-pars$Amatrix%*%Na+ai)*cutoff(Na/(1e-6)) #  na*(ba-alpha.a%*%na+ai)*cutoff(na/(1e-8)) #population dynamics
  dndt_p<- Np*(bp-pars$Pmatrix%*%Np+aj)*cutoff(Np/(1e-6))  #population dynamics
  dudt_A<- pars$trait_evolve*gvarA/varA*(bar_ba+ bi) #mean trait dynamics
  dudt_P<- pars$trait_evolve*gvarP/varP*(bar_bp+ bj) #mean trait dynamics
  dudt_vA<- pars$var_evolve*1/2*(gvarA/varA)^2*(xar_ba + vi)*cutoff(gvarA/(1e-7))
  dudt_vP<- pars$var_evolve*1/2*(gvarP/varP)^2*(xar_bp + vj)*cutoff(gvarP/(1e-7))
  
  ## return equations by first flattening them back into a single vector
  return(list(c(dndt_a, dndt_p,dudt_A,dudt_P,dudt_vA,dudt_vP)))
}


#na, #np : vector of animal and plant densities
#a_index: animal species i
# w : gaussian width of interaction
#degree.animal: degree of the animal
# P : phenotypic variance

type_2_animals<-function(m,sigma,w,h,np,na,mut.strength,a_index,
                         points,mat,degree.animal,Pa){
  temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
  z1<-gauss.quad.prob(points, dist = "normal", mu=m$ma[a_index], sigma =sqrt(sigma$sa[a_index]))$nodes #z of animals
  w1<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$ma[a_index],sigma =sqrt(sigma$sa[a_index]))$weights #p_i(z)
  
  z2<-matrix(0,nrow=length(np),ncol=points)
  w2<-matrix(0,nrow=length(np),ncol=points)
  numer_a<-denom_a<-numer_m<-denom_m<-numer_g<-denom_g<-matrix(0,nrow=points,ncol=length(np))
  N_strength<-m_strength<-g_strength<-numeric()
  
  
  for(j in 1:points){  
    for(k in 1:length(np)){
      z2[k,]<-gauss.quad.prob(points, dist = "normal", mu=m$mp[k], sigma =sqrt(sigma$sp[k]))$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w2[k,]<-gauss.quad.prob(points, dist = "normal", 
                              mu=m$mp[k],sigma =sqrt(sigma$sp[k]))$weights #pj(z'')
      
      numer_a[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_a[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_m[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum((z1[j]-m$ma[a_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_m[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_g[j,k]<- np[k]*mat[k,a_index]*1/Pa[a_index]^2*(mut.strength/degree.animal)*sum(((z1[j]-m$ma[a_index])^2-Pa[a_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_g[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
    }
    
    N_strength[j] <-(sum(numer_a[j,])/(1+h*sum(denom_a[j,])))*w1[j]
    m_strength[j] <-(sum(numer_m[j,])/(1+h*sum(denom_m[j,])))*w1[j]
    g_strength[j] <-(sum(numer_g[j,])/(1+h*sum(denom_g[j,])))*w1[j]
  }
  
  
  #  temp2[i]<-sum(np*(mut.strength/degree.animal)*exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]/(1+h*sum(denom[i,])))
  # f <-  exp(-(z1[i]- z2)^2/w^2) 
  #j1[i]<- sum(na*(mut.strength)*f/(1+h*np*(mut.strength/degree.animal)*f)^2*w2*w1[i])
  #gvar[i]<-1/P^2*sum(np*(mut.strength/degree.animal)*((z1[i]-m$ma)^2-P)*exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]/(1+h*exp(-(z1[i]- z2)^2/w^2)*np))
  #dat2[i]<- sum(np*(mut.strength/degree.animal)*(z1[i]-m$ma[index])*(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i])/(1+h*sum(denom[i,])))
  
  
  
  G = sum(N_strength)
  B = sum(m_strength) 
  V = sum(g_strength)
  
  
  
  return(list(G= G, B = B,V=V))
  
}

type_2_plants<-function(m,sigma,w,h,np,na,mut.strength,p_index,
                        points,mat,degree.plant,Pa){
  temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
  
  z1<-gauss.quad.prob(points, dist = "normal", mu=m$mp[p_index], sigma =sqrt(sigma$sp[p_index]))$nodes #z of animals
  w1<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$mp[p_index],sigma =sqrt(sigma$sp[p_index]))$weights #p_i(z)
  
  z2<-matrix(0,nrow=length(na),ncol=points)
  w2<-matrix(0,nrow=length(na),ncol=points)
  numer_a<-denom_a<-numer_m<-denom_m<-numer_g<-denom_g<-matrix(0,nrow=points,ncol=length(na))
  N_strength<-m_strength<-g_strength<-numeric()
  
  for(j in 1:points){  
    for(k in 1:length(na)){
      z2[k,]<-gauss.quad.prob(points, dist = "normal", mu=m$ma[k], sigma =sqrt(sigma$sa[k]))$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w2[k,]<-gauss.quad.prob(points, dist = "normal", 
                              mu=m$ma[k],sigma =sqrt(sigma$sa[k]))$weights #pj(z'')
      
      numer_a[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_a[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_m[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum((z1[j]-m$mp[p_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_m[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_g[j,k]<- na[k]*mat[p_index,k]*1/Pa[p_index]^2*(mut.strength/degree.plant)*sum(((z1[j]-m$mp[p_index])^2-Pa[p_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_g[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
    }
    
    N_strength[j] <-(sum(numer_a[j,])/(1+h*sum(denom_a[j,])))*w1[j]
    m_strength[j] <-(sum(numer_m[j,])/(1+h*sum(denom_m[j,])))*w1[j]
    g_strength[j] <-(sum(numer_g[j,])/(1+h*sum(denom_g[j,])))*w1[j]
  }
  
  
  G = sum(N_strength)
  B = sum(m_strength) 
  V = sum(g_strength)
  
  
  
  return(list(G= G, B = B,V=V))
  
}




type_2_animals<-function(m,sigma,w,h,np,na,mut.strength,a_index,
                         points,mat,degree.animal,Pa){
  temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
  z1<-gauss.quad.prob(points, dist = "normal", mu=m$ma[a_index], sigma =sqrt(sigma$sa[a_index]))$nodes #z of animals
  w1<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$ma[a_index],sigma =sqrt(sigma$sa[a_index]))$weights #p_i(z)
  
  z2<-matrix(0,nrow=length(np),ncol=points)
  w2<-matrix(0,nrow=length(np),ncol=points)
  numer_a<-denom_a<-numer_m<-denom_m<-numer_g<-denom_g<-matrix(0,nrow=points,ncol=length(np))
  N_strength<-m_strength<-g_strength<-numeric()
  
  
  for(j in 1:points){  
    for(k in 1:length(np)){
      z2[k,]<-gauss.quad.prob(points, dist = "normal", mu=m$mp[k], sigma =sqrt(sigma$sp[k]))$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w2[k,]<-gauss.quad.prob(points, dist = "normal", 
                              mu=m$mp[k],sigma =sqrt(sigma$sp[k]))$weights #pj(z'')
      
      numer_a[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_a[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_m[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum((z1[j]-m$ma[a_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_m[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_g[j,k]<- np[k]*mat[k,a_index]*1/Pa[a_index]^2*(mut.strength/degree.animal)*sum(((z1[j]-m$ma[a_index])^2-Pa[a_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_g[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
    }
    
    N_strength[j] <-(sum(numer_a[j,])/(1+h*sum(denom_a[j,])))*w1[j]
    m_strength[j] <-(sum(numer_m[j,])/(1+h*sum(denom_m[j,])))*w1[j]
    g_strength[j] <-(sum(numer_g[j,])/(1+h*sum(denom_g[j,])))*w1[j]
  }
  
  
  #  temp2[i]<-sum(np*(mut.strength/degree.animal)*exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]/(1+h*sum(denom[i,])))
  # f <-  exp(-(z1[i]- z2)^2/w^2) 
  #j1[i]<- sum(na*(mut.strength)*f/(1+h*np*(mut.strength/degree.animal)*f)^2*w2*w1[i])
  #gvar[i]<-1/P^2*sum(np*(mut.strength/degree.animal)*((z1[i]-m$ma)^2-P)*exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]/(1+h*exp(-(z1[i]- z2)^2/w^2)*np))
  #dat2[i]<- sum(np*(mut.strength/degree.animal)*(z1[i]-m$ma[index])*(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i])/(1+h*sum(denom[i,])))
  
  
  
  G = sum(N_strength)
  B = sum(m_strength) 
  V = sum(g_strength)
  
  
  
  return(list(G= G, B = B,V=V))
  
}

type_2_plants<-function(m,sigma,w,h,np,na,mut.strength,p_index,
                        points,mat,degree.plant,Pa){
  temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
  
  z1<-gauss.quad.prob(points, dist = "normal", mu=m$mp[p_index], sigma =sqrt(sigma$sp[p_index]))$nodes #z of animals
  w1<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$mp[p_index],sigma =sqrt(sigma$sp[p_index]))$weights #p_i(z)
  
  z2<-matrix(0,nrow=length(na),ncol=points)
  w2<-matrix(0,nrow=length(na),ncol=points)
  numer_a<-denom_a<-numer_m<-denom_m<-numer_g<-denom_g<-matrix(0,nrow=points,ncol=length(na))
  N_strength<-m_strength<-g_strength<-numeric()
  
  for(j in 1:points){  
    for(k in 1:length(na)){
      z2[k,]<-gauss.quad.prob(points, dist = "normal", mu=m$ma[k], sigma =sqrt(sigma$sa[k]))$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w2[k,]<-gauss.quad.prob(points, dist = "normal", 
                              mu=m$ma[k],sigma =sqrt(sigma$sa[k]))$weights #pj(z'')
      
      numer_a[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_a[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_m[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum((z1[j]-m$mp[p_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_m[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_g[j,k]<- na[k]*mat[p_index,k]*1/Pa[p_index]^2*(mut.strength/degree.plant)*sum(((z1[j]-m$mp[p_index])^2-Pa[p_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_g[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
    }
    
    N_strength[j] <-(sum(numer_a[j,])/(1+h*sum(denom_a[j,])))*w1[j]
    m_strength[j] <-(sum(numer_m[j,])/(1+h*sum(denom_m[j,])))*w1[j]
    g_strength[j] <-(sum(numer_g[j,])/(1+h*sum(denom_g[j,])))*w1[j]
  }
  
  
  G = sum(N_strength)
  B = sum(m_strength) 
  V = sum(g_strength)
  
  
  
  return(list(G= G, B = B,V=V))
  
}



type_2_animals<-function(m,sigma,w,h,np,na,mut.strength,a_index,
                         points,mat,degree.animal,Pa){
  temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
  z1<-gauss.quad.prob(points, dist = "normal", mu=m$ma[a_index], sigma =sqrt(sigma$sa[a_index]))$nodes #z of animals
  w1<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$ma[a_index],sigma =sqrt(sigma$sa[a_index]))$weights #p_i(z)
  
  z2<-matrix(0,nrow=length(np),ncol=points)
  w2<-matrix(0,nrow=length(np),ncol=points)
  numer_a<-denom_a<-numer_m<-denom_m<-numer_g<-denom_g<-matrix(0,nrow=points,ncol=length(np))
  N_strength<-m_strength<-g_strength<-numeric()
  
  
  for(j in 1:points){  
    for(k in 1:length(np)){
      z2[k,]<-gauss.quad.prob(points, dist = "normal", mu=m$mp[k], sigma =sqrt(sigma$sp[k]))$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w2[k,]<-gauss.quad.prob(points, dist = "normal", 
                              mu=m$mp[k],sigma =sqrt(sigma$sp[k]))$weights #pj(z'')
      
      numer_a[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_a[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_m[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum((z1[j]-m$ma[a_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_m[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_g[j,k]<- np[k]*mat[k,a_index]*1/Pa[a_index]^2*(mut.strength/degree.animal)*sum(((z1[j]-m$ma[a_index])^2-Pa[a_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_g[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
    }
    
    N_strength[j] <-(sum(numer_a[j,])/(1+h*sum(denom_a[j,])))*w1[j]
    m_strength[j] <-(sum(numer_m[j,])/(1+h*sum(denom_m[j,])))*w1[j]
    g_strength[j] <-(sum(numer_g[j,])/(1+h*sum(denom_g[j,])))*w1[j]
  }
  
  
  #  temp2[i]<-sum(np*(mut.strength/degree.animal)*exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]/(1+h*sum(denom[i,])))
  # f <-  exp(-(z1[i]- z2)^2/w^2) 
  #j1[i]<- sum(na*(mut.strength)*f/(1+h*np*(mut.strength/degree.animal)*f)^2*w2*w1[i])
  #gvar[i]<-1/P^2*sum(np*(mut.strength/degree.animal)*((z1[i]-m$ma)^2-P)*exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]/(1+h*exp(-(z1[i]- z2)^2/w^2)*np))
  #dat2[i]<- sum(np*(mut.strength/degree.animal)*(z1[i]-m$ma[index])*(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i])/(1+h*sum(denom[i,])))
  
  
  
  G = sum(N_strength)
  B = sum(m_strength) 
  V = sum(g_strength)
  
  
  
  return(list(G= G, B = B,V=V))
  
}

type_2_plants<-function(m,sigma,w,h,np,na,mut.strength,p_index,
                        points,mat,degree.plant,Pa){
  temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
  
  z1<-gauss.quad.prob(points, dist = "normal", mu=m$mp[p_index], sigma =sqrt(sigma$sp[p_index]))$nodes #z of animals
  w1<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$mp[p_index],sigma =sqrt(sigma$sp[p_index]))$weights #p_i(z)
  
  z2<-matrix(0,nrow=length(na),ncol=points)
  w2<-matrix(0,nrow=length(na),ncol=points)
  numer_a<-denom_a<-numer_m<-denom_m<-numer_g<-denom_g<-matrix(0,nrow=points,ncol=length(na))
  N_strength<-m_strength<-g_strength<-numeric()
  
  for(j in 1:points){  
    for(k in 1:length(na)){
      z2[k,]<-gauss.quad.prob(points, dist = "normal", mu=m$ma[k], sigma =sqrt(sigma$sa[k]))$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w2[k,]<-gauss.quad.prob(points, dist = "normal", 
                              mu=m$ma[k],sigma =sqrt(sigma$sa[k]))$weights #pj(z'')
      
      numer_a[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_a[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_m[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum((z1[j]-m$mp[p_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_m[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      numer_g[j,k]<- na[k]*mat[p_index,k]*1/Pa[p_index]^2*(mut.strength/degree.plant)*sum(((z1[j]-m$mp[p_index])^2-Pa[p_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_g[j,k]<- na[k]*mat[p_index,k]*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
    }
    
    N_strength[j] <-(sum(numer_a[j,])/(1+h*sum(denom_a[j,])))*w1[j]
    m_strength[j] <-(sum(numer_m[j,])/(1+h*sum(denom_m[j,])))*w1[j]
    g_strength[j] <-(sum(numer_g[j,])/(1+h*sum(denom_g[j,])))*w1[j]
  }
  
  
  G = sum(N_strength)
  B = sum(m_strength) 
  V = sum(g_strength)
  
  
  
  return(list(G= G, B = B,V=V))
  
}

fitness_f_P<-function(z,Np,Na, pars,Temp, muA,muP, varA, varP){
  
  b <- pars$gi/pars$bw*exp(-(Temp-z)^2/(2*pars$bw)) 
  
  m2.temp<-list(ma=muA,mp=muP)
  sigma2<-list(sa=varA,sp=varP)
  mut_temp<-azz<-numeric()
  for(i in 1:length(z)){
    azz[i] <-   sum(pars$w_c/sqrt(pars$w_c^2+ 2*varP)*exp(-(z[i]-muP)^2/(pars$w_c^2+2*varP))*Np) 
  }
  mut_temp<- type_2_plants_s1fig(m=m2.temp,sigma=sigma2,w=pars$w,h=0.25,np=Np,na=Na,
                                 mut.strength=pars$mut.strength,
                                 points=length(z),mat=pars$matrix,degree.plant=mean(pars$dgplants),
                                 Pa=varP,z = z)$G
  f_g <- b - mean(Np) - azz + mut_temp
  
  return(list(z=z,F_g=f_g))  
}

fitness_f_A<-function(z,Np,Na, pars,Temp,muA,muP, varA, varP){
  
  b <- pars$gi/pars$bw*exp(-(Temp-z)^2/(2*pars$bw)) 
  
  m2.temp<-list(ma=muA,mp=muP)
  sigma2<-list(sa=varA,sp=varP)
  mut_temp<-azz<-numeric()
  for(i in 1:length(z)){
    azz[i] <-   sum(pars$w_c/sqrt(pars$w_c^2+ 2*varA)*exp(-(z[i]-muA)^2/(pars$w_c^2+2*varA))*Na) 
  }
  mut_temp<- type_2_animals_s1fig(m=m2.temp,sigma=sigma2,w=pars$w,h=0.25,np=Np,na=Na,
                                  mut.strength=pars$mut.strength,
                                  points=length(z),mat=pars$matrix,degree.animal=mean(pars$dganimals),
                                  Pa=varA,z = z)$G
  
  
  
  
  f_g <- b - mean(Na) - azz + mut_temp
  
  return(list(z=z,F_g=f_g))  
}





type_2_animals_s1fig<-function(m,sigma,w,h,np,na,mut.strength,
                               points,mat,degree.animal,Pa,z){
  temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
  z1<- z 
  z2<-matrix(0,nrow=length(np),ncol=points)
  w2<-matrix(0,nrow=length(np),ncol=points)
  numer_a<-denom_a<-numer_m<-denom_m<-numer_g<-denom_g<-matrix(0,nrow=points,ncol=length(np))
  N_strength<-m_strength<-g_strength<-numeric()
  
  
  for(j in 1:points){  
    for(k in 1:length(np)){
      z2[k,]<-gauss.quad.prob(points, dist = "normal", mu=m$mp[k], sigma =sqrt(sigma$sp[k]))$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w2[k,]<-gauss.quad.prob(points, dist = "normal", 
                              mu=m$mp[k],sigma =sqrt(sigma$sp[k]))$weights #pj(z'')
      
      numer_a[j,k]<- np[k]*sum(mat[k,])*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_a[j,k]<- np[k]*sum(mat[k,])*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      # numer_m[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum((z1[j]-m$ma[a_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      #  denom_m[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      #   numer_g[j,k]<- np[k]*mat[k,a_index]*1/Pa[a_index]^2*(mut.strength/degree.animal)*sum(((z1[j]-m$ma[a_index])^2-Pa[a_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      #  denom_g[j,k]<- np[k]*mat[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
    }
    
    N_strength[j] <-(sum(numer_a[j,])/(1+h*sum(denom_a[j,])))
    # m_strength[j] <-(sum(numer_m[j,])/(1+h*sum(denom_m[j,])))
    #  g_strength[j] <-(sum(numer_g[j,])/(1+h*sum(denom_g[j,])))
  }
  return(list(G= N_strength))
  
}

type_2_plants_s1fig<-function(m,sigma,w,h,np,na,mut.strength,
                              points,mat,degree.plant,Pa,z){
  temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
  
  z1<- z #gauss.quad.prob(points, dist = "normal", mu=m$mp[p_index], sigma =sqrt(sigma$sp[p_index]))$nodes #z of animals
  #w1<- gauss.quad.prob(points, dist = "normal", 
  #            mu=m$mp[p_index],sigma =sqrt(sigma$sp[p_index]))$weights #p_i(z)
  
  z2<-matrix(0,nrow=length(na),ncol=points)
  w2<-matrix(0,nrow=length(na),ncol=points)
  numer_a<-denom_a<-numer_m<-denom_m<-numer_g<-denom_g<-matrix(0,nrow=points,ncol=length(na))
  N_strength<-m_strength<-g_strength<-numeric()
  
  for(j in 1:points){  
    for(k in 1:length(na)){
      z2[k,]<-gauss.quad.prob(points, dist = "normal", mu=m$ma[k], sigma =sqrt(sigma$sa[k]))$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w2[k,]<-gauss.quad.prob(points, dist = "normal", 
                              mu=m$ma[k],sigma =sqrt(sigma$sa[k]))$weights #pj(z'')
      
      numer_a[j,k]<- na[k]*sum(mat[,k])*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      denom_a[j,k]<- na[k]*sum(mat[,k])*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      #   numer_m[j,k]<- na[k]*sum(mat[,k])*(mut.strength/degree.plant)*sum((z1[j]-m$mp[p_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      #  denom_m[j,k]<- na[k]*sum(mat[,k])*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
      #     numer_g[j,k]<- na[k]*sum(mat[,k])*1/Pa[p_index]^2*(mut.strength/degree.plant)*sum(((z1[j]-m$mp[p_index])^2-Pa[p_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      #    denom_g[j,k]<- na[k]*sum(mat[,k])*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
      
    }
    
    N_strength[j] <-(sum(numer_a[j,])/(1+h*sum(denom_a[j,])))
  }
  
  
  #G = sum(N_strength)
  #B = sum(m_strength) 
  #V = sum(g_strength)
  
  
  
  return(list(G=N_strength))
  
}


eqs_type2_new<- function(time, state, pars) {
  A <- dim(pars$matrix)[2]  ## number of animal species
  P <-dim(pars$matrix)[1]
  Np<-state[(A+1):(A+P)]
  Na<-state[1:A]
  s <- pars$sig ## species' trait standard deviations
  ma<-state[(A+P+1):(A+P+A)]
  mp<-state[(A+P+A+1):(A+P+A+P)]
  gvarA<-state[(A+P+A+P+1):(A+P+A+P+A)]
  gvarP<-state[(A+P+A+P+A+1):(A+P+A+P+A+P)]
  
  s <- pars$sig ## species' trait standard deviations
  alpha.a<-pars$Amatrix ## alpha matrix
  alpha.p<-pars$Pmatrix
  
  aij<-bij<-vij<-matrix(0, nrow=A,ncol=P) 
  aji<-bji<-vji<-matrix(0, nrow=P,ncol=A) 
  muA<-ma
  muP<-mp
  aj<-bj<-ai<-bi<-vj<-vi<-numeric()
  #*cutoff(gvarA/(1e-7))
  varA <- ( (gvarA + pars$envA))
  varP <- ( (gvarP + pars$envP))
  
  noise_f_P<-approxfun(pars$noise_T$import[,1],method="linear", rule=2)
  Temp<- replicate(1,noise_f_P(time)) ##pars$Temp# 
  
  dmA <- outer(muA, muA, FUN="-") ## difference matrix of trait means
  dmP <- outer(muP, muP, FUN="-") ## difference matrix of trait means
  
  svA <- outer(varA, varA, FUN="+") ## sum matrix of trait variances
  svP <- outer(varP, varP, FUN="+") ## sum matrix of trait variances
  
  alphaA <- exp(-dmA^2/(2*svA+pars$w_c^2))*pars$w_c/sqrt(2*svA+pars$w_c^2) ## alpha matrix
  alphaP <-  exp(-dmP^2/(2*svP+pars$w_c^2))*pars$w_c/sqrt(2*svP+pars$w_c^2) ## alpha matrix
  
  diag(alphaA)<-0.5
  diag(alphaP)<-0.5
  
  betaA <- alphaA*2*varA*(-dmA)/(2*svA+pars$w_c^2)^1.5 ## beta matrix
  betaP <- alphaP*2*varP*(-dmP)/(2*svP+pars$w_c^2)^1.5 ## beta matrix
  diag(betaA)<-0
  diag(betaP)<-0
  
  xetaA  <- alphaA/(2*svA+pars$w_c^2)^2*(dmA^2-(2*svA+pars$w_c^2))
  xetaP <- alphaP/(2*svP+pars$w_c^2)^2*(dmP^2-(2*svP+pars$w_c^2))
  diag(xetaA)<-0
  diag(xetaP)<-0
  
  
  
  ba<- pars$gi/(pars$bw)*(pars$bw)/(sqrt((pars$bw)^2+varA))*exp(-(Temp- muA)^2/(2*(pars$bw)^2+varA)) - pars$ki
  bar_ba<- pars$gi/(pars$bw)*exp(-(Temp- muA)^2/(2*(pars$bw)^2+varA))*(varA)*(pars$bw)*(Temp -muA)/((pars$bw)^2 + varA)^1.5
  xar_ba<- (2*pars$gi*(2*muA^2-4*Temp*muA+2*Temp^2-2*pars$bw^2-varA)*exp(-(Temp-muA)^2/(2*pars$bw^2+varA)))/(sqrt(pars$bw^2+varA)*(2*pars$bw^2+varA)^2)
  
  bp<- pars$gi/(pars$bw)*(pars$bw)/(sqrt((pars$bw)^2+ varP))*exp(-(Temp- muP)^2/(2*(pars$bw)^2+varP)) - pars$ki
  bar_bp<-pars$gi/(pars$bw)*exp(-(Temp- muP)^2/(2*(pars$bw)^2+varP))*(varP)*(pars$bw)*(Temp - muP)/((pars$bw)^2 + varP)^1.5 
  xar_bp<- (2*pars$gi*(2*muP^2-4*Temp*muP+2*Temp^2-2*pars$bw^2-varP)*exp(-(Temp-muP)^2/(2*pars$bw^2+varP)))/(sqrt(pars$bw^2+varP)*(2*pars$bw^2+varP)^2)
  
  for(r in 1:A){
    
    
    m.temp<-list(ma=muA,mp=muP)
    sigma1<-list(sa=varA,sp=varP)
    temp1<-type_2_animals(m=m.temp,sigma=sigma1,w=pars$w,h=0.1,np=Np,na=Na,mut.strength=pars$mut.strength,a_index=r,
                          points=5,mat=pars$matrix,degree.animal=pars$dganimals[r],Pa=varA)
    
    
    ai[r]<-temp1$G
    bi[r]<-temp1$B
    vi[r]<-temp1$V
  }
  for(k in 1:P){
    
    m2.temp<-list(ma=muA,mp=muP)
    sigma2<-list(sa=varA,sp=varP)
    temp2<-type_2_plants(m=m2.temp,sigma=sigma2,w=pars$w,h=0.1,np=Np,na=Na,
                         mut.strength=pars$mut.strength,p_index=k,
                         points=5,mat=pars$matrix, 
                         degree.plant =pars$dgplants[k],
                         Pa=varP)
    
    aj[k]<-temp2$G
    bj[k]<-temp2$B
    vj[k]<-temp2$V
  }
  
  #pars$Amatrix
  #pars$Pmatrix
  dndt_a<- Na*(ba-alphaA%*%Na+ai)*cutoff(Na/(1e-6)) #  na*(ba-alpha.a%*%na+ai)*cutoff(na/(1e-8)) #population dynamics
  dndt_p<- Np*(bp-alphaP%*%Np+aj)*cutoff(Np/(1e-6))  #population dynamics
  dudt_A<- pars$trait_evolve*gvarA/varA*(bar_ba- betaA%*%Na+ bi) #mean trait dynamics
  dudt_P<- pars$trait_evolve*gvarP/varP*(bar_bp -betaP%*%Np+ bj) #mean trait dynamics
  dudt_vA<- pars$var_evolve*1/2*(gvarA/varA)^2*(xar_ba - xetaA%*%Na + vi)*cutoff(gvarA/(1e-7))
  dudt_vP<- pars$var_evolve*1/2*(gvarP/varP)^2*(xar_bp - xetaP%*%Np + vj)*cutoff(gvarP/(1e-7))
  
  ## return equations by first flattening them back into a single vector
  return(list(c(dndt_a, dndt_p,dudt_A,dudt_P,dudt_vA,dudt_vP)))
}


create_matrix<-function(SA,SP){
  web <- matrix(rbinom(SA*SP, 1, prob=1),nrow=SA,ncol =SP)
  return(web)
  
}



cluster_run_func<-function(params, ic, tmax ){
  
  OUT<-ode(func=eqs_type2_new, y=ic, parms=params,
           times=seq(0, tmax, by=50)) %>% 
    organize_results(pars = params) 
  
  
  #print(OUT %>% plot_density())
  
  
  (Na_r<-(OUT  %>% filter(time == tmax,type=="N"))$v)
  (Np_r<-(OUT %>% filter(time == tmax,type=="P"))$v)
  (mua_r<-(OUT  %>% filter(time == tmax,type=="muA"))$v)
  (mup_r<-(OUT %>% filter(time == tmax,type=="muP"))$v)
  (va_r<-(OUT  %>% filter(time == tmax,type=="sA"))$v)
  (vp_r<-(OUT %>% filter(time == tmax,type=="sP"))$v)
  
  Na_r[(Na_r < 1e-7)] <- 0
  Np_r[(Np_r < 1e-7)] <- 0
  
  A<-length(Na_r)
  P<-length(Np_r)
  na<-ic[1:length(Na_r)]
  np<-ic[(length(Na_r)+1) : (length(Na_r)+length(Np_r) )]
  
  N_initial<-c(na,np)   
  N_dat_min<-(OUT  %>% filter( time > 100 & time < tmax ,type=="N") %>% select(species, v) %>%  
                group_by(species) %>% summarise(min_n = min(v)))$min_n
  P_dat_min<-(OUT  %>% filter( time > 100 & time < tmax, type=="P") %>% select(species, v) %>%  
                group_by(species) %>% summarise(min_p = min(v)))$min_p
  max_Na_t <-(OUT  %>% filter( time == tmax, type=="N"))$v 
  max_Np_t <-(OUT  %>% filter( time == tmax, type=="P"))$v 
  
  N_dat_min[(N_dat_min< 1e-7)]<-0
  P_dat_min[(P_dat_min< 1e-7)]<-0
  
  max_Na_t[(max_Np_t< 1e-7)]<-0
  max_Np_t[(max_Np_t< 1e-7)]<-0
  
  
  N_min <-c(N_dat_min,P_dat_min)
  N_final <-c(max_Na_t,max_Np_t)
  
  index_N<- which(N_min < N_final & N_min <= N_initial)
  fraction_evo_recue<-length(index_N)/(A+P)
  N_evorescue<-N_initial
  N_evorescue[index_N] <- 1
  N_evorescue[which(N_evorescue != 1)] <- 0
  
  
  index_N_declines<- which(N_initial>N_min & N_min >= N_final & N_final < 1e-3)
  index_stasis <- which(round(N_initial,1) == round(N_min,1) & round(N_min,1) == round(N_final,1))
  index_growth<-which(N_initial < N_final & N_min >= N_initial & N_final>=N_min) 
  
  fraction_decline <- length(index_N_declines)/(A+P)
  fraction_stasis <- length(index_stasis)/(A+P)
  fraction_growth<- length(index_growth)/(A+P)
  
  
  N_decline <- N_initial
  N_decline[index_N_declines]<-1
  N_decline[which(N_decline!=1)]<-0
  
  
  N_stasis<-N_initial
  N_stasis[index_stasis]<-1
  N_stasis[which(N_stasis!=1)]<-0
  
  N_growth<-N_initial
  N_growth[index_growth]<-1
  N_growth[which(N_growth!=1)]<-0
  
  
  final_density<-c(Na_r,Np_r)
  g<-params$matrix
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] 
  proportion_evo_rescue<-fraction_evo_recue
  na<-ic[1:length(Na_r)]
  np<-ic[(length(Na_r)+1) : (length(Na_r)+length(Np_r) )]
  muA<-ic[(Aspecies+Plantspecies+1):(Aspecies+Plantspecies+Aspecies)]
  muP<-ic[(Aspecies+Plantspecies+Aspecies+1):(Aspecies+Plantspecies+Aspecies+Plantspecies)]
  Pa<-ic[(Aspecies+Plantspecies+Aspecies+Plantspecies+1):(Aspecies+Plantspecies+Aspecies+Plantspecies+Aspecies)]
  Pp<-ic[(Aspecies+Plantspecies+Aspecies+Plantspecies+Aspecies+1):(Aspecies+Plantspecies+Aspecies+Plantspecies+Aspecies+Plantspecies)]
  
  
  #quantifying competition strength:
  dmA <- outer(muA, muA, FUN="-") ## difference matrix of trait means
  dmP <- outer(muP, muP, FUN="-") ## difference matrix of trait means
  
  svA <- outer(Pa, Pa, FUN="+") ## sum matrix of trait variances
  svP <- outer(Pp, Pp, FUN="+") ## sum matrix of trait variances
  
  
  comp_str_A<- exp(-dmA^2/(2*svA+params$w_c^2))*params$w_c/sqrt(2*svA+params$w_c^2)
  comp_str_P<- exp(-dmP^2/(2*svP+params$w_c^2))*params$w_c/sqrt(2*svP+params$w_c^2)
  cumulative_comp_A<-comp_str_A%*%na
  cumulative_comp_P<-comp_str_P%*%np
  cumulative_comp<-c(cumulative_comp_A,cumulative_comp_P)
  
  #trait lag
  trait_lag_A <- params$Temp - muA  
  trait_lag_P <- params$Temp - muP
  trait_lag<-c(trait_lag_A,trait_lag_P)
  
  #Absolute trait value at equilibrium
  eq_muA<-muA
  eq_muP<-muP
  eq_mu<-c(muA,muP)
  
  
  #if(index_top_g1 > dim(g)[1]){
  # g[,(index_top_g1-Aspecies)]
  #}
  
  
  
  network_metric <- networklevel(params$matrix)
  graph_mat<-graph_from_biadjacency_matrix(t(g),weighted = T)
  
  extinctions <- length(which( final_density < 1e-3))/(Aspecies+Plantspecies)
  change_in_density <- sum(final_density)-sum(c(na,np))
  richness <- length(which(final_density>1e-3) )/(Aspecies+Plantspecies)
  change_richness <- length(which(final_density>1e-3) ) -length(which(c(na,np)>1e-3) ) 
  community_biomass<- sum(final_density)
  change_community_biomass<- sum(final_density) - sum(c(na,np))
  NODF <- nestedness_NODF(g)
  connectance <- Connectance(g)
  network_size <- (Aspecies+Plantspecies)
  modularity<-network_metric[[7]]
  avg_betweeness <- mean(igraph::betweenness(graph_mat))
  sd_betweeness <- sd(igraph::betweenness(graph_mat))
  
  sa_sp<-OUT %>% filter(type %in%c("sA","sP"))
  dd<-sa_sp %>%
    group_by(type, species) %>%
    summarize(max_v = max(v, na.rm = TRUE))
  
  result <- OUT %>% filter(type %in%c("muA","muP")) %>%     # Filter only relevant types
    mutate(diff = abs(v - 28)) %>%                # Calculate absolute difference from 28
    filter(diff <= 2) %>%                         # Keep rows with a difference of 1 or less
    group_by(type,species) %>%                            # Group by `type`
    slice_min(time, with_ties = FALSE) %>%        # Get the minimum time for each type
    select(type, species, time, v)   
  
  degree.animals<-params$degree[1:Aspecies]
  degree.plants<-params$degree[(Aspecies+1):(Aspecies+Plantspecies)]
  degree_a<- degree.animals[(result %>% filter(type =="muA"))$species]
  degree_p<-degree.plants[(result %>% filter(type =="muP"))$species]
  
  result$degree<-c(degree_a,degree_p)
  
  species_index<-params$species_index
  final_density_1<-final_density
  final_density_1[which(final_density_1 < 1e-4)]<-0
  final_density_1[which(final_density_1 > 1e-4)]<-1
  extinct_binomial<-final_density_1
  proportion_survived<-sum(final_density_1)/(Aspecies+Plantspecies)
  rownames(g) <- as.character(seq(1:nrow(g)))
  Ls<-c(letters,LETTERS, seq(500,800,1))
  colnames(g) <- Ls[c(1:ncol(g))]
  index_names_of_g<-c(Ls[c(1:ncol(g))],as.character(seq(1:nrow(g))))   
  target_generalist<-index_names_of_g[species_index]
  min_no_path_to_generalist<-numeric()
  for(i in 1:(Aspecies+Plantspecies)){
    
    min_no_path_to_generalist[i]<- min_no_path_to_target_node(graph_matrix = g,start_node = index_names_of_g[i],
                                                              target_node = target_generalist)
    
    
  }
  # making sure the index of the generalist path to itself is zero
  min_no_path_to_generalist[species_index] <- 0
  #min_no_path_to_generalist[is.na(min_no_path_to_generalist)]<- Aspecies+Plantspecies
  
  min_no_path_to_generalist[species_index] <- 0
  min_path_a<- min_no_path_to_generalist[1:Aspecies]
  min_path_p<- min_no_path_to_generalist[(Aspecies+1):(Aspecies+Plantspecies)]
  
  result$min_path_generalist<-c(min_path_a[result$species[result$type=="muA"]],min_path_p[result$species[result$type=="muP"]])
  
  percentage_change_biomass<- ((c(Na_r,Np_r) - c(na,np))/c(na,np))*100
  
  species_level_data <-as.data.frame(cbind(
    rep(seq(1,(nrow(g)+ncol(g)),1)), #number of species
    (c(mua_r, mup_r) - c(muA,muP)),
    (c(va_r,vp_r) - c(Pa,Pp)),
    c(params$envA,params$envP),
    as.numeric(c(Na_r,Np_r )),
    as.numeric(proportion_survived),
    as.numeric(fraction_evo_recue),
    as.numeric(fraction_decline),
    as.numeric(fraction_growth),
    as.numeric(fraction_stasis),
    dd$max_v,
    cumulative_comp,
    trait_lag,
    N_decline,
    N_stasis,
    N_growth,
    N_evorescue,
    params$mut.strength[1],
    params$ki,
    params$Temp,
    params$w,
    params$w_c,
    percentage_change_biomass,
    igraph::betweenness(graph_mat),
    as.numeric(min_no_path_to_generalist),
    rep(as.character(params$web.name),each=(Aspecies+Plantspecies)),
    params$degree,
    as.character(params$variation),
    as.numeric(rep(nestedness_NODF(params$matrix), each=((Aspecies+Plantspecies)) )),
    as.numeric(rep(Connectance(params$matrix), each=((Aspecies+Plantspecies)) )),
    as.numeric(rep( (Aspecies+Plantspecies),each=((Aspecies+Plantspecies)))),
    as.character(params$h2),
    rep(modularity,each=(Aspecies+Plantspecies))))
  
  
  
  colnames(species_level_data)<-c("Species",
                                  "change_in_mean_trait",
                                  "change_in_mean_gvariance",
                                  "environmental_variance",
                                  "density", 
                                  "proportion_persisted",
                                  "proportion_rescued",
                                  "proportion_declined",
                                  "proportion_growth",
                                  "proportion_stasis",
                                  "max_genvar",
                                  "competition_strength",
                                  "trait_lag",
                                  "decline",
                                  "stasis",
                                  "growth",
                                  "evorescue",
                                  "mutualism_strength", 
                                  "mortality",
                                  "Temperature_shift",
                                  "width_mutualism",
                                  "width_competition",
                                  "percent_change_biomass",
                                  "betweeness_centrality",
                                  "minimum_path_generalist",
                                  "webname",
                                  "Degree",
                                  "variation",
                                  "Nestedness", 
                                  "Connectance",
                                  "Network_size",
                                  "h2",
                                  "modularity")
  
  # species_data<-rbind(species_data, species_level_data)
  
  final_output<- list(output=result, species_level_data=species_level_data)
  return(final_output)
  
}
# function for evaluating tipping points for each species

species.tipping.points<-function(){
  no.of.plants<-length(dat$Plants[1,])
  no.of.animals<-length(dat$Animals[1,])
  mean.animal.abundance<-mean.plant.abundance<-numeric()
  for( i in 1:no.of.plants){
    mean.plant.abundance[i]<-mean(dat$Plants[300:1000,i],na.rm=T)
  }
  for(j in 1:no.of.animals){
    mean.animal.abundance[j]<-mean(dat$Animals[300:1000,i],na.rm=T)
  }
  
  
}



## Organize simulation results into tidy table
## Input:
## - sol: output produced by the function ode()
## - pars: list of parameters, with the following elements:
## Output:
## - a tibble with columns: time, species, n (density), m (trait mean),
##   sigma
organize_results <- function(sol, pars) {
  S <- length(pars$sig) ## number of species
  A<-dim(pars$matrix)[2] # no. of animals
  P<-dim(pars$matrix)[1] # no. of plants
  temp<- sol %>% as.data.frame %>% as_tibble ## convert to tibble
  ## name the first column "time"
  temp<- temp 
  names(temp)[1] <- "time"
  names(temp)[2:(A+1)] <- paste0("N_", 1:(A)) ## name abundance columns (n_k)
  
  names(temp)[(A+2):(A+1+P)] <- paste0("P_", 1:P) ## name trait mean columns
  names(temp)[(A+P+2):(A+P+A+1)] <-  paste0("muA_", 1:(A))
  names(temp)[(A+P+A+2):(A+P+A+P+1)] <- paste0("muP_", 1:(P))
  names(temp)[(A+P+A+P+2):(A+P+A+P+A+1)] <- paste0("sA_", 1:(A))
  names(temp)[(A+P+A+P+A+2):(A+P+A+P+A+P+1)]<-paste0("sP_", 1:(P))
  
  temp <- temp %>%
    tidyr::gather("variable", "v", 2:ncol(temp)) %>% ## normalize the data
    tidyr::separate(variable, c("type", "species"), sep="_") %>%
    #spread(type, v) %>% ## separate columns for animal densities n and plant densities m
    dplyr::select(time, type, species,v) %>% ## rearrange columns
    mutate(species=as.integer(species), w=pars$gamma,
           Nestedness=pars$nestedness, Connectance=pars$C,
           theta=pars$theta,Web.name=pars$web.name) ## add params
  return(as_tibble(temp))
}


## Plot time series of densities, time series of trait values, and
## snapshot of the trait distributions at time = moment
## Input:
## - dat: data generated by organize_results()
## - moment: time at which trait distribution should be plotted
## - limits: a vector of two entries (x_low, x_high) for the x-axis limits
## - res: number of evenly spaced sampling points along the trait axis
##               for the trait distribution plot
## Output:
## - a ggplot2 plot with three panels in one column: abundance time series,
##   trait value time seties, and snapshot of trait distribution
plot_all <- function(dat, moment=0, limits=c(-1, 1), res=1001) {
  plot_grid(plot_density(dat), ncol=1, align="hv") %>%
    return
}



## Plot species densities through time
## Input:
## - dat: data generated by organize_results()
## Output:
## - a ggplot2 plot
## used to produce figure 1.
plot_density<- function(dat) {
  dat %>%
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=1.25) +
    scale_y_continuous(name="population density") +
    theme(legend.position="none") + facet_wrap(.~type,scales = "free") %>%
    return
}

## Plot species densities through time
## Input:
## - dat: data generated by organize_results()
## Output:
## - a ggplot2 plot
## used to produce figure 1.
plot_density_timeseries<- function(dat) {
  dat %>% filter(type == c("N","P")) %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=1.25) +
    scale_y_continuous(name="population density") +
    theme_cowplot()+
    theme(legend.position="none") + 
    facet_wrap(.~type,scales = "free") %>%
    return
}

plot_density_timeseries_genvar<- function(dat) {
  dat %>% filter(type == c("sA","sP")) %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=1.25) +
    scale_y_continuous(name="Genetic variance") +
    theme_cowplot()+
    theme(legend.position="none") + 
    facet_wrap(.~type,scales = "free") %>%
    return
}

#lay out function for multiple plots
lay_out = function(...) {    
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                                           layout.pos.col = x[[i]][[3]]))
  }
} 



#multiplot of ggplot2 figures with a common shared legend. Code taken from :https://rpubs.com/sjackman/grid_arrange_shared_legend
grid_arrange_shared_legend <- function(..., ncol, nrow, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol =ncol , nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}


# plots the density distribution  at a particular timepoint. This function was used to produce figure 1.
# Na: Abundance of animals at equilibrium
# Np: Abundance of plants at equilibrium
# m: mean traits at equilibrium
# sigma: variance of traits
# moment: mean
# limits: limits of the mean trait axis which in the study are -1,1
# code adapted from Barabas & D'Andrea 2016 Eco. Letts.

plot_snapshot <- function(Na, Np, m, sigma, Temp, limits=c(-1, 1), res=1001) {
  S_a <- length(Na) ## number of species
  S_p <- length(Np)
  ma<- m[1:(S_a)]
  mp<- m[(S_a+1):(S_a+S_p)]
  sigma_a <-sigma[1:(S_p)]
  sigma_p <- sigma[(S_a+1):(S_a+S_p)]
  traitaxis <- seq(limits[1], limits[2], l=res) ## sampling the trait axis
  #snap <- dat %>% filter(time==moment) %>% select(-time) ## time = moment
  traits_a <- expand.grid(species=1:S_a, trait=traitaxis) %>% as_tibble ## trait table
  traits_p <- expand.grid(species=(S_a+1):(S_a+S_p), trait=traitaxis) %>% as_tibble ## trait table
  
  traits_a["density"] <- 0 ## add column for population densities
  traits_p["density"] <- 0
  
  for (i in 1:S_a) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_a$density[(traits_a$species==i)] <- Na[i]*
      dnorm(traits_a$trait[(traits_a$species==i)], ma[i], sigma_a[i]) ## times density
  }
  traits_a$density[traits_a$density<max(traits_a$density)/1e3] <- NA
  
  for (i in 1:S_p) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_p$density[(traits_p$species==((S_a+1)+i))] <- Np[i]*dnorm(traits_p$trait[(traits_p$species==((S_a+1)+i))], 
                                                                     mp[i], sigma_p[i]) ## times density
  }
  traits_p$density[traits_p$density<max(traits_p$density)/1e3] <- NA
  
  landscape <- tibble(trait = traitaxis) %>% # for plotting intrinsic rates
    mutate(r= (1.15/2*exp(-(Temp-trait)^2/(2*2^2)))-0.2) %>%      
    mutate(r=ifelse(r<=0, NA, r)) %>%
    mutate(r=r*max(c(traits_a$density,traits_p$density),na.rm = T))
  
  
  traits<-data.frame(rbind(traits_a,traits_p), 
                     species_group=c(rep("Animals", nrow(traits_a)),
                                     rep("Plants", nrow(traits_p))))
  
  ggplot(traits) + ## generate plot
    geom_line(aes(x=trait, y=density, colour=factor(species)), na.rm=TRUE) +
    geom_ribbon(aes(x=trait, ymin=0, ymax=density, fill=factor(species)),
                alpha=0.15, colour=NA)+scale_fill_viridis_d(alpha = 1)+
    facet_wrap(.~species_group, nrow = 2,scales = "free")+
    theme_classic()+
    theme(legend.title = element_text(size = 14), 
          legend.position = "right", panel.background = element_blank(), 
          axis.text = element_text(colour = "black", size = 14), 
          axis.title = element_text(size = 14), 
          legend.text = element_text(size = 14), legend.key = element_blank(),
          strip.text.x = element_text(size= 14))+
    geom_line(data=landscape, aes(x=trait, y=r), linetype="dashed",
              colour="black", alpha=0.75, na.rm=TRUE) +
    scale_x_continuous(name="trait value", limits=limits) +
    scale_y_continuous(name="density", limits=c(0, NA)) +
    scale_color_viridis_d(alpha=1)+
    theme(legend.position="none") %>%
    return 
}


plot_snapshot_S1 <- function(Na, Np, m, sigma, Temp, pars,limits=c(-1, 1), res=301) {
  S_a <- length(Na) ## number of species
  S_p <- length(Np)
  ma<- m[1:(S_a)]
  mp<- m[(S_a+1):(S_a+S_p)]
  sigma_a <-sigma[1:(S_p)]
  sigma_p <- sigma[(S_a+1):(S_a+S_p)]
  traitaxis <- seq(limits[1], limits[2], l=res) ## sampling the trait axis
  #snap <- dat %>% filter(time==moment) %>% select(-time) ## time = moment
  traits_a <- expand.grid(species=1:S_a, trait=traitaxis) %>% as_tibble ## trait table
  traits_p <- expand.grid(species=(S_a+1):(S_a+S_p), trait=traitaxis) %>% as_tibble ## trait table
  
  traits_a["density"] <- 0 ## add column for population densities
  traits_p["density"] <- 0
  Na<-Na/(sum(Na))
  Np<-Np/sum(Np)
  for (i in 1:S_a) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_a$density[(traits_a$species==i)] <- Na[i]*
      dnorm(traits_a$trait[(traits_a$species==i)], ma[i], sigma_a[i]) ## times density
  }
  traits_a$density[traits_a$density<max(traits_a$density)/1e3] <- NA
  
  for (i in 1:S_p) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_p$density[(traits_p$species==((S_a+1)+i))] <- Np[i]*dnorm(traits_p$trait[(traits_p$species==((S_a+1)+i))], 
                                                                     mp[i], sigma_p[i]) ## times density
  }
  traits_p$density[traits_p$density<max(traits_p$density)/1e3] <- NA
  
  f_P<-fitness_f_P(z = traitaxis,Np = Np,Na = Na,pars = pars,Temp=Temp,muA = ma,muP = mp,
                   varA =sigma_a,varP=sigma_p)$F_g
  f_A<-fitness_f_A(z = traitaxis,Np = Np,Na = Na,pars = pars,Temp=Temp,muA = ma,muP = mp,
                   varA =sigma_a,varP = sigma_p)$F_g
  
  landscape_a <- tibble(trait = traitaxis) %>% # for plotting intrinsic rates
    mutate(r= f_A)
  
  landscape_p <- tibble(trait = traitaxis) %>% # for plotting intrinsic rates
    mutate(r= f_P)
  
  landscape<-data.frame(rbind(landscape_a,landscape_p), 
                        species_group=c(rep("Animals", nrow(landscape_a)),
                                        rep("Plants", nrow(landscape_p))))
  
  traits<-data.frame(rbind(traits_a,traits_p), 
                     species_group=c(rep("Animals", nrow(traits_a)),
                                     rep("Plants", nrow(traits_p))))
  
  ggplot(traits) + ## generate plot
    geom_line(aes(x=trait, y=density, colour=factor(species)), na.rm=TRUE) +
    geom_ribbon(aes(x=trait, ymin=0, ymax=density, fill=factor(species)),
                alpha=0.15, colour=NA)+scale_fill_viridis_d(alpha = 1)+
    facet_wrap(.~species_group, nrow = 2,scales = "free")+
    theme_classic()+
    theme(legend.title = element_text(size = 14), 
          legend.position = "right", panel.background = element_blank(), 
          axis.text = element_text(colour = "black", size = 14), 
          axis.title = element_text(size = 14), 
          legend.text = element_text(size = 14), legend.key = element_blank(),
          strip.text.x = element_text(size= 14))+
    geom_line(data=landscape, aes(x=trait, y=r),
              colour="firebrick", alpha=1, na.rm=TRUE) +
    geom_hline(yintercept = 0,linetype="dashed", color="black")+
    scale_x_continuous(name="trait value", limits=limits) +
    scale_y_continuous(sec.axis = sec_axis(~.*1, name="Fitness")) +
    scale_color_viridis_d(alpha=1)+
    theme(legend.position="none") %>%
    return 
}



plot_snapshot_per_cap_growth <- function(Na, Np, m, sigma, Temp, pars,limits=c(-1, 1), res=301) {
  S_a <- length(Na) ## number of species
  S_p <- length(Np)
  ma<- m[1:(S_a)]
  mp<- m[(S_a+1):(S_a+S_p)]
  sigma_a <-sigma[1:(S_p)]
  sigma_p <- sigma[(S_a+1):(S_a+S_p)]
  traitaxis <- seq(limits[1], limits[2], l=res) ## sampling the trait axis
  #snap <- dat %>% filter(time==moment) %>% select(-time) ## time = moment
  traits_a <- expand.grid(species=1:S_a, trait=traitaxis) %>% as_tibble ## trait table
  traits_p <- expand.grid(species=(S_a+1):(S_a+S_p), trait=traitaxis) %>% as_tibble ## trait table
  
  traits_a["density"] <- 0 ## add column for population densities
  traits_p["density"] <- 0
  Na<-Na/(sum(Na))
  Np<-Np/sum(Np)
  for (i in 1:S_a) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_a$density[(traits_a$species==i)] <- Na[i]*
      dnorm(traits_a$trait[(traits_a$species==i)], ma[i], sigma_a[i]) ## times density
  }
  traits_a$density[traits_a$density<max(traits_a$density)/1e3] <- NA
  
  for (i in 1:S_p) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_p$density[(traits_p$species==((S_a+1)+i))] <- Np[i]*dnorm(traits_p$trait[(traits_p$species==((S_a+1)+i))], 
                                                                     mp[i], sigma_p[i]) ## times density
  }
  traits_p$density[traits_p$density<max(traits_p$density)/1e3] <- NA
  
  f_P<-fitness_f_P(z = traitaxis,Np = Np,Na = Na,pars = pars,Temp=Temp,muA = ma,muP = mp,
                   varA =sigma_a,varP=sigma_p)$F_g
  f_A<-fitness_f_A(z = traitaxis,Np = Np,Na = Na,pars = pars,Temp=Temp,muA = ma,muP = mp,
                   varA =sigma_a,varP = sigma_p)$F_g
  
  landscape_a <- tibble(trait = traitaxis) %>% # for plotting intrinsic rates
    mutate(r= f_A)
  
  landscape_p <- tibble(trait = traitaxis) %>% # for plotting intrinsic rates
    mutate(r= f_P)
  
  landscape<-data.frame(rbind(landscape_a,landscape_p), 
                        species_group=c(rep("Animals", nrow(landscape_a)),
                                        rep("Plants", nrow(landscape_p))))
  
  traits<-data.frame(rbind(traits_a,traits_p), 
                     species_group=c(rep("Animals", nrow(traits_a)),
                                     rep("Plants", nrow(traits_p))))
  
  ggplot(landscape) + ## generate plot
    #geom_line(aes(x=trait, y=density, colour=factor(species)), na.rm=TRUE) +
    # geom_ribbon(aes(x=trait, ymin=0, ymax=density, fill=factor(species)),
    #            alpha=0.15, colour=NA)
    geom_line(data=landscape, aes(x=trait, y=r),
              colour="firebrick", alpha=1, na.rm=TRUE) +
    scale_fill_viridis_d(alpha = 1)+
    facet_wrap(.~species_group, nrow = 2,scales = "free")+
    theme_classic()+
    theme(legend.title = element_text(size = 14), 
          legend.position = "right", panel.background = element_blank(), 
          axis.text = element_text(colour = "black", size = 14), 
          axis.title = element_text(size = 14), 
          legend.text = element_text(size = 14), legend.key = element_blank(),
          strip.text.x = element_text(size= 14))+
    geom_hline(yintercept = 0,linetype="dashed", color="black")+
    scale_x_continuous(name="trait value", limits=limits) +
    scale_y_continuous( name="Per-capita growth") +
    scale_color_viridis_d(alpha=1)+
    theme(legend.position="none") %>%
    return 
}


# plot species' density distributions in trait space (in 1D)
# Input:
# - snap: data generated by "organize_results()" in "solve_eqs.R", restricted
#         to a single point in time (i.e., a snapshot of the dynamics)
# - limits: vector with three entries: lower limit of trait axis, upper limit
#           of trait axis, and upper limit of y-axis. If this last entry is set
#           to NA, the y-axis will be scaled automatically
# - res: resolution of the trait axis, the number of grid points it gets
#        subdivided into to generate the plot
# - b: a function, giving the intrinsic growth rate as a function of phenotype z
# Output:
# - a plot of the species' density distributions


#computes the raw NODF taken from Song et al 2017 J. Animal Ecology
#input: web = mutualistic network
#output: raw NODF of the given network
nestedness_NODF <- function(web){
  web[web > 0] = 1
  SA <- nrow(web)
  SP <- ncol(web)
  N <- t(web) %*% web
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SP,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SP)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele == 0] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)] = 1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n1 <- sum(nes)
  
  N <- web %*% t(web)
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SA,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SA)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele ==0 ] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)]=1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n2 <- sum(nes)
  out <- 2*(n1 + n2) / (SA*(SA-1)+SP*(SP-1))
  return(out)
}




# function for competition coefficients within a guild.
# competitive interactions  were scaled by the total number of species within a guild as Dakos & Bascompte 2014 PNAS.
# matrix: network of interactions which are 0 or 1. 
# strength: average competition strength


mat.comp<-function(matrix,degree.animals,degree.plants,max_strength){
  Aspecies<- dim(matrix)[2]
  Plantspecies<- dim(matrix)[1]
  
  Amatrix<-matrix(runif(Aspecies^2, max_strength, 2*max_strength), nrow=Aspecies, ncol = Aspecies)
  diag(Amatrix)<-0.5
  #diag(Amatrix)<-  diag(Amatrix) #/degree.animals
  
  Pmatrix<-matrix(runif(Plantspecies^2, max_strength, 2*max_strength), nrow=Plantspecies, ncol = Plantspecies)
  diag(Pmatrix)<-0.5
  #diag(Pmatrix)<-2-  diag(Pmatrix)/degree.plants
  out<-return(list(Amatrix=Amatrix,Pmatrix=Pmatrix))
  
}

# measures connectance of a web network
Connectance<-function(web)
{
  return(sum(web)/(ncol(web)*nrow(web)))}



trait.matching<-function(mA,mP,adj.mat_1,gamma){
  tm<-numeric()
  for(i in 1:nrow(adj.mat_1)){
    tm[i] <- mean(adj.mat_1[i,]*exp(-(mA-mP[i])^2)/gamma)
    
  }
  return(tm=mean(tm))
}


community.w.mean<-function(mA,mP,Na,Np){
  
  cwm_plants <-sum(mA*Na/sum(Na))
  cwm_animals <-sum(mP*Np/sum(Np))
  m<-c(mA,mP)
  N<-c(Na,Np)
  cwm_community <-sum(m*N/sum(N))
  
  return(list(cwm_plants=cwm_plants,
              cwm_animals=cwm_animals,
              cwm_community=cwm_community))
}


#adjmat = adjacency matrix of which species interacts with whom. Ones and zeros.
#Na,Np = vector of animal and plant species abundance at evolutionary equilibrium
#Aspecies, Pspecies = no. of animal and plant species
#comp_mat_p = competition matrix of plants
#comp_mat_a = competition matrix for animals
# Ja, Jp are the matrix that includes the mutualistic interaction after the derivative
# halpha is the matrix of higherorder interactions which is of species x species dimension
jacobian_mat<-function(Na,Np,Aspecies,Pspecies, comp_mat_a,comp_mat_p, Ja, Jp){
  #the goal is to create a species x species jacobian matrix
  
  #creating the jacobian for the competition matrix : -aij*Ni, Ni is the equilibrium abundances of all species
  
  J<-matrix(0,nrow=(Aspecies+Pspecies),ncol=(Aspecies+Pspecies))
  J_m<-array(NA,dim=c( (Aspecies+Pspecies), (Aspecies+Pspecies)))
  Aij<-array(NA,dim=c( (Aspecies+Pspecies), (Aspecies+Pspecies)))
  
  Aij<-cbind(rbind(comp_mat_a, matrix(0,nrow=Pspecies,ncol=Aspecies)),
             rbind(matrix(0,nrow=Aspecies,ncol=Pspecies), comp_mat_p))
  
  J_m<-rbind(cbind(matrix(0,nrow=Aspecies,ncol=Aspecies), Ja),
             cbind(Jp,matrix(0,nrow=Pspecies,ncol=Pspecies)))
  
  
  
  N<-rbind(Na,Np)
  
  #creating the jacobian for the mutualism part of the model
  
  
  #creating the jacobian for the competition part of the model
  Jc<-Jm<-array(NA,dim=c((Aspecies+Pspecies),(Aspecies+Pspecies)))
  # Jm<-J_hoi<-matrix(NA,nrow=(Aspecies+Pspecies), ncol = (Aspecies+Pspecies))
  
  for(i in 1:(Aspecies+Pspecies)){
    for(j in 1:(Aspecies+Pspecies)){
      
      Jc[i,j] <- -Aij[i,j]*N[i]
      Jm[i,j] <- J_m[i,j]*N[i]
      
    }
  }
  
  
  
  #Jacobian matrix derived from the combination of competition, mutualism and hoi interaction
  #
  RH<-Avg.Robust<-dominant.eg<-numeric()
  
  
  Jij<- Jc + J_m
  
  
  eJ  <- eigen(Jij, only.values=TRUE)$values # eigenvalues
  dominant.eg<-eJ[1]
  Avg.Robust <- exp(mean(log(abs(eJ))))
  RH <- exp(sd(log(abs(eJ))))  
  
  
  
  return(list(Jac=Jij, eJ=eJ,Avg.Robust=Avg.Robust,RH=RH, dominant.eg=dominant.eg))
  
  
}


