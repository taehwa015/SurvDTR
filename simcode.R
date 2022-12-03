######################################################
# Title: Accountable survival contrast-learning for  # 
#        optimal dynamic treatment regime            #
# Maintainer: Taehwa Choi                            #
######################################################
data_generator_dtr2=function(n,p,c0,logit=c("t","f","r"))
{
  x=as.data.frame(matrix(runif(n*p,-2,2),n,p))
  colnames(x)=paste0("x",1,1:p)
  x2=with(x,runif(n,min(x11),max(x11)))
  
  if(logit=="t") pii=with(x,expit(x12-0.6*x13))
  if(logit=="f") pii=with(x,expit(x12-0.6*x13-0.4*x13^2))
  if(logit=="r") pii=0.5
  A1=with(x,rbinom(n,1,pii))
  if(logit=="t") pii2=with(x,expit(x2/2))
  if(logit=="f") pii2=with(x,expit(x2/2-0.3*x2^2))
  if(logit=="r") pii2=0.5
  A2=rbinom(n,1,pii2)
  
  T1=with(x,1.5+0.5*x11+A1*(x12-0.5)+log(rexp(n)))
  T2=with(x,1.5+0.5*x11+A1*(x12-0.5)+A2*(x2-0.5)+log(rexp(n)))
  T2opt=((x2-0.5>0)-A2)*(x2-0.5)
  if(cdist=="exp") C=rexp(n,c0)
  if(cdist=="con") C=rexp(n,c0*with(x,exp(x11-x12)))
  if(cdist=="unf") C=runif(n,10,200)
  eta=(T1<C)*1
  T=case_when(eta==1~exp(T2),eta==0~exp(T1+T2opt))
  Y=pmin(T,C)
  delta=ifelse(T<C,1,0)
  
  data.frame(Y,delta,A1,A2,x,x2,eta)
}


discriminator_dtr2=function(A1,A2)
{
  T1=with(xt,1.5+0.5*x11+A1*(x12-0.5)+log(rexp(N)))
  T2=with(xt,1.5+0.5*x11+A1*(x12-0.5)+A2*(x2t-0.5)+log(rexp(N)))
  T2opt=((x2t-0.5>0)-A2)*(x2t-0.5)
  if(cdist=="exp") C=rexp(N,c0)
  if(cdist=="con") C=rexp(N,c0*with(xt,exp(x11-x12)))
  eta=(T1<C)*1
  T=case_when(eta==1~exp(T2),eta==0~exp(T1+T2opt))
  Y=pmin(T,C)
  delta=ifelse(T<C,1,0)
  summary(survfit(Surv(Y,delta)~1),times=t0)$surv
}


true_generator_dtr2=function(n=N,p=3)
{
  x=as.data.frame(matrix(runif(n*p,-2,2),n,p))
  colnames(x)=paste0("x",1,1:p)
  x2=with(x,runif(n,min(x11),max(x11)))
  A1=with(x,1*(0.9*x12-0.436>=0))
  A2=1*(x2-0.5>=0)
  
  T1=with(x,1.5+0.5*x11+A1*(x12-0.5)+log(rexp(n)))
  T2=with(x,1.5+0.5*x11+A1*(x12-0.5)+A2*(x2-0.5)+log(rexp(n)))
  T2opt=((x2-0.5>0)-A2)*(x2-0.5)
  if(cdist=="exp") C=rexp(n,c0)
  if(cdist=="unf") C=runif(n,10,200)
  if(cdist=="con") C=rexp(n,c0*with(x,exp(x11-x12)))
  eta=(T1<C)*1
  T=case_when(eta==1~exp(T2),eta==0~exp(T1+T2opt))
  Y=pmin(T,C)
  delta=ifelse(T<C,1,0)
  summary(survfit(Surv(Y,delta)~1),times=t0)$surv
}

pkg=c("survival","pseudo","MASS","tictoc",'dplyr','lpSolve',
      "e1071","cmprsk","glmnet","progress")
sapply(pkg,library,character.only=T)
expit=function(x)exp(x)/(1+exp(x))
n=500;p=10
cdist=c("exp","con")[1]
t0=3
c0=exp(-3.6)
logit="f"

#Testdata
set.seed(1)
N=50000
testdat=data_generator_dtr2(N,p,c0,logit)
mean(testdat$delta)
mean(testdat$eta)
xt=data.frame(testdat[,5:(p+4)])
colnames(xt)=paste0("x",1,1:p)
x2t=as.matrix(testdat$x2)
etat=testdat$eta

true1=as.numeric(+0.9*xt[,2]-0.436>=0)
true2=as.numeric(x2t-0.5>=0)
true_generator_dtr2()
discriminator_dtr2(true1,true2)

source("fn_prior.R")
nsim=1000
per1=percen=cra=cr1=cr2=test=testhat=NULL
pb=progress_bar$new(total=nsim)
for(i in 1:nsim)
{
  pb$tick()
  data=data_generator_dtr2(n,p,c0,logit)
  Y=data$Y;delta=data$delta;A1=data$A1;A2=data$A2;eta=data$eta
  mean(delta==0)
  x2=as.matrix(data$x2)
  x=as.matrix(data[,5:(p+4)])
  
  # Stage 2
  X0=x;X1=x2
  svm_mat=cbind(x,x2)
  A=A2
  psurv=drop(pseudosurv(Y,delta,t0)$pseudo)
  padd=data$padd=ifelse(psurv==0,1e-8,psurv)
  Q=Q_model(x0=cbind(X0,A1*X0),x1=X1,y=psurv,trt=A,pen="non")
  
  if (logit=="r") {
    pi1=mean(A==1)
  } else {
    pi1=fitted(glm(as.factor(A2)~x2,family="binomial",data=data))
  }
  
  ipw=as.numeric(with(data,((A==1)*1/pi1-(A!=1)*1/(1-pi1))*padd))
  dr=drop(with(data,ipw-((A==1)*1-pi1)/pi1*Q$Q1-((A!=1)*1-pi1)/(1-pi1)*Q$Q2))
  
  st2_non_ipw=bicsvm2(x=svm_mat,y=Sign(ipw),weight=pmax(abs(ipw*eta),1e-5),type="l2")
  st2_non_dr=bicsvm2(x=svm_mat,y=Sign(dr),weight=pmax(abs(dr*eta),1e-5),type="l2")
  st2_pen_ipw=bicsvm(x=svm_mat,y=Sign(ipw),weight=abs(ipw*eta))
  st2_pen_dr=bicsvm(x=svm_mat,y=Sign(dr),weight=abs(dr*eta))
  
  A2_ni=ifelse(as.numeric(predict(st2_non_ipw))<=1,0,1)
  A2_nd=ifelse(as.numeric(predict(st2_non_dr))<=1,0,1)
  A2_pi=predict(st2_pen_ipw,svm_mat)*1
  A2_pd=predict(st2_pen_dr,svm_mat)*1
  
  new.Y0=psurv+eta*(0-A)*(drop(Q$Q1-Q$Q2))
  new.Y1=psurv+eta*(1-A)*(drop(Q$Q1-Q$Q2))
  new.Y_ni=psurv+eta*(A2_ni-A)*(drop(Q$Q1-Q$Q2))
  new.Y_nd=psurv+eta*(A2_nd-A)*(drop(Q$Q1-Q$Q2))
  new.Y_pi=psurv+eta*(A2_pi-A)*(drop(Q$Q1-Q$Q2))
  new.Y_pd=psurv+eta*(A2_pd-A)*(drop(Q$Q1-Q$Q2))
  
  # Stage 1
  X0=as.matrix(x[,-(2)])
  X1=as.matrix(x[,(2)])
  A=A1
  Q00=Q_model(x0=X0,x1=X1,y=(new.Y0),trt=A,pen="non")
  Q11=Q_model(x0=X0,x1=X1,y=(new.Y1),trt=A,pen="non")
  Qni=Q_model(x0=X0,x1=X1,y=(new.Y_ni),trt=A,pen="non")
  Qnd=Q_model(x0=X0,x1=X1,y=(new.Y_nd),trt=A,pen="non")
  Qpi=Q_model(x0=X0,x1=X1,y=(new.Y_pi),trt=A,pen="non")
  Qpd=Q_model(x0=X0,x1=X1,y=(new.Y_pd),trt=A,pen="non")
  if (logit=="r") {
    pi1=mean(A1==1)
  } else {
    pi1=fitted(glm(as.factor(A1)~x[,2:3],family="binomial",data=data))
  }
  ipw1=as.numeric((A/pi1-(1-A)/(1-pi1))*(new.Y_ni+abs(min(new.Y_ni)-1e-5)))
  ipw2=as.numeric((A/pi1-(1-A)/(1-pi1))*(new.Y_pi+abs(min(new.Y_pi)-1e-5)))
  dr1=drop((A/pi1-(1-A)/(1-pi1))*(new.Y_nd)-
             ((A)-pi1)/pi1*Qnd$Q1-((A)-pi1)/(1-pi1)*Qnd$Q2)
  dr2=drop((A/pi1-(1-A)/(1-pi1))*(new.Y_pd)-
             ((A)-pi1)/pi1*Qpd$Q1-((A)-pi1)/(1-pi1)*Qpd$Q2)
  
  st1_non_ipw=bicsvm2(x=x,y=Sign(ipw1),weight=abs(ipw1),type="l2")
  st1_non_dr=bicsvm2(x=x,y=Sign(dr1),weight=abs(dr1),type="l2")
  st1_pen_ipw=bicsvm(x=x,y=Sign(ipw2),weight=abs(ipw2))
  st1_pen_dr=bicsvm(x=x,y=Sign(dr2),weight=abs(dr2))
  A1_ni=ifelse(as.numeric(predict(st1_non_ipw))<=1,0,1)
  A1_nd=ifelse(as.numeric(predict(st1_non_dr))<=1,0,1)
  A1_pi=predict(st1_pen_ipw,x)*1
  A1_pd=predict(st1_pen_dr,x)*1
  
  testhat=rbind(testhat,
                colMeans(cbind(new.Y0,new.Y1,new.Y_ni,new.Y_nd,new.Y_pi,new.Y_pd)+
                           (cbind(0,1,A1_ni,A1_nd,A1_pi,A1_pd)-A)*
                           (cbind(Q00$Q1-Q00$Q2,Q11$Q1-Q11$Q2,Qni$Q1-Qni$Q2,
                                  Qnd$Q1-Qnd$Q2,Qpi$Q1-Qpi$Q2,Qpd$Q1-Qpd$Q2))))
  
  # New data
  # Stage 1
  newA1_ipw=ifelse(as.numeric(predict(st1_non_ipw,xt))<=1,0,1)
  newA1_dr=ifelse(as.numeric(predict(st1_non_dr,xt))<=1,0,1)
  newA1_pip=predict(st1_pen_ipw,xt)*1
  newA1_pdr=predict(st1_pen_dr,xt)*1
  
  # Stage 2
  newA2_ipw=ifelse(as.numeric(predict(st2_non_ipw,cbind(xt,x2t)))<=1,0,1)
  newA2_dr=ifelse(as.numeric(predict(st2_non_dr,cbind(xt,x2t)))<=1,0,1)
  newA2_pip=predict(st2_pen_ipw,cbind(xt,x2t))*1
  newA2_pdr=predict(st2_pen_dr,cbind(xt,x2t))*1
  
  # Result
  est_trt1=cbind(0,1,newA1_ipw,newA1_dr,newA1_pip,newA1_pdr)
  est_trt2=cbind(0,1,newA2_ipw,newA2_dr,newA2_pip,newA2_pdr)
  cr1=rbind(cr1,colMeans(est_trt1==true1))
  cra=rbind(cra,colMeans((est_trt1[etat==1,]==(matrix(1,nrow(est_trt1),ncol(est_trt1))*true1)[etat==1,])*
                           (est_trt2[etat==1,]==(matrix(1,nrow(est_trt2),ncol(est_trt2))*true2)[etat==1,])))
  
  dsc1=discriminator_dtr2(A1=0,A2=0)
  dsc2=discriminator_dtr2(A1=1,A2=1)
  dsc3=discriminator_dtr2(A1=newA1_ipw,A2=newA2_ipw)
  dsc4=discriminator_dtr2(A1=newA1_dr,A2=newA2_dr)
  dsc5=discriminator_dtr2(A1=newA1_pip,A2=newA2_pip)
  dsc6=discriminator_dtr2(A1=newA1_pdr,A2=newA2_pdr)
  test=rbind(test,c(dsc1,dsc2,dsc3,dsc4,dsc5,dsc6))
  
  tb=cbind(colMeans(test),colMeans(testhat),colMeans(cr1),colMeans(cra))
  colnames(tb)=c("etrue","est","cr1","cra")
  rownames(tb)=c("all0","all1","unpen-ipw","unpen-dr","pen-ipw","pen-dr")
  percen=rbind(percen,mean(delta==0))
  per1=rbind(per1,mean(delta==1))
  if (i>1) {
    print(round(tb,2));
    print(paste0(round(mean(percen),2)*100,"%","Censor&",
                 round(mean(per1),2)*100,"%","Cause-1"))
  }
}
