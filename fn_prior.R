######################################################
# Title: Accountable survival contrast-learning for  # 
#        optimal dynamic treatment regime            #
# Maintainer: Taehwa Choi                            #
######################################################
svmLP=function(x,y,lambda=1,weight=1){
  if(nrow(x)!=length(y))stop("length(y)mustmatchnrow(x)")
  weight=rep_len(weight,nrow(x))
  y.num=y
  opt=lp(direction="min",
          objective.in=c(rep(lambda,2L*ncol(x)),0,0,weight),
          const.mat=cbind(y.num*x,-y.num*x,y.num,-y.num,diag(nrow(x))),
          const.dir=">=",
          const.rhs=1
  )
  u=opt$solution[seq(1L,length.out=ncol(x))]
  v=opt$solution[seq(ncol(x)+1L,length.out=ncol(x))]
  b=opt$solution[2*ncol(x)+1]-opt$solution[2*ncol(x)+2]
  w=c(b,u-v)
  class(w)="svmLP"
  w
}

svmscad=function(x,y,weight,lambda){
  init=svmLP(x=x,y=y,weight=weight,lambda=0)
  penval=pen(x=init[-1],lam=lambda)
  
  if(nrow(x)!=length(y))stop("length(y)mustmatchnrow(x)")
  weight=rep_len(weight,nrow(x))
  y.num=y
  opt=lp(direction="min",
          objective.in=c(rep(penval,2),0,0,weight),
          const.mat=cbind(y.num*x,-y.num*x,y.num,-y.num,diag(nrow(x))),
          const.dir=">=",
          const.rhs=1
  )
  u=opt$solution[seq(1L,length.out=ncol(x))]
  v=opt$solution[seq(ncol(x)+1L,length.out=ncol(x))]
  b=opt$solution[2*ncol(x)+1]-opt$solution[2*ncol(x)+2]
  w=c(b,u-v)
  class(w)="svmLP"
  w
}

bicsvm=function(x,y,weight,ngrid=21){
  grid.lambda=2^((-(ngrid-1)/2):((ngrid-1)/2))
  est=sapply(grid.lambda,svmscad,x=x,y=y,weight=weight)
  bicval=0
  for(i in 1:length(grid.lambda)){
    bicval[i]=log(sum(weight*pmax(1-y*cbind(1,x)%*%est[,i],0)))
  }
  svmscad(x=x,y=y,weight=weight,lambda=grid.lambda[which.min(bicval)])
}

predict.svmLP=function(object,x,...){
  f=as.matrix(cbind(1,x))%*%object
  y=f>0
  return(y)
}

pen=function(x,lam,a=3.7)
{
  lam*(abs(x)<=lam)+pmax(a*lam-abs(x),0)/(a-1)*(abs(x)>lam)
}

nonsvm=function(x,y,lambda,weight){
  mod=svm(x,y,type="C",kernel="linear",cost=weight*lambda)
  list(est=c(-mod$rho,c(mod$coefs)%*%mod$SV))
}

bicsvm2=function(x,y,weight,nfold=5,ngrid=11,type=c("scad","l2")){
  y=c(y);weight=c(weight)
  n=nrow(x)
  gridlam=2^((-(ngrid-1)/2):((ngrid-1)/2))
  foldid=sample(nfold,n,replace=TRUE)
  V=0
  if (type=="scad"){
    for(i in 1:ngrid){
      fit=pensvm(x=x,y=y,weight=weight,lambda=gridlam[i])
      V[i]=n*log(mean((weight*pmax(0,1-cbind(1,x)%*%fit$est*y))))+log(n)*sum(fit$est[-1]!=0)
    }
    lam=gridlam[which.min(V)]
    fit=pensvm(x=x,y=y,weight=weight,lambda=lam)
  } else if (type=="l2"){
    for(i in 1:ngrid){
      fit=nonsvm(x=x,y=y,weight=weight,lambda=gridlam[i])
      V[i]=n*log(mean((weight*pmax(0,1-cbind(1,x)%*%fit$est*y))))+log(n)*sum(fit$est[-1]!=0)
    }
    lam=gridlam[which.min(V)]
    fit=svm(x=x,y=y,type="C",kernel="linear",cost=weight*lam)
  }
  fit
}

Q_model=function(x0,x1,y,trt,pen=c("non","lasso"))
{
  Xmat=as.matrix(cbind(x0,x1*trt))
  if(pen=="lasso")
  {
    ind1=(trt==1);ind2=(trt!=1)
    Xmat1=Xmat[ind1,];Y1=y[ind1]
    Xmat2=Xmat[ind2,];Y2=y[ind2]
    Q1=cbind(1,Xmat)%*%as.numeric(predict(cv.glmnet(Xmat1,Y1,nfolds=5),s="lambda.min",type="coeff"))
    Q2=cbind(1,Xmat)%*%as.numeric(predict(cv.glmnet(Xmat2,Y2,nfolds=5),s="lambda.min",type="coeff"))
    
  }else{
    id1=(trt==1);id2=(trt!=1)
    lmfit1=coef(lm(y~Xmat))
    lmfit1=ifelse(is.na(lmfit1),0,lmfit1)
    
    Q1=cbind(1,as.matrix(cbind(x0,x1*1)))%*%lmfit1
    Q2=cbind(1,as.matrix(cbind(x0,x1*0)))%*%lmfit1
    
  }
  return(list(Q1=Q1,Q2=Q2))
}

expit=function(x)exp(x)/(1+exp(x))

weight_fun=function(y,pi1,trt,Q,type)
{
  ipw=y/ifelse(A==1,pi1,1-pi1)
  dr=(A/pi1-(1-A)/(1-pi1))*y-(A-pi1)/pi1*Q$Q1-(A-pi1)/(1-pi1)*Q$Q2
  ifelse(type=="ipw",return(ipw),return(dr))
}

Sign=function(x)ifelse(x>0,1,-1)

