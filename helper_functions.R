fptcdf=function(z,x0max,chi,driftrate,sddrift) {
  zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu ; xx=chiminuszu-x0max
  chizu=chiminuszu/zs ; chizumax=xx/zs
  tmp1=zs*(dnorm(chizumax)-dnorm(chizu))
  tmp2=xx*pnorm(chizumax)-chiminuszu*pnorm(chizu)
  1+(tmp1+tmp2)/x0max
}

fptpdf=function(z,x0max,chi,driftrate,sddrift) {
  zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu
  chizu=chiminuszu/zs ; chizumax=(chiminuszu-x0max)/zs
  (driftrate*(pnorm(chizu)-pnorm(chizumax)) + 
      sddrift*(dnorm(chizumax)-dnorm(chizu)))/x0max
}

allrtCDF=function(t,x0max,chi,drift,sdI) {
  # Generates CDF for all RTs irrespective of response.
  N=length(drift) # Number of responses.
  tmp=array(dim=c(length(t),N))
  for (i in 1:N) tmp[,i]=fptcdf(z=t,x0max=x0max,chi=chi,driftrate=drift[i],sddrift=sdI)
  1-apply(1-tmp,1,prod)
}

n1PDF=function(t,x0max,chi,drift,sdI) {
  # Generates defective PDF for responses on node #1. 
  N=length(drift) # Number of responses.
  if (N>2) {
    tmp=array(dim=c(length(t),N-1))
    for (i in 2:N) tmp[,i-1]=fptcdf(z=t,x0max=x0max,chi=chi,driftrate=drift[i],sddrift=sdI)
    G=apply(1-tmp,1,prod)
  } else {
    G=1-fptcdf(z=t,x0max=x0max,chi=chi,driftrate=drift[2],sddrift=sdI)
  }
  G*fptpdf(z=t,x0max=x0max,chi=chi,driftrate=drift[1],sddrift=sdI)
}

n1CDF=function(t,x0max,chi,drift,sdI) {
  # Generates defective CDF for responses on node #1. 
  outs=numeric(length(t)) ; bounds=c(0,t)
  for (i in 1:length(t)) {
    tmp="error"
    repeat {
      if (bounds[i]>=bounds[i+1]) {outs[i]=0;break}
      tmp=try(integrate(f=n1PDF,lower=bounds[i],upper=bounds[i+1],
                        x0max=x0max,chi=chi,drift=drift,sdI=sdI)$value,silent=T)
      if (is.numeric(tmp)) {outs[i]=tmp;break}
      # Try smart lower bound.
      if (bounds[i]<=0) {
        bounds[i]=max(c((chi-0.98*x0max)/(max(mean(drift),drift[1])+2*sdI),0))
        next
      }
      # Try smart upper bound.
      if (bounds[i+1]==Inf) {
        bounds[i+1]=0.02*chi/(mean(drift)-2*sdI)
        next
      }
      stop("Error in n1CDF that I could not catch.")
    }
  }
  cumsum(outs)
}

n1mean=function(x0max,chi,drift,sdI) {
  # Generates mean RT for responses on node #1. 
  pc=n1CDF(Inf,x0max,chi,drift,sdI)
  fn=function(t,x0max,chi,drift,sdI,pc) t*n1PDF(t,x0max,chi,drift,sdI)/pc
  tmp=integrate(f=fn,lower=0,upper=100*chi,x0max=x0max,chi=chi,pc=pc,
                drift=drift,sdI=sdI)$value
  list(mean=tmp,p=pc)
}

actCDF=function(z,t,x0max,chi,drift,sdI) {
  # CDF for the distribution (over trials) of activation values at time t.
  zat=(z-x0max)/t ; zt=z/t ; sdi2=2*(sdI^2)
  exp1=exp(-((zt-drift)^2)/sdi2)
  exp2=exp(-((zat-drift)^2)/sdi2)
  tmp1=t*sdI*(exp1-exp2)/sqrt(2*pi)
  tmp2=pnorm(zat,mean=drift,sd=sdI)
  tmp3=pnorm(zt,mean=drift,sd=sdI)
  (tmp1+(x0max-z+drift*t)*tmp2+(z-drift*t)*tmp3)/x0max
}

actPDF=function(z,t,x0max,chi,drift,sdI) {
  # CDF for the distribution (over trials) of activation values at time t.
  tmp1=pnorm((z-x0max)/t,mean=drift,sd=sdI)
  tmp2=pnorm(z/t,mean=drift,sd=sdI)
  (-tmp1+tmp2)/x0max
}  

lbameans=function(Is,sdI,x0max,Ter,chi) {
  # Ter should be a vector of length ncond, the others atomic, 
  # except Is which is ncond x 2.
  ncond=length(Is)/2
  outm<-outp<-array(dim=c(ncond,2))
  for (i in 1:ncond) {
    for (j in 1:2) {
      tmp=n1mean(x0max,chi,drift=Is[i,switch(j,1:2,2:1)],sdI)
      outm[i,j]=tmp$mean+Ter[i]
      outp[i,]=tmp$p
    }}
  list(mns=c(outm[1:ncond,2],outm[ncond:1,1]),ps=c(outp[1:ncond,2],outp[ncond:1,1]))
}

deadlineaccuracy=function(t,x0max,chi,drift,sdI,guess=.5,meth="noboundary") {
  # Works out deadline accuracy, using one of three 
  # methods:
  #   - noboundary = no implicity boundaries.
  #   - partial =  uses implicit boundaries, and partial information.
  #   - nopartial = uses implicit boundaries and guesses otherwise.
  meth=match.arg(meth,c("noboundary","partial","nopartial"))
  noboundaries=function(t,x0max,chi,drift,sdI,ulimit=Inf) {
    # Probability of a correct response in a deadline experiment
    # at times t, with no implicit boundaries.
    N=length(drift)
    tmpf=function(x,t,x0max,chi,drift,sdI) {
      if (N>2) {
        tmp=array(dim=c(length(x),N-1))
        for (i in 2:N) tmp[,i-1]=actCDF(x,t,x0max,chi,drift[i],sdI)
        G=apply(tmp,1,prod)*actPDF(x,t,x0max,chi,drift[1],sdI)
      } else {
        G=actCDF(x,t,x0max,chi,drift[2],sdI)*actPDF(x,t,x0max,chi,drift[1],sdI)
      }}
    outs=numeric(length(t))
    for (i in 1:length(t)) {
      if (t[i]<=0) {
        outs[i]=.5
      } else {
        outs[i]=integrate(f=tmpf,lower=-Inf,upper=ulimit,t=t[i],
                          x0max=x0max,drift=drift,sdI=sdI)$value
      }
    }
    outs
  }
  if (meth=="noboundary") {
    noboundaries(t,x0max,chi,drift,sdI,ulimit=Inf)
  } else { 
    pt=n1CDF(t=t,x0max=x0max,chi=chi,drift=drift,sdI=sdI)
    pa=allrtCDF(t=t,x0max=x0max,chi=chi,drift=drift,sdI=sdI)
    pguess=switch(meth,"nopartial"=guess*(1-pa),"partial"=
                    noboundaries(t=t,x0max=x0max,chi=chi,drift=drift,sdI=sdI,ulimit=chi))
    pt+pguess  
  }
}


# Like pnorm, punif etc evaluates the CDF of the linear BA model at a 
# pre-specified set of quantiles. If input "qs" is null instead
# returns predicted quantiles given by qps.
pqlba = function(Is,sdI,x0max,Ter,chi,qs,qps=seq(.1,.9,.2)) {
  # Input checks.
  if (length(Is)!=2) stop("Not two drifts in pqlba.")
  dop = (!is.null(qs)) # Switch for return p-values or q-values.
  if (dop) {
    nq=dim(qs)[1] # qs should be quants x responses.
    if (length(dim(qs))!=2) stop("Crazy q-values in pqba.")
    if ((dim(qs)[2])!=2) stop("Wrong number of response choices in q-values for pqlba.")
  }
  # First get probability of each response. 
  out=list(p=numeric(2))
  out$pfail=prod(pnorm(-Is/sdI))
  out$p[1]=n1CDF(t=Inf,x0max=x0max,chi=chi,drift=Is,sdI=sdI)
  out$p[2]=n1CDF(t=Inf,x0max=x0max,chi=chi,drift=Is[2:1],sdI=sdI)
  #  out$p[2]=1-out$p[1]
  if (dop) {
    # Calculate probability masses in inter-quantile ranges.
    out$predp=array(dim=c(nq+1,2))
    for (i in 1:2) {
      tmpI=switch(i,Is,Is[2:1]) 
      tmp=n1CDF(t=qs[,i]-Ter,x0max=x0max,chi=chi,drift=tmpI,sdI=sdI)
      out$predp[,i]=diff(c(0,tmp,out$p[i]))
    }
  } else {
    # Calculate predicted quantile values.
    out$predq=array(dim=c(length(qps),2))
    tmpf=function(t,drift,sdI,x0max,chi,p) 
      n1CDF(t=t,x0max=x0max,chi=chi,drift=drift,sdI=sdI)-p
    for (i in 1:2) {
      tmpI=switch(i,Is,Is[2:1]) 
      for (j in 1:length(qps)) {
        interval=switch((Ter<5)+1,c(1,1e4),c(1e-3,10)) # Sec or MSEC.
        tmp=uniroot(f=tmpf,interval=interval,
                    x0max=x0max,chi=chi,drift=tmpI,sdI=sdI,p=qps[j]*out$p[i])
        out$predq[j,i]=tmp$root+Ter
      }
      out$p[i]=n1CDF(t=Inf,x0max=x0max,chi=chi,drift=tmpI,sdI=sdI)
    }
  }
  out$p=out$p+out$pfail/2
  out
}

getpreds=function(I1,I2,sdI,x0max,Ter,chi,pred=F,dat,nc,qps=seq(.1,.9,.2)) {
  # If pred=F, returns chi squared. If pred=T, returns a list
  # like "dat" with predicted quantiles and probabilities.
  # This function expects the parameters (I1,I2,x0max,chi,Ter,sdI)
  # to be 3 arrays (stim1/stim2 x accuracy/neutral/speed).
  mod=list()
  if (pred) {
    mod=list(p=array(dim=c(nc)),q=array(dim=c(length(qps),2,nc)))
  } else {
    mod=list(p=array(dim=c(nc)),q=array(dim=c(length(qps)+1,2,nc)))
  }
  for (j in 1:nc) {
    tmp=pqlba(Is=c(I1[j],I2[j]),sdI=sdI[j],x0max=x0max[j],Ter=Ter[j],chi=chi[j],
              qs=switch(pred+1,dat$q[,,j],NULL))
    mod$p[j]=tmp$p[1]
    mod$q[,,j]=switch(pred+1,tmp$predp,tmp$predq)
  }
  mod
}

obj=function(x,dat,nc,pred=F,qps=seq(.1,.9,.2),trace=0) {
  numits<<-numits+1
  sdI=rep(exp(x[1]),nc) ; x0max=rep(exp(x[2]),nc) ; Ter=rep(exp(x[3]),nc)
  chi=rep(exp(x[4]),nc)+x0max ; I1=(x[4+1:nc])
  I2=1-I1
  preds=getpreds(I1=I1,I2=I2,sdI=sdI,x0max=x0max,Ter=Ter,chi=chi,
                 pred=pred,dat=dat,qps=qps,nc=nc)
  if (pred==F) {
    tmp=-sum(dat$pb*log(pmax(preds$q,1e-10)))
    if (trace>0) if ((numits%%trace)==0) {
      names(x)=NULL ; print(c(exp(x[1:4]),x[4+1:nc],tmp),2)}
    tmp
  } else {
    preds
  }
}

fitter=function(dat,nc,maxit=seq(1000,9000,1000),qps=seq(.1,.9,.2),trace=0) {
  fit=length(maxit)
  I1=qnorm((1-dat$p),sd=.3*sqrt(2))*2
  I1=.5+.5*I1
  tmp=sum1=.5
  Ter=min(dat$q)*.9 ;  x0max=mean(dat$q[4,,]-dat$q[2,,])*(4*tmp)
  chi=(Ter/4) ; sdI=0.3
  par=c(log(c(sdI,x0max,Ter,chi)),I1) ; attr(par,"obj")=NA
  for (fitnum in 1:fit) { # Do the fitting or not.
    numits<<-0
    #return(par)
    tmp=optim(fn=obj,par=par,control=list(maxit=maxit[fitnum],parscale=par),
              dat=dat,qps=qps,trace=trace,nc=nc)
    out<-par<-tmp$par
    attr(out,"obj")=tmp$value
  }
  out
}


fit_data <- function(rawdata){  
  
  #set equal condition for all trials; we are NOT constraining condition
  rawdata$condition <- 1
  
  #set the quantiles used in the QMPE fitting
  qps=seq(.1,.9,.2)
  #which RT values in the data that should be discarded (smaller than the first element, larger than the second)
  trim=c(180,10000)
  
  #rawdata <- fulldata[fulldata$ID==ID,]
  #colnames(rawdata) <- c("ID","correct","rt","difficulty")
  use=(rawdata$rt>trim[1])&(rawdata$rt<trim[2])
  #set up a list of names for variables for later use
  nms=list(c("err","crct"),c("easy"))
  dims=unlist(lapply(nms,length))
  #q gets the values of the quantile values (set by qps above) for correct and error responses for each of the three difficulties
  q=tapply(rawdata$rt[use],list(rawdata$correct[use],rawdata$condition[use]),quantile,probs=qps)
  q=array(unlist(q),dim=c(length(qps),dims),dimnames=c(list(qps),nms))
  #n gets the number of correct and error responses in each difficulty condition
  n=tapply(rawdata$rt[use],list(rawdata$correct[use],rawdata$condition[use]),length)
  #p gets the proportion of correct responses for each difficulty level
  p=tapply(rawdata$correct[use],list(rawdata$condition[use]),mean)
  #gives names to previous variables
  dimnames(n)=nms ; dimnames(p)=nms[-1]
  #calculates how many observations make up each quantile value for correct and error responses in each of the three difficulties
  pb=array(rep(n,each=length(qps)+1),dim=c(length(qps)+1,dim(n)),
           dimnames=c(list(NULL),dimnames(n)))*c(.1,.2,.2,.2,.2,.1)
  #use p, n, q and pb to make a list
  data=list(p=p,n=n,q=q,pb=pb)
  #define the number of drift rate parameters to use
  ndrifts=1
  
  #here is the call to fit the LBA to the model. It requires three arguments - dat, the data, maxit, the number of iterations the fitting algorithm should take, and nc, the number of drift rates to be estimated
  fit=fitter(dat=data,maxit=c(100,200,500,1000,2000,5000),nc=ndrifts)
  #put the results of the fit into a variable called pars
  pars=c(exp(fit[1:4]),fit[5:7])
  #the fitting function gives back (b-A) instead of b. To make the output look like that given in Brown and Heathcote (2008) we add A to b-A to get b
  pars[4]=pars[4]+pars[2]
  #the fitting algorithm returns 1 minus the drift rates for the correct responses. Change these so that they are the actual drift rates
  pars[5:7]=1-pars[5:7]
  fit[5:7]=1-fit[5:7]
  #give names to the parameters
  names(pars)=c("s","A","ter","b","v1","v2","v3")
  #print the parameters onscreen
  #print(pars,3)
  
  return(pars[1:5])
}