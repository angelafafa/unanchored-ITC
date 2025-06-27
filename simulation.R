library(WeightIt)
library(copula)
library(WeightIt)
library(copula)
#main simulation of four model specifications for bias
#linear=TRUE, outcome="nbinary";linear=TRUE, outcome="binary"
#linear=FALSE,outcome="nbinary";linear=FALSE,outcome="binary"
sim2_norm=function(ORgood=1,PSgood=1,n=2000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=FALSE,shift1=FALSE,shift2=FALSE,shift_param=0.1,shift_param1=15,shift_param2=0.5){
  #coef for outcome
  theta<-c(2,0.6,1.75)
  #Generate data
  if(linear){
    #Scenario 1 -- when X2=G=g~bernoulli
    cutoff=0.99
    al<-c(-1.75,0.3,0.5)
    G1=G0=g0=g1=rbinom(n,1,0.25)
    #x0=x1<-runif(n,0,10)
    x0=x1<-rnorm(n,5,sqrt(100/12))
    #x1=rnorm(n,5,2)
    error=rnorm(n)
    pi.vec0=pi.vec1<-1/(1+exp(-cbind(1,x1,g1) %*% al))
    Z0=Z1<-rbinom(n,1,pi.vec1)
    #Shift in S1
    if(shift){
      if(shift1){
        #covariate shift for untreated population Z=0
        x0=rnorm(n,shift_param,shift_param1)
      }
      if(shift2){
        g0=G0=rbinom(n,1,shift_param2)
      }
      #treatment assignment after shifting
      pi.vec0<-1/(1+exp(-cbind(1,x0,g0) %*% al))
      Z0<-rbinom(n,1,pi.vec0)
    }
    #Y[Z==1]=Y1[Z==1]
    #Y[Z==0]=Y0[Z==0]
    #Y<-cbind(Z,Z*x,G)%*%theta+x+error
    Y1<-cbind(1,x1,G1)%*%theta+x1+error
    Y0<-cbind(0,0,G0)%*%theta+x0+error
    Y1_0<-cbind(0,0,G1)%*%theta+x1+error
  }else{
    #Scenario 2 -- when X2=G1~bernoulli and X3=g1~Normal
    cutoff=0.5
    #al<-c(-1.75,0.3,0.01)
    al<-c(-1.75,0.3,0.1)
    #Instrumental variable X3=g0=g1
    #g0=g1=rnorm(n,mean = 56,sd=12)
    g0=g1=runif(n,0,10)
    G0=G1=rbinom(n,1,0.25)
    #x0=x1<-runif(n,0,10)
    x0=x1<-rnorm(n,mean=5,sd=sqrt(100/12))
    error=rnorm(n)
    pi.vec0=pi.vec1<-1/(1+exp(-cbind(1,x1,g1) %*% al))
    Z0=Z1<-rbinom(n,1,pi.vec1)
    if(shift){
      if(shift1){
        #strength of IV
        al[3]=shift_param
      }
      if(shift2){
        #G0=rbinom(n,1,shift_param2)
        x0=runif(n,shift_param2,shift_param2+10)
        #g0=rnorm(n,mean = 60,sd=10)
      }
      pi.vec0<-1/(1+exp(-cbind(1,x0,g0) %*% al))
      Z0<-rbinom(n,1,pi.vec0)
    }
    #Y<-cbind(Z,Z*x,G)%*%theta+x-0.2*x^2+error
    Y1<-cbind(1,x1,G1)%*%theta+x1-0.2*x1^2+error
    Y0<-cbind(0,0,G0)%*%theta+x0-0.2*x0^2+error
    Y1_0<-cbind(0,0,G1)%*%theta+x1-0.2*x1^2+error
  }
    Y=c(Y1[Z1==1][1:n_A],Y0[Z0==0][1:n_B])
    Z=c(rep(1,n_A),rep(0,n_B))
    x=c(x1[Z1==1][1:n_A],x0[Z0==0][1:n_B])
    g=c(g1[Z1==1][1:n_A],g0[Z0==0][1:n_B])
    G=c(G1[Z1==1][1:n_A],G0[Z0==0][1:n_B])
    att=mean(Y1[Z1==1][1:n_A])-mean(Y1_0[Z1==1][1:n_A])
    
    #For binary outcome
    if(outcome=="binary"){
      p=exp(Y)/(1+exp(Y))
      y=ifelse(p>cutoff,1,0)
      p0=exp(Y1_0)/(1+exp(Y1_0))
      y_0=ifelse(p0>cutoff,1,0)
      att=mean(y[1:n_A])-mean(y_0[Z1==1][1:n_A])
      Y=y
    }
    data=data.frame(Y,Z,x,x^2,x^4,Z*x,g,g^2,G,G^2)
    dat0 <- dat <- data.frame(Y=data$Y,Z=data$Z,x=data$x,g=data$g,G=data$G); dat0$Z<-0
    dat0_1=dat[dat$Z==1,];dat0_1$Z<-0
    dat0_0=dat[dat$Z==0,]
    
    if(PSgood==1){
      pi.hat<-fitted(glm(Z~x+g,family=binomial,data=data))
    }
    else{
      if(linear){pi.hat<-fitted(glm(Z~x,family=binomial,data=data))}
      else{pi.hat<-fitted(glm(Z~g,family=binomial,data=data))}
      }
    data$w.hat<-pi.hat/(1-pi.hat)
    data$ps=ifelse(data$Z==1,1,data$w.hat)
    #IPW
    sr.att <-sum((data$Z-(1-data$Z)*data$w.hat)*data$Y)/sum(data$Z)
    #Stabilized IPW normalizing weights
    sr_att=sum(data$ps*data$Z*data$Y)/sum(data$ps*data$Z)-sum(data$ps*(1-data$Z)*data$Y)/sum(data$ps*(1-data$Z))
    #Outcome regression (one unified model with treatment indicator)
    
    #MAIC
    objfn <- function(a1, X){
      sum(exp(X %*% a1))
    }
    gradfn <- function(a1, X){
      colSums(sweep(X, 1, exp(X %*% a1), "*"))
    }
    mp1.mean=apply(data[data$Z==1,c("x","x.2","g","g.2","G","G.2")],2,mean)
    mp1.sd=apply(data[data$Z==1,c("x","x.2","g","g.2","G","G.2")],2,sd)
    if(linear){
      if(PSgood==1) {
        #center both mean and standard deviation
        X.EM.0_cen=data.matrix(data.frame(x_cen=data[data$Z==0,c("x")]-mp1.mean[[1]],
                                          x.2_cen=data[data$Z==0,c("x.2")]-((mp1.mean[[1]])^2+(mp1.sd[[1]])^2),
                                          g_cen=data[data$Z==0,c("g")]-mp1.mean[[3]],
                                          g.2_cen=data[data$Z==0,c("g.2")]-((mp1.mean[[3]])^2+(mp1.sd[[3]])^2)
        ))
        W2 <- weightit(Z ~ x+g, data = data,method = "ebal", estimand = "ATT")
        opt1 <- optim(par = c(0,0,0,0), fn = objfn, gr = gradfn, X = X.EM.0_cen, method = "BFGS")
      }
      else{
        X.EM.0_cen=data.matrix(data.frame(
          x_cen=data[data$Z==0,c("x")]-mp1.mean[[1]],
          x.2_cen=data[data$Z==0,c("x.2")]-((mp1.mean[[1]])^2+(mp1.sd[[1]])^2)
        ))
        W2 <- weightit(Z ~ x, data = data,method = "ebal", estimand = "ATT")
        opt1 <- optim(par = c(0,0), fn = objfn, gr = gradfn, X = X.EM.0_cen, method = "BFGS")
      }
    }else{
      if(PSgood==1) {
        #center both mean and standard deviation
        X.EM.0_cen=data.matrix(data.frame(x_cen=data[data$Z==0,c("x")]-mp1.mean[[1]],
                                          x.2_cen=data[data$Z==0,c("x.2")]-((mp1.mean[[1]])^2+(mp1.sd[[1]])^2),
                                          x_2_cen=data[data$Z==0,c("x.2")]-mp1.mean[[2]],
                                          x_2.2_cen=data[data$Z==0,c("x.4")]-((mp1.mean[[2]])^2+(mp1.sd[[2]])^2),
                                          g_cen=data[data$Z==0,c("G")]-mp1.mean[[5]],
                                          g.2_cen=data[data$Z==0,c("G.2")]-((mp1.mean[[5]])^2+(mp1.sd[[5]])^2)
        ))
        W2 <- weightit(Z ~ x+I(x^2)+G, data = data,method = "ebal", estimand = "ATT")
        opt1 <- optim(par = c(0,0,0,0,0,0), fn = objfn, gr = gradfn, X = X.EM.0_cen, method = "BFGS")
      }
      else{
        X.EM.0_cen=data.matrix(data.frame(g_cen=data[data$Z==0,c("G")]-mp1.mean[[5]],
                                          g.2_cen=data[data$Z==0,c("G.2")]-((mp1.mean[[5]])^2+(mp1.sd[[5]])^2)
        ))
        W2 <- weightit(Z ~ G, data = data,method = "ebal", estimand = "ATT")
        opt1 <- optim(par = c(0,0), fn = objfn, gr = gradfn, X = X.EM.0_cen, method = "BFGS")
      }
    }
    a1 <- opt1$par
    #weights of MAIC
    wt <- exp(X.EM.0_cen %*% a1)/sum(exp(X.EM.0_cen %*% a1))
    maic_att=mean(data[data$Z==1,"Y"])-sum(data[data$Z==0,"Y"]*wt)
    #balance weight-based methods
    wt_eb=W2$weights[data$Z==0]/(sum(W2$weights[data$Z==0]))
    eb_att=mean(data[data$Z==1,"Y"])-sum(wt_eb*data[data$Z==0,"Y"])
    
    if(linear){
      if(ORgood==1){
        if(outcome=="binary"){
          mu0 <- predict(glm(Y~Z*x+G,data=dat,family=binomial()),dat0,type="response")
          fit_all1=predict(glm(Y~Z*x+G,data=data,family=binomial()),newdata = dat0_1,type="response")
          fit_all0=predict(glm(Y~Z*x+G,data=data,family=binomial()),newdata = dat0_0,type="response")
          #STC
          fit=glm(Y~x+G,data=data[data$Z==0,],family=binomial())
          #STC_park
          fit_park=glm(Y~x+G,data=data[data$Z==0,],family=binomial(),weights = wt_eb)
          # simulate covariates for the AgD study
          ipd_cor <- cor(data[data$Z==0, c("x", "G")])
          ipd_copula <-
            copula::normalCopula(copula::P2p(ipd_cor),
                                 dim = ncol(ipd_cor),
                                 dispstr = "un")
     
          myMvd_AgD <-
            mvdc(
              copula = ipd_copula,
              margins = c("norm","binom"),
              paramMargins = list(list(mean = mean(data[data$Z==1,"x"]), sd = sd(data[data$Z==1,"x"])),
                                  list(size = 1, prob = mean(data[data$Z==1,"G"])))
            )
          mycov_AgD <- data.frame(rMvdc(sum(data$Z==1), myMvd_AgD))
          colnames(mycov_AgD)=c("x", "G")
          mycov_AgD$Z=1
        }
        #regression adjustment w/ one unified model
        else{
        mu0 <- predict(lm(Y~Z*x+G,data=dat),dat0)
        fit_all1=predict(lm(Y~Z*x+G,data=data),newdata = dat0_1)
        fit_all0=predict(lm(Y~Z*x+G,data=data),newdata = dat0_0)
        #STC
        fit=lm(Y~x+G,data=data[data$Z==0,])
        fit_park=lm(Y~x+G,data=data[data$Z==0,],weights = wt_eb)}
        }
      else{
        if(outcome=="binary"){
          mu0 <- predict(glm(Y~Z*x,data=dat,family=binomial()),dat0,type="response")
          fit_all1=predict(glm(Y~Z*x,data=data,family=binomial()),newdata = dat0_1,type="response")
          fit_all0=predict(glm(Y~Z*x,data=data,family=binomial()),newdata = dat0_0,type="response")
          #STC
          fit=glm(Y~x,data=data[data$Z==0,],family=binomial())
          fit_park=glm(Y~x,data=data[data$Z==0,],family=binomial(),weights = wt_eb)
          # simulate covariates for the AgD study
          mycov_AgD <- data.frame(Z=rep(1,sum(data$Z==1)), x=rnorm(sum(data$Z==1),mean = mean(data[data$Z==1,"x"]), sd = sd(data[data$Z==1,"x"])))
        }
        else{
        mu0 <- predict(lm(Y~Z*x,data=dat),dat0)
        fit_all1=predict(lm(Y~Z*x,data=data),newdata = dat0_1)
        fit_all0=predict(lm(Y~Z*x,data=data),newdata = dat0_0)
        fit=lm(Y~x,data=data[data$Z==0,])
        fit_park=lm(Y~x,data=data[data$Z==0,],weights = wt_eb)
       }
        }
    }else{
      if(ORgood==1){
        if(outcome=="binary"){
          mu0 <- predict(glm(Y~Z*x+I(x^2)+G,data=dat,family = binomial()),dat0,type="response")
          fit_all1=predict(glm(Y~Z*x+I(x^2)+G,data=data,family = binomial()),newdata = dat0_1,type="response")
          fit_all0=predict(glm(Y~Z*x+I(x^2)+G,data=data,family = binomial()),newdata = dat0_0,type="response")
          #STC
          fit=glm(Y~x+I(x^2)+G,data=data[data$Z==0,],family = binomial())
          fit_park=glm(Y~x+I(x^2)+G,data=data[data$Z==0,],family = binomial(),weights = wt_eb)
          # simulate covariates for the AgD study
          ipd_cor <- cor(data[data$Z==0, c("x", "x.2","G")])
          ipd_copula <-
            copula::normalCopula(copula::P2p(ipd_cor),
                                 dim = ncol(ipd_cor),
                                 dispstr = "un")
          myMvd_AgD <-
            mvdc(
              copula = ipd_copula,
              margins = c("norm","norm","binom"),
              paramMargins = list(list(mean = mean(data[data$Z==1,"x"]), sd = sd(data[data$Z==1,"x"])),
                                  list(mean = mean(data[data$Z==1,"x.2"]), sd = sd(data[data$Z==1,"x.2"])),
                                  list(size = 1, prob = mean(data[data$Z==1,"G"])))
            )
          mycov_AgD <- data.frame(rMvdc(sum(data$Z==1), myMvd_AgD))
          colnames(mycov_AgD)=c("x","x.2","G")
          mycov_AgD$Z=1
        }
        else{mu0 <- predict(lm(Y~Z*x+I(x^2)+G,data=dat),dat0)
        fit_all1=predict(lm(Y~Z*x+I(x^2)+G,data=data),newdata = dat0_1)
        fit_all0=predict(lm(Y~Z*x+I(x^2)+G,data=data),newdata = dat0_0)
        #STC
        fit=lm(Y~x+I(x^2)+G,data=data[data$Z==0,])
        fit_park=lm(Y~x+I(x^2)+G,data=data[data$Z==0,],weights = wt_eb)
       }
        }
      else{
        if(outcome=="binary"){
          mu0 <- predict(glm(Y~Z+G,data=dat,family = binomial()),dat0,type="response")
          fit_all1=predict(glm(Y~Z+G,data=data,family = binomial()),newdata = dat0_1,type="response")
          fit_all0=predict(glm(Y~Z+G,data=data,family = binomial()),newdata = dat0_0,type="response")
          #STC
          fit=glm(Y~G,data=data[data$Z==0,],family = binomial())
          fit_park=glm(Y~G,data=data[data$Z==0,],family = binomial(),weights = wt_eb)
          # simulate covariates for the AgD study
          mycov_AgD <- data.frame(Z=rep(1,sum(data$Z==1)), G=rbinom(n=sum(data$Z==1),size=1,prob = mean(data[data$Z==1,"G"])))
        }
        else{mu0 <- predict(lm(Y~Z+G,data=dat),dat0)
        fit_all1=predict(lm(Y~Z+G,data=data),newdata = dat0_1)
        fit_all0=predict(lm(Y~Z+G,data=data),newdata = dat0_0)
        fit=lm(Y~G,data=data[data$Z==0,])
        fit_park=lm(Y~G,data=data[data$Z==0,],weights = wt_eb)
       }
        }
    }
    OR.att=mean(data[data$Z==1,"Y"])-mean(fit_all1)
    ## Doubly-robust estimation
    #AIPW
    dr.att.norm=(mean(data[data$Z==1,"Y"])-mean(fit_all1))-sum(data[data$Z==0,"ps"]*(data[data$Z==0,"Y"]-fit_all0))/sum(data$ps[data$Z==0])
    if(outcome=="binary"){
      mu0_stc <- predict(fit,newdata=mycov_AgD,type="response")
      stc_pred=mean(mu0_stc)
      mu0_stc_park <- predict(fit_park,newdata=mycov_AgD,type="response")
      stc_pred_park=mean(mu0_stc_park)
      #unweighted AIPW using separate regression model
      aipw_att=(mean(data[data$Z==1,"Y"])-mean(predict(fit,newdata=data[data$Z==1,],type="response")))-sum(data$ps[data$Z==0]*(data[data$Z==0,"Y"]-predict(fit,type="response")))/sum(data$Z)
      #AIPW using separate regression model
      aipw_norm_att=(mean(data[data$Z==1,"Y"])-mean(predict(fit,newdata=data[data$Z==1,],type="response")))-sum(data$ps[data$Z==0]*(data[data$Z==0,"Y"]-predict(fit,type="response")))/sum(data$ps[data$Z==0])
      #double robust estimator for aggregate data
      dr_att=(mean(data[data$Z==1,"Y"])-stc_pred)-sum(wt*(data[data$Z==0,"Y"]-predict(fit,type="response")))
      dr_att2=(mean(data[data$Z==1,"Y"])-stc_pred)-sum(wt_eb*(data[data$Z==0,"Y"]-predict(fit,type="response")))
    }else{
      mu0_stc <- predict(fit,newdata=data[data$Z==1,])
      stc_pred=mean(mu0_stc)
      mu0_stc_park <- predict(fit_park,newdata=data[data$Z==1,])
      stc_pred_park=mean(mu0_stc_park)
      dr_att2=(mean(data[data$Z==1,"Y"])-stc_pred)-sum(wt_eb*(data[data$Z==0,"Y"]-predict(fit)))
    }
    stc_att=mean(data[data$Z==1,"Y"])-stc_pred
    park_att=mean(data[data$Z==1,"Y"])-stc_pred_park
    summary.bias=c(sr.att-att,sr_att-att,eb_att-att,maic_att-att,OR.att-att,stc_att-att,dr.att.norm-att,dr_att2-att,park_att-att)
    return(summary.bias)
}
#bias for table 1
seeds <- 1:1000
count <- 0
result1_4=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=1,PSgood=1,n=3000,n_A=200,n_B=200,linear=FALSE,outcome="binary",shift=FALSE)})
seeds <- 1:1000
count <- 0
result2_4=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=0,PSgood=1,n=3000,n_A=200,n_B=200,linear=FALSE,outcome="binary",shift=FALSE)})
seeds <- 1:1000
count <- 0
result3_4=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=1,PSgood=0,n=3000,n_A=200,n_B=200,linear=FALSE,outcome="binary",shift=FALSE)})
seeds <- 1:1000
count <- 0
result4_4=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=0,PSgood=0,n=3000,n_A=200,n_B=200,linear=FALSE,outcome="binary",shift=FALSE)})
seeds <- 1:2000
count <- 0
result1_2=replicate(2000, { count <<- count + 1; set.seed(seeds[count]); sim2(ORgood=1,PSgood=1,n=3000,n_A=1000,n_B=1000,linear=TRUE,outcome="binary",shift=FALSE)})
seeds <- 1:2000
count <- 0
result2_2=replicate(2000, { count <<- count + 1; set.seed(seeds[count]); sim2(ORgood=0,PSgood=1,n=3000,n_A=1000,n_B=1000,linear=TRUE,outcome="binary",shift=FALSE)})
seeds <- 1:2000
count <- 0
result3_2=replicate(2000, { count <<- count + 1; set.seed(seeds[count]); sim2(ORgood=1,PSgood=0,n=3000,n_A=1000,n_B=1000,linear=TRUE,outcome="binary",shift=FALSE)})
seeds <- 1:2000
count <- 0
result4_2=replicate(2000, { count <<- count + 1; set.seed(seeds[count]); sim2(ORgood=0,PSgood=0,n=3000,n_A=1000,n_B=1000,linear=TRUE,outcome="binary",shift=FALSE)})
df3_2=data.frame(rbind(apply(result1_2, 1, mean),apply(result2_2, 1, mean),apply(result3_2, 1, mean),apply(result4_2, 1, mean)))
df3_2_sd=data.frame(rbind(apply(result1_2, 1, sd),apply(result2_2, 1, sd),apply(result3_2, 1, sd),apply(result4_2, 1, sd)))
df3_2$N="n=1000"
df3_2$outcome="Binary (linear)"
df3_2$Misspecification=c("Neither","OR","PS","Both")
colnames(df3_2)=c("MAIC1","MAIC2","IPW1","IPW2","RA","STC","AIPW","DR","DR (park)","N","Misspecification")
df3_2_sd$Misspecification=c("Neither","OR","PS","Both")
colnames(df3_2_sd)=c("MAIC1","MAIC2","IPW1","IPW2","RA","STC","AIPW","DR","DR (park)","Misspecification")
rownames(df3_2)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
seeds <- 1:1000
count <- 0
result1_4=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=1,PSgood=1,n=3000,n_A=500,n_B=500,linear=FALSE,outcome="binary",shift=FALSE)})
seeds <- 1:1000
count <- 0
result2_4=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=0,PSgood=1,n=3000,n_A=500,n_B=500,linear=FALSE,outcome="binary",shift=FALSE)})
seeds <- 1:1000
count <- 0
result3_4=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=1,PSgood=0,n=3000,n_A=500,n_B=500,linear=FALSE,outcome="binary",shift=FALSE)})
seeds <- 1:1000
count <- 0
result4_4=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=0,PSgood=0,n=3000,n_A=500,n_B=500,linear=FALSE,outcome="binary",shift=FALSE)})
df1_4=data.frame(rbind(apply(result1_4, 1, mean),apply(result2_4, 1, mean),apply(result3_4, 1, mean),apply(result4_4, 1, mean)))
df1_4_sd=data.frame(rbind(apply(result1_4, 1, sd),apply(result2_4, 1, sd),apply(result3_4, 1, sd),apply(result4_4, 1, sd)))
df1_4$N="n=200 per study"
df1_4$Misspecification=c("Neither","OR","PS","Both")
colnames(df1_4)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","N","Misspecification")
df1_4_sd$Misspecification=c("Neither","OR","PS","Both")
colnames(df1_4_sd)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","Misspecification")
rownames(df1_4)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
df2_4=data.frame(rbind(apply(result1_4, 1, mean),apply(result2_4, 1, mean),apply(result3_4, 1, mean),apply(result4_4, 1, mean)))
df2_4_sd=data.frame(rbind(apply(result1_4, 1, sd),apply(result2_4, 1, sd),apply(result3_4, 1, sd),apply(result4_4, 1, sd)))
df2_4$N="n=500 per study"
df2_4$Misspecification=c("Neither","OR","PS","Both")
colnames(df2_4)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","N","Misspecification")
df2_4_sd$Misspecification=c("Neither","OR","PS","Both")
colnames(df2_4_sd)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","Misspecification")
rownames(df2_4)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")

df3_4=data.frame(rbind(apply(result1_4, 1, mean),apply(result2_4, 1, mean),apply(result3_4, 1, mean),apply(result4_4, 1, mean)))
df3_4_sd=data.frame(rbind(apply(result1_4, 1, sd),apply(result2_4, 1, sd),apply(result3_4, 1, sd),apply(result4_4, 1, sd)))
df3_4$N="n=1000 per study"
df3_4$Misspecification=c("Neither","OR","PS","Both")
colnames(df3_4)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","N","Misspecification")
df3_4_sd$Misspecification=c("Neither","OR","PS","Both")
colnames(df3_4_sd)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","Misspecification")
rownames(df3_4)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")



df_1=rbind(df1_1,df2_1,df3_1)
df_1$outcome="Linear"
df_2=rbind(df1_2,df2_2,df3_2)
df_2$outcome="Binary (linear)"
df_3=rbind(df1_3,df2_3,df3_3)
df_3$outcome="Nonlinear"
df_4=rbind(df1_4,df2_4,df3_4)
df=rbind(df_1,df_2,df_3,df_4)
df_4$outcome="Binary (nonlinear)"
df_1_sd=rbind(df1_1_sd,df2_1_sd,df3_1_sd)
df_2_sd=rbind(df1_2_sd,df2_2_sd,df3_2_sd)
df_3_sd=rbind(df1_3_sd,df2_3_sd,df3_3_sd)
df_4_sd=rbind(df1_4_sd,df2_4_sd,df3_4_sd)
df_sd=rbind(df_1_sd,df_2_sd,df_3_sd,df_4_sd)
library(reshape2)
library(ggplot2)
save(df,file="bias_3ss.RData")
save(df_sd,file="bias_sd_3ss.RData")
# Melt to long format
df_long <- melt(df,
                id.vars = c("Misspecification", "N", "outcome"),
                variable.name = "Method",
                value.name = "Estimate")
df_long_sd <- melt(df_sd,
                id.vars = c("Misspecification"),
                variable.name = "Method",
                value.name = "SE")
df_long$CI_lower <- df_long$Estimate - 1.96 * df_long_sd$SE
df_long$CI_upper <- df_long$Estimate + 1.96 * df_long_sd$SE

df_long$N=ifelse(df_long$N=="n=200","n=200 per study",ifelse(df_long$N=="n=500","n=500 per study","n=1000 per study"))
  
df_long$N=factor(df_long$N,levels = c("n=200 per study","n=500 per study","n=1000 per study"))
#n_A=n_B=200,n_A=n_B=500,n_A=n_B=1000

ggplot(df_long, aes(x = Method, y = Estimate, color = Misspecification)) +
  #geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  #geom_bar(stat = "identity", position = "dodge") +
  geom_point(position = position_dodge(width=0.9),size=2)+
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2, position = position_dodge(width = 0.9)) +
  facet_grid(outcome ~ N,scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  theme_minimal(base_size = 13) +
  labs(x = "Method", y = "Bias") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))
#bootstrap SE
sim2.BB=function(BB=50,seed=1,ORgood=1,PSgood=1,n=3000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=FALSE,shift1=FALSE,shift2=FALSE,shift_param=0.5,shift_param1=0.1,shift_param2=0.5){
  set.seed(seed)
  out=NULL
  #coef for outcome
  theta<-c(2,0.6,1.75)
  #att <- 5.749209
  #att=2.312511 if x~unif(0,1)
  #Generate data
  if(linear){
    #Scenario 1 -- when X2=G=g~bernoulli
    cutoff=0.99
    al<-c(-1.75,0.3,0.5)
    G1=G0=g0=g1=rbinom(n,1,0.25)
    #x0=x1<-runif(n,0,10)
    x0=x1<-rnorm(n,5,sqrt(100/12))
    error=rnorm(n)
    pi.vec0=pi.vec1<-1/(1+exp(-cbind(1,x1,g1) %*% al))
    Z0=Z1<-rbinom(n,1,pi.vec1)
    if(shift){
      #covariate shift for untreated population Z=0
      if(shift1){
        x0=runif(n,shift_param,shift_param1)
      }
      if(shift2){
        g0=G0=rbinom(n,1,shift_param2)
      }
      #treatment assignment after shifting
      pi.vec0<-1/(1+exp(-cbind(1,x0,g0) %*% al))
      Z0<-rbinom(n,1,pi.vec0)
    }
    #Y[Z==1]=Y1[Z==1]
    #Y[Z==0]=Y0[Z==0]
    #Y<-cbind(Z,Z*x,G)%*%theta+x+error
    Y1<-cbind(1,x1,G1)%*%theta+x1+error
    Y0<-cbind(0,0,G0)%*%theta+x0+error
    Y1_0<-cbind(0,0,G1)%*%theta+x1+error
  }else{
    #Scenario 2 -- when X2=G1~bernoulli and X3=g1~Normal
    cutoff=0.5
    al<-c(-1.75,0.3,0.01)
    #Instrumental variable X3=g0=g1
    g0=g1=runif(n,0,10)
    #g0=g1=rnorm(n,mean = 56,sd=12)
    G0=G1=rbinom(n,1,0.25)
    #x0=x1<-runif(n,0,10)
    x0=x1<-rnorm(n,mean=5,sd=sqrt(100/12))
    error=rnorm(n)
    pi.vec0=pi.vec1<-1/(1+exp(-cbind(1,x1,g1) %*% al))
    Z0=Z1<-rbinom(n,1,pi.vec1)
    if(shift){
      if(shift1){
        #strength of IV
        al[3]=shift_param
      }
      if(shift2){
        #G0=rbinom(n,1,shift_param2)
        x0=rnorm(n,mean=shift_param2,sd=sqrt(100/12))
      }
      pi.vec0<-1/(1+exp(-cbind(1,x0,g0) %*% al))
      Z0<-rbinom(n,1,pi.vec0)
    }
    
    #Y<-cbind(Z,Z*x,G)%*%theta+x-0.2*x^2+error
    Y1<-cbind(1,x1,G1)%*%theta+x1-0.2*x1^2+error
    Y0<-cbind(0,0,G0)%*%theta+x0-0.2*x0^2+error
    Y1_0<-cbind(0,0,G1)%*%theta+x1-0.2*x1^2+error
  }
  Y=c(Y1[Z1==1][1:n_A],Y0[Z0==0][1:n_B])
  Z=c(rep(1,n_A),rep(0,n_B))
  x=c(x1[Z1==1][1:n_A],x0[Z0==0][1:n_B])
  g=c(g1[Z1==1][1:n_A],g0[Z0==0][1:n_B])
  G=c(G1[Z1==1][1:n_A],G0[Z0==0][1:n_B])
  att=mean(Y1[Z1==1][1:n_A])-mean(Y1_0[Z1==1][1:n_A])
  
  #For binary outcome
  if(outcome=="binary"){
    p=exp(Y)/(1+exp(Y))
    y=ifelse(p>cutoff,1,0)
    p0=exp(Y1_0)/(1+exp(Y1_0))
    y_0=ifelse(p0>cutoff,1,0)
    #att=mean(y[Z==1])-mean(y_0[Z1==1])
    att=mean(y[1:n_A])-mean(y_0[Z1==1][1:n_A])
    Y=y
  }
  data1=data.frame(Y,Z,x,x^2,x^4,Z*x,g,g^2,G,G^2)
  bootstrap_index0=replicate(BB,sample(1:sum(data1$Z==0), replace=TRUE))
  bootstrap_index1=replicate(BB,sample(1:sum(data1$Z==1), replace=TRUE))
  ctl=data1[data1$Z==0,]
  trt=data1[data1$Z==1,]
  if(outcome=="binary"){
    var_trt=mean(trt$Y)*(1-mean(trt$Y))/nrow(trt)
    var_ctl=mean(ctl$Y)*(1-mean(ctl$Y))/nrow(ctl)
  }else{var_trt=var(trt$Y)/nrow(trt)
  var_ctl=var(ctl$Y)/nrow(ctl)}
  for(i in 1:BB){
    #Generate data
    ctl_boot=ctl[bootstrap_index0[,i],]
    trt_boot=trt[bootstrap_index1[,i],]
    data=rbind(ctl_boot,trt_boot)
    #data_AD=rbind(ctl_boot,trt)
    data_AD=data
    if(linear){
      if(PSgood==1) {pi.hat<-fitted(glm(Z~x+g,family=binomial,data=data))}
      else{pi.hat<-fitted(glm(Z~x,family=binomial,data=data))}
    }else{
      if(PSgood==1) {pi.hat<-fitted(glm(Z~x+G,family=binomial,data=data))}
      else{pi.hat<-fitted(glm(Z~G,family=binomial,data=data))}
    }
    data$w.hat<-pi.hat/(1-pi.hat)
    data$ps=ifelse(data$Z==1,1,data$w.hat)
    #IPW
    sr.att <-sum((data$Z-(1-data$Z)*data$w.hat)*data$Y)/sum(data$Z) 
    #Stabilized IPW normalizing weights
    sr_att=sum(data$ps*data$Z*data$Y)/sum(data$ps*data$Z)-sum(data$ps*(1-data$Z)*data$Y)/sum(data$ps*(1-data$Z))
    #Outcome regression (one unified model with treatment indicator)
    dat0 <- dat <- data.frame(Y=data$Y,Z=data$Z,x=data$x,g=data$g,G=data$G); dat0$Z<-0
    dat0_1=dat[dat$Z==1,];dat0_1$Z<-0
    dat0_0=dat[dat$Z==0,]
    
    #MAIC
    objfn <- function(a1, X){
      sum(exp(X %*% a1))
    }
    gradfn <- function(a1, X){
      colSums(sweep(X, 1, exp(X %*% a1), "*"))
    }
    mp1.mean=apply(data_AD[data_AD$Z==1,c("x","x.2","g","g.2","G","G.2")],2,mean)
    mp1.sd=apply(data_AD[data_AD$Z==1,c("x","x.2","g","g.2","G","G.2")],2,sd)
    if(linear){
      if(PSgood==1) {
        #center both mean and standard deviation
        X.EM.0_cen=data.matrix(data.frame(x_cen=data_AD[data_AD$Z==0,c("x")]-mp1.mean[[1]],
                                          x.2_cen=data_AD[data_AD$Z==0,c("x.2")]-((mp1.mean[[1]])^2+(mp1.sd[[1]])^2),
                                          g_cen=data_AD[data_AD$Z==0,c("g")]-mp1.mean[[3]],
                                          g.2_cen=data_AD[data_AD$Z==0,c("g.2")]-((mp1.mean[[3]])^2+(mp1.sd[[3]])^2)
        ))
        W2 <- weightit(Z ~ x+g, data = data_AD,method = "ebal", estimand = "ATT")
        opt1 <- optim(par = c(0,0,0,0), fn = objfn, gr = gradfn, X = X.EM.0_cen, method = "BFGS")
      }
      else{
        X.EM.0_cen=data.matrix(data.frame(
          x_cen=data_AD[data_AD$Z==0,c("x")]-mp1.mean[[1]],
          x.2_cen=data_AD[data_AD$Z==0,c("x.2")]-((mp1.mean[[1]])^2+(mp1.sd[[1]])^2)
        ))
        W2 <- weightit(Z ~ x, data = data_AD,method = "ebal", estimand = "ATT")
        opt1 <- optim(par = c(0,0), fn = objfn, gr = gradfn, X = X.EM.0_cen, method = "BFGS")
      }
    }else{
      if(PSgood==1) {
        #center both mean and standard deviation
        X.EM.0_cen=data.matrix(data.frame(x_cen=data_AD[data_AD$Z==0,c("x")]-mp1.mean[[1]],
                                          x.2_cen=data_AD[data_AD$Z==0,c("x.2")]-((mp1.mean[[1]])^2+(mp1.sd[[1]])^2),
                                          x_2_cen=data_AD[data_AD$Z==0,c("x.2")]-mp1.mean[[2]],
                                          x_2.2_cen=data_AD[data_AD$Z==0,c("x.4")]-((mp1.mean[[2]])^2+(mp1.sd[[2]])^2),
                                          g_cen=data_AD[data_AD$Z==0,c("G")]-mp1.mean[[5]],
                                          g.2_cen=data_AD[data_AD$Z==0,c("G.2")]-((mp1.mean[[5]])^2+(mp1.sd[[5]])^2)
        ))
        W2 <- weightit(Z ~ x+I(x^2)+G, data = data_AD,method = "ebal", estimand = "ATT")
        opt1 <- optim(par = c(0,0,0,0,0,0), fn = objfn, gr = gradfn, X = X.EM.0_cen, method = "BFGS")
      }
      else{
        X.EM.0_cen=data.matrix(data.frame(g_cen=data_AD[data_AD$Z==0,c("G")]-mp1.mean[[5]],
                                          g.2_cen=data_AD[data_AD$Z==0,c("G.2")]-((mp1.mean[[5]])^2+(mp1.sd[[5]])^2)
        ))
        W2 <- weightit(Z ~ G, data = data_AD,method = "ebal", estimand = "ATT")
        opt1 <- optim(par = c(0,0), fn = objfn, gr = gradfn, X = X.EM.0_cen, method = "BFGS")
      }
    }
    
    a1 <- opt1$par
    #weights of MAIC
    wt <- exp(X.EM.0_cen %*% a1)/sum(exp(X.EM.0_cen %*% a1))
    maic_att=mean(data_AD[data_AD$Z==1,"Y"])-sum(data_AD[data_AD$Z==0,"Y"]*wt)
    maic_est=sum(data_AD[data_AD$Z==0,"Y"]*wt)
    #balance weight-based methods
    wt_eb=W2$weights[data_AD$Z==0]/(sum(W2$weights[data_AD$Z==0]))
    eb_att=mean(data_AD[data_AD$Z==1,"Y"])-sum(wt_eb*data_AD[data_AD$Z==0,"Y"])
    eb_est=sum(wt_eb*data_AD[data_AD$Z==0,"Y"])
    
    if(linear){
      if(ORgood==1){
        if(outcome=="binary"){
          mu0 <- predict(glm(Y~Z*x+G,data=dat,family=binomial()),dat0,type="response")
          fit_all1=predict(glm(Y~Z*x+G,data=data,family=binomial()),newdata = dat0_1,type="response")
          fit_all0=predict(glm(Y~Z*x+G,data=data,family=binomial()),newdata = dat0_0,type="response")
          #STC
          fit=glm(Y~x+G,data=data_AD[data_AD$Z==0,],family=binomial())
          fit_park=glm(Y~x+G,data=data_AD[data_AD$Z==0,],family=binomial(),weights = wt_eb)
          # simulate covariates for the AgD study
          ipd_cor <- cor(data_AD[data_AD$Z==0, c("x", "G")])
          ipd_copula <-
            copula::normalCopula(copula::P2p(ipd_cor),
                                 dim = ncol(ipd_cor),
                                 dispstr = "un")
          myMvd_AgD <-
            mvdc(
              copula = ipd_copula,
              margins = c("norm","binom"),
              paramMargins = list(list(mean = mean(data_AD[data_AD$Z==1,"x"]), sd = sd(data_AD[data_AD$Z==1,"x"])),
                                  list(size = 1, prob = mean(data_AD[data_AD$Z==1,"G"])))
            )
          mycov_AgD <- data.frame(rMvdc(sum(data_AD$Z==1), myMvd_AgD))
          colnames(mycov_AgD)=c("x", "G")
          mycov_AgD$Z=1
        }else{
          #regression adjustment w/ one unified model
          mu0 <- predict(lm(Y~Z*x+G,data=dat),dat0)
          fit_all1=predict(lm(Y~Z*x+G,data=data),newdata = dat0_1)
          fit_all0=predict(lm(Y~Z*x+G,data=data),newdata = dat0_0)
          #STC
          fit=lm(Y~x+G,data=data_AD[data_AD$Z==0,])
          fit_park=lm(Y~x+G,data=data_AD[data_AD$Z==0,],weights = wt_eb)
          
          #regression adjustment using two models
          #mu=lm(Y~x+G,data=dat0_0)
        }
      }
      else{
        if(outcome=="binary"){
          mu0 <- predict(glm(Y~Z*x,data=dat,family=binomial()),dat0,type="response")
          fit_all1=predict(glm(Y~Z*x,data=data,family=binomial()),newdata = dat0_1,type="response")
          fit_all0=predict(glm(Y~Z*x,data=data,family=binomial()),newdata = dat0_0,type="response")
          #STC
          fit=glm(Y~x,data=data_AD[data_AD$Z==0,],family=binomial())
          fit_park=glm(Y~x,data=data_AD[data_AD$Z==0,],family=binomial(),weights = wt_eb)
          # simulate covariates for the AgD study
          mycov_AgD <- data.frame(Z=rep(1,sum(data_AD$Z==1)), x=rnorm(sum(data_AD$Z==1),mean = mean(data_AD[data_AD$Z==1,"x"]), sd = sd(data_AD[data_AD$Z==1,"x"])))
        }
        else{mu0 <- predict(lm(Y~Z*x,data=dat),dat0)
        fit_all1=predict(lm(Y~Z*x,data=data),newdata = dat0_1)
        fit_all0=predict(lm(Y~Z*x,data=data),newdata = dat0_0)
        fit=lm(Y~x,data=data_AD[data_AD$Z==0,])
        fit_park=lm(Y~x,data=data_AD[data_AD$Z==0,],weights = wt_eb)
        }
      }
    }else{
      if(ORgood==1){
        if(outcome=="binary"){
          mu0 <- predict(glm(Y~Z*x+I(x^2)+G,data=dat,family = binomial()),dat0,type="response")
          fit_all1=predict(glm(Y~Z*x+I(x^2)+G,data=data,family = binomial()),newdata = dat0_1,type="response")
          fit_all0=predict(glm(Y~Z*x+I(x^2)+G,data=data,family = binomial()),newdata = dat0_0,type="response")
          #STC
          fit=glm(Y~x+I(x^2)+G,data=data_AD[data_AD$Z==0,],family = binomial())
          fit_park=glm(Y~x+I(x^2)+G,data=data_AD[data_AD$Z==0,],family = binomial(),weights = wt_eb)
          # simulate covariates for the AgD study
          ipd_cor <- cor(data_AD[data_AD$Z==0, c("x", "x.2","G")])
          ipd_copula <-
            copula::normalCopula(copula::P2p(ipd_cor),
                                 dim = ncol(ipd_cor),
                                 dispstr = "un")
          myMvd_AgD <-
            mvdc(
              copula = ipd_copula,
              margins = c("norm","norm","binom"),
              paramMargins = list(list(mean = mean(data_AD[data_AD$Z==1,"x"]), sd = sd(data_AD[data_AD$Z==1,"x"])),
                                  list(mean = mean(data_AD[data_AD$Z==1,"x.2"]), sd = sd(data_AD[data_AD$Z==1,"x.2"])),
                                  list(size = 1, prob = mean(data_AD[data_AD$Z==1,"G"])))
            )
          mycov_AgD <- data.frame(rMvdc(sum(data_AD$Z==1), myMvd_AgD))
          colnames(mycov_AgD)=c("x","x.2","G")
          mycov_AgD$Z=1
        }
        else{mu0 <- predict(lm(Y~Z*x+I(x^2)+G,data=dat),dat0)
        fit_all1=predict(lm(Y~Z*x+I(x^2)+G,data=data),newdata = dat0_1)
        fit_all0=predict(lm(Y~Z*x+I(x^2)+G,data=data),newdata = dat0_0)
        #STC
        fit=lm(Y~x+I(x^2)+G,data=data_AD[data_AD$Z==0,])
        fit_park=lm(Y~x+I(x^2)+G,data=data_AD[data_AD$Z==0,],weights = wt_eb)
        }
      }
      else{
        if(outcome=="binary"){
          mu0 <- predict(glm(Y~Z+G,data=dat,family = binomial()),dat0,type="response")
          fit_all1=predict(glm(Y~Z+G,data=data,family = binomial()),newdata = dat0_1,type="response")
          fit_all0=predict(glm(Y~Z+G,data=data,family = binomial()),newdata = dat0_0,type="response")
          #STC
          fit=glm(Y~G,data=data_AD[data_AD$Z==0,],family = binomial())
          fit_park=glm(Y~G,data=data_AD[data_AD$Z==0,],family = binomial(),weights = wt_eb)
          # simulate covariates for the AgD study
          mycov_AgD <- data.frame(Z=rep(1,sum(data_AD$Z==1)), G=rbinom(n=sum(data_AD$Z==1),size=1,prob = mean(data_AD[data_AD$Z==1,"G"])))
        }
        else{mu0 <- predict(lm(Y~Z+G,data=dat),dat0)
        fit_all1=predict(lm(Y~Z+G,data=data),newdata = dat0_1)
        fit_all0=predict(lm(Y~Z+G,data=data),newdata = dat0_0)
        fit=lm(Y~G,data=data_AD[data_AD$Z==0,])
        fit_park=lm(Y~G,data=data_AD[data_AD$Z==0,],weights = wt_eb)
        }
      }
    }
    #OR.att=mean(data[data$Z==1,"Y"])-mean(mu0[Z==1])
    #OR.att=mean(data[data$Z==1,"Y"])-mean(mu0[dat$Z==1])
    OR.att=mean(data[data$Z==1,"Y"])-mean(fit_all1)
    OR_est=mean(fit_all1)

    #Stabilized AIPW (one unified model)
    dr.att.norm=(mean(data[data$Z==1,"Y"])-mean(fit_all1))-sum(data[data$Z==0,"ps"]*(data[data$Z==0,"Y"]-fit_all0))/sum(data$ps[data$Z==0])

    if(outcome=="binary"){
      mu0_stc <- predict(fit,newdata=mycov_AgD,type="response")
      stc_pred=mean(mu0_stc)
      mu0_stc_park <- predict(fit_park,newdata=mycov_AgD,type="response")
      stc_pred_park=mean(mu0_stc_park)
      #unweighted AIPW using separate regression model
      aipw_att=(mean(data[data$Z==1,"Y"])-mean(predict(fit,newdata=data[data$Z==1,],type="response")))-sum(data$ps[data$Z==0]*(data[data$Z==0,"Y"]-predict(fit,type="response")))/sum(data$Z)
      #AIPW using separate regression model
      aipw_norm_att=(mean(data[data$Z==1,"Y"])-mean(predict(fit,newdata=data[data$Z==1,],type="response")))-sum(data$ps[data$Z==0]*(data[data$Z==0,"Y"]-predict(fit,type="response")))/sum(data$ps[data$Z==0])
      #double robust estimator for aggregate data
      dr_att=(stc_pred)+sum(wt_eb*(data_AD[data_AD$Z==0,"Y"]-predict(fit,type="response")))
      dr_att2=(mean(data_AD[data_AD$Z==1,"Y"])-stc_pred)-sum(wt_eb*(data_AD[data_AD$Z==0,"Y"]-predict(fit,type="response")))
      dr_att_half=sum(wt_eb*(data_AD[data_AD$Z==0,"Y"]-predict(fit,type="response")))
    }else{
      #mu0_ra=predict(mu,newdata=dat[dat$Z==1,])
      #ra_pred=mean(mu0_ra)
      mu0_stc <- predict(fit,newdata=data_AD[data_AD$Z==1,])
      stc_pred=mean(mu0_stc)
      mu0_stc_park <- predict(fit_park,newdata=data_AD[data_AD$Z==1,])
      stc_pred_park=mean(mu0_stc_park)
      #unweighted AIPW using separate regression model
      aipw_att=(mean(data[data$Z==1,"Y"])-mean(mu0_stc))-sum(data$ps[data$Z==0]*(data[data$Z==0,"Y"]-predict(fit)))/sum(data$Z)
      #AIPW using separate regression model
      aipw_norm_att=(mean(data[data$Z==1,"Y"])-mean(mu0_stc))-sum(data$ps[data$Z==0]*(data[data$Z==0,"Y"]-predict(fit)))/sum(data$ps[data$Z==0])
      #double robust estimator for aggregate data
      #dr_att=(stc_pred)+sum(wt_eb*(data_AD[data_AD$Z==0,"Y"]-predict(fit)))
      dr_att2=(mean(data_AD[data_AD$Z==1,"Y"])-stc_pred)-sum(wt_eb*(data_AD[data_AD$Z==0,"Y"]-predict(fit)))
      #dr_att_half=sum(wt_eb*(data_AD[data_AD$Z==0,"Y"]-predict(fit)))
    }
    #ra_att=mean(dat[dat$Z==1,"Y"])-ra_pred
    stc_att=mean(data_AD[data_AD$Z==1,"Y"])-stc_pred
    stc_att_park=mean(data_AD[data_AD$Z==1,"Y"])-stc_pred_park
    #double robust estimator for aggregate data
    summary.vec=c(sr.att,sr_att,eb_att,maic_att,OR.att,stc_att,dr.att.norm,dr_att2,stc_att_park)
    out=rbind(out,summary.vec)
  }
  cp=ifelse((att>=apply(out,2,quantile,probs=0.025,na.rm = TRUE))&(att<=apply(out,2,quantile,probs=0.975, na.rm = TRUE)),1,0)
  ci.length=apply(out,2,quantile,probs=0.975,na.rm = TRUE)-apply(out,2,quantile,probs=0.025,na.rm = TRUE)
  #print(att);print(ci.length);print(att);#return(out)
  print(seed)
  return(c(apply(out,2,sd),cp,ci.length))
}
bootSE1_1_1=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=1,PSgood=1,n=3000,n_A=1000,n_B=1000,linear=TRUE,outcome="nbinary"))
bootSE0_1_1=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=0,PSgood=1,n=3000,n_A=1000,n_B=1000,linear=TRUE,outcome="nbinary"))
bootSE1_0_1=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=1,PSgood=0,n=3000,n_A=1000,n_B=1000,linear=TRUE,outcome="nbinary"))
bootSE0_0_1=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=0,PSgood=0,n=3000,n_A=1000,n_B=1000,linear=TRUE,outcome="nbinary"))
bootSE1_1_2=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=1,PSgood=1,n=3000,n_A=1000,n_B=1000,linear=TRUE,outcome="binary"))
bootSE0_1_2=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=0,PSgood=1,n=3000,n_A=1000,n_B=1000,linear=TRUE,outcome="binary"))
bootSE1_0_2=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=1,PSgood=0,n=3000,n_A=1000,n_B=1000,linear=TRUE,outcome="binary"))
bootSE0_0_2=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=0,PSgood=0,n=3000,n_A=1000,n_B=1000,linear=TRUE,outcome="binary"))
bootSE1_1_3=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=1,PSgood=1,n=3000,n_A=1000,n_B=1000,linear=FALSE,outcome="nbinary"))
bootSE0_1_3=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=0,PSgood=1,n=3000,n_A=1000,n_B=1000,linear=FALSE,outcome="nbinary"))
bootSE1_0_3=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=1,PSgood=0,n=3000,n_A=1000,n_B=1000,linear=FALSE,outcome="nbinary"))
bootSE0_0_3=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=0,PSgood=0,n=3000,n_A=1000,n_B=1000,linear=FALSE,outcome="nbinary"))
bootSE1_1_4=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=1,PSgood=1,n=3000,n_A=1000,n_B=1000,linear=FALSE,outcome="binary"))
bootSE0_1_4=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=0,PSgood=1,n=3000,n_A=1000,n_B=1000,linear=FALSE,outcome="binary"))
bootSE1_0_4=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=1,PSgood=0,n=3000,n_A=1000,n_B=1000,linear=FALSE,outcome="binary"))
bootSE0_0_4=sapply(1:500, function(i) sim2.BB(BB=100,seed=i,ORgood=0,PSgood=0,n=3000,n_A=1000,n_B=1000,linear=FALSE,outcome="binary"))




df3_4_iv3=data.frame(rbind(apply(result1_4, 1, mean),apply(result2_4, 1, mean),apply(result3_4, 1, mean),apply(result4_4, 1, mean)))
df3_4_iv3_sd=data.frame(rbind(apply(result1_4, 1, sd),apply(result2_4, 1, sd),apply(result3_4, 1, sd),apply(result4_4, 1, sd)))
df3_4_iv3$iv_strength=1
#df3_4_iv3$N="n=1000 per study"
#df3_4_iv3$outcome="Binary (nonlinear)"
df3_4_iv3$Misspecification=c("Neither","OR","PS","Both")
colnames(df3_4_iv3)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","iv_strength")
df3_4_iv3_sd$Misspecification=c("Neither","OR","PS","Both")
colnames(df3_4_iv3_sd)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","Misspecification")
rownames(df3_4_iv3)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")

#bias+CI width+coverage probability
c(df3_4[,c(1,2,3,4,5,6,7,8,9)][1,1],round(apply(bootSE1_1_4, 1, mean),3)[c(19,10)],df3_4[,c(1,2,3,4,5,6,7,8,9)][3,1],round(apply(bootSE1_0_4, 1, mean),3)[c(19,10)])#IPW1
c(df3_4[,c(1,2,3,4,5,6,7,8,9)][1,2],round(apply(bootSE1_1_4, 1, mean),3)[c(20,11)],df3_4[,c(1,2,3,4,5,6,7,8,9)][3,2],round(apply(bootSE1_0_4, 1, mean),3)[c(20,11)])#IPW2
c(df3_4[,c(1,2,3,4,5,6,7,8,9)][1,3],round(apply(bootSE1_1_4, 1, mean),3)[c(21,12)],df3_4[,c(1,2,3,4,5,6,7,8,9)][3,3],round(apply(bootSE1_0_4, 1, mean),3)[c(21,12)])#MAIC1
c(df3_4[,c(1,2,3,4,5,6,7,8,9)][1,4],round(apply(bootSE1_1_4, 1, mean),3)[c(22,13)],df3_4[,c(1,2,3,4,5,6,7,8,9)][3,4],round(apply(bootSE1_0_4, 1, mean),3)[c(22,13)])#MAIC2
#RA
c(df3_4[,c(1,2,3,4,5,6,7,8,9)][2,5],round(apply(bootSE0_1_4, 1, mean),3)[c(23,14)],df3_4[,c(1,2,3,4,5,6,7,8,9)][1,5],round(apply(bootSE1_1_4, 1, mean),3)[c(23,14)])#RA
c(df3_4[,c(1,2,3,4,5,6,7,8,9)][2,6],round(apply(bootSE0_1_4, 1, mean),3)[c(24,15)],df3_4[,c(1,2,3,4,5,6,7,8,9)][1,6],round(apply(bootSE1_1_4, 1, mean),3)[c(24,15)])#STC
c(df3_4[,c(1,2,3,4,5,6,7,8,9)][1,7],round(apply(bootSE1_1_4, 1, mean),3)[c(25,16)],df3_4[,c(1,2,3,4,5,6,7,8,9)][2,7],round(apply(bootSE0_1_4, 1, mean),3)[c(25,16)],df3_4[,c(1,2,3,4,5,6,7,8,9)][3,7],round(apply(bootSE1_0_4, 1, mean),3)[c(25,16)],df3_4[,c(1,2,3,4,5,6,7,8,9)][4,7],round(apply(bootSE0_0_4, 1, mean),3)[c(25,16)])#AIPW
c(df3_4[,c(1,2,3,4,5,6,7,8,9)][1,8],round(apply(bootSE1_1_4, 1, mean),3)[c(26,17)],df3_4[,c(1,2,3,4,5,6,7,8,9)][2,8],round(apply(bootSE0_1_4, 1, mean),3)[c(26,17)],df3_4[,c(1,2,3,4,5,6,7,8,9)][3,8],round(apply(bootSE1_0_4, 1, mean),3)[c(26,17)],df3_4[,c(1,2,3,4,5,6,7,8,9)][4,8],round(apply(bootSE0_0_4, 1, mean),3)[c(26,17)])#DR
c(df3_4[,c(1,2,3,4,5,6,7,8,9)][1,9],round(apply(bootSE1_1_4, 1, mean),3)[c(27,18)],df3_4[,c(1,2,3,4,5,6,7,8,9)][2,9],round(apply(bootSE0_1_4, 1, mean),3)[c(27,18)],df3_4[,c(1,2,3,4,5,6,7,8,9)][3,9],round(apply(bootSE1_0_4, 1, mean),3)[c(27,18)],df3_4[,c(1,2,3,4,5,6,7,8,9)][4,9],round(apply(bootSE0_0_4, 1, mean),3)[c(27,18)])#DR (Park)

#Coverage plot
cv_index=10:18
cv_index=19:27
cv3_1=data.frame(rbind(apply(bootSE1_1_1[cv_index,], 1, mean),apply(bootSE0_1_1[cv_index,], 1, mean),apply(bootSE1_0_1[cv_index,], 1, mean),apply(bootSE0_0_1[cv_index,], 1, mean)))
cv3_2=data.frame(rbind(apply(bootSE1_1_2[cv_index,], 1, mean),apply(bootSE0_1_2[cv_index,], 1, mean),apply(bootSE1_0_2[cv_index,], 1, mean),apply(bootSE0_0_2[cv_index,], 1, mean)))
cv3_3=data.frame(rbind(apply(bootSE1_1_3[cv_index,], 1, mean),apply(bootSE0_1_3[cv_index,], 1, mean),apply(bootSE1_0_3[cv_index,], 1, mean),apply(bootSE0_0_3[cv_index,], 1, mean)))
cv3_4=data.frame(rbind(apply(bootSE1_1_4[cv_index,], 1, mean),apply(bootSE0_1_4[cv_index,], 1, mean),apply(bootSE1_0_4[cv_index,], 1, mean),apply(bootSE0_0_4[cv_index,], 1, mean)))
colnames(cv3_1)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)")
colnames(cv3_2)=colnames(cv3_3)=colnames(cv3_4)=colnames(cv3_1)
cv3_1$Misspecification=cv3_2$Misspecification=cv3_3$Misspecification=cv3_4$Misspecification=c("Neither","OR","PS","Both")
rownames(cv3_1)=rownames(cv3_2)=rownames(cv3_3)=rownames(cv3_4)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
cv3_1$N=cv3_2$N=cv3_3$N=cv3_4$N="n=1000 per study"
#n=500
cv2_1=data.frame(rbind(apply(bootSE1_1_1_na500_nb500[cv_index,], 1, mean),apply(bootSE0_1_1_na500_nb500[cv_index,], 1, mean),apply(bootSE1_0_1_na500_nb500[cv_index,], 1, mean),apply(bootSE0_0_1_na500_nb500[cv_index,], 1, mean)))
cv2_2=data.frame(rbind(apply(bootSE1_1_2_na500_nb500[cv_index,], 1, mean),apply(bootSE0_1_2_na500_nb500[cv_index,], 1, mean),apply(bootSE1_0_2_na500_nb500[cv_index,], 1, mean),apply(bootSE0_0_2_na500_nb500[cv_index,], 1, mean)))
cv2_3=data.frame(rbind(apply(bootSE1_1_3_na500_nb500[cv_index,], 1, mean),apply(bootSE0_1_3_na500_nb500[cv_index,], 1, mean),apply(bootSE1_0_3_na500_nb500[cv_index,], 1, mean),apply(bootSE0_0_3_na500_nb500[cv_index,], 1, mean)))
cv2_4=data.frame(rbind(apply(bootSE1_1_4_na500_nb500[cv_index,], 1, mean),apply(bootSE0_1_4_na500_nb500[cv_index,], 1, mean),apply(bootSE1_0_4_na500_nb500[cv_index,], 1, mean),apply(bootSE0_0_4_na500_nb500[cv_index,], 1, mean)))
colnames(cv2_1)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)")
colnames(cv2_2)=colnames(cv2_3)=colnames(cv2_4)=colnames(cv2_1)
cv2_1$Misspecification=cv2_2$Misspecification=cv2_3$Misspecification=cv2_4$Misspecification=c("Neither","OR","PS","Both")
rownames(cv2_1)=rownames(cv2_2)=rownames(cv2_3)=rownames(cv2_4)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
cv2_1$N=cv2_2$N=cv2_3$N=cv2_4$N="n=500 per study"
#N=200
cv1_1=data.frame(rbind(apply(bootSE1_1_1_na200_nb200[cv_index,], 1, mean),apply(bootSE0_1_1_na200_nb200[cv_index,], 1, mean),apply(bootSE1_0_1_na200_nb200[cv_index,], 1, mean),apply(bootSE0_0_1_na200_nb200[cv_index,], 1, mean)))
cv1_2=data.frame(rbind(apply(bootSE1_1_2_na200_nb200[cv_index,], 1, mean),apply(bootSE0_1_2_na200_nb200[cv_index,], 1, mean),apply(bootSE1_0_2_na200_nb200[cv_index,], 1, mean),apply(bootSE0_0_2_na200_nb200[cv_index,], 1, mean)))
cv1_3=data.frame(rbind(apply(bootSE1_1_3_na200_nb200[cv_index,], 1, mean),apply(bootSE0_1_3_na200_nb200[cv_index,], 1, mean),apply(bootSE1_0_3_na200_nb200[cv_index,], 1, mean),apply(bootSE0_0_3_na200_nb200[cv_index,], 1, mean)))
cv1_4=data.frame(rbind(apply(bootSE1_1_4_na200_nb200[cv_index,], 1, mean),apply(bootSE0_1_4_na200_nb200[cv_index,], 1, mean),apply(bootSE1_0_4_na200_nb200[cv_index,], 1, mean),apply(bootSE0_0_4_na200_nb200[cv_index,], 1, mean)))
colnames(cv1_1)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)")
colnames(cv1_2)=colnames(cv1_3)=colnames(cv1_4)=colnames(cv1_1)
cv1_1$Misspecification=cv1_2$Misspecification=cv1_3$Misspecification=cv1_4$Misspecification=c("Neither","OR","PS","Both")
rownames(cv1_1)=rownames(cv1_2)=rownames(cv1_3)=rownames(cv1_4)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
cv1_1$N=cv1_2$N=cv1_3$N=cv1_4$N="n=200 per study"
cv_1=rbind(cv1_1,cv2_1,cv3_1)
cv_1$outcome="Linear"
cv_2=rbind(cv1_2,cv2_2,cv3_2)
cv_2$outcome="Binary (linear)"
cv_3=rbind(cv1_3,cv2_3,cv3_3)
cv_3$outcome="Nonlinear"
cv_4=rbind(cv1_4,cv2_4,cv3_4)
cv_4$outcome="Binary (nonlinear)"
cv=rbind(cv_1,cv_2,cv_3,cv_4)

# Melt to long format
cv_long <- melt(cv,
                id.vars = c("Misspecification", "N", "outcome"),
                variable.name = "Method",
                value.name = "Coverage")#ciLength
cv_long$N=factor(cv_long$N,levels = c("n=200 per study","n=500 per study","n=1000 per study"))
filtered_cv_long <- cv_long %>%
  filter(!(Method %in% c("IPW1", "IPW2", "MAIC1", "MAIC2", "RA", "STC") &
             Misspecification %in% c("Neither", "Both")))
ggplot(filtered_cv_long, aes(x = Method, y = Coverage, fill = Misspecification)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
  #geom_bar(stat = "identity", position = "dodge") +
  #geom_point(position = position_dodge(width=0.9),size=2)+
  #geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
  #              width = 0.2, position = position_dodge(width = 0.9)) +
  facet_grid(outcome ~ N,scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  theme_minimal(base_size = 13) +
  labs(x = "Method", y = "C.I. Length") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))

cv5_1=data.frame(rbind(apply(bootSE1_1_1_na500_nb200[cv_index,], 1, mean),apply(bootSE0_1_1_na500_nb200[cv_index,], 1, mean),apply(bootSE1_0_1_na500_nb200[cv_index,], 1, mean),apply(bootSE0_0_1_na500_nb200[cv_index,], 1, mean)))
cv5_2=data.frame(rbind(apply(bootSE1_1_2_na500_nb200[cv_index,], 1, mean),apply(bootSE0_1_2_na500_nb200[cv_index,], 1, mean),apply(bootSE1_0_2_na500_nb200[cv_index,], 1, mean),apply(bootSE0_0_2_na500_nb200[cv_index,], 1, mean)))
cv5_3=data.frame(rbind(apply(bootSE1_1_3_na500_nb200[cv_index,], 1, mean),apply(bootSE0_1_3_na500_nb200[cv_index,], 1, mean),apply(bootSE1_0_3_na500_nb200[cv_index,], 1, mean),apply(bootSE0_0_3_na500_nb200[cv_index,], 1, mean)))
cv5_4=data.frame(rbind(apply(bootSE1_1_4_na500_nb200[cv_index,], 1, mean),apply(bootSE0_1_4_na500_nb200[cv_index,], 1, mean),apply(bootSE1_0_4_na500_nb200[cv_index,], 1, mean),apply(bootSE0_0_4_na500_nb200[cv_index,], 1, mean)))
colnames(cv5_1)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)")
colnames(cv5_2)=colnames(cv5_3)=colnames(cv5_4)=colnames(cv5_1)
cv5_1$Misspecification=cv5_2$Misspecification=cv5_3$Misspecification=cv5_4$Misspecification=c("Neither","OR","PS","Both")
rownames(cv5_1)=rownames(cv5_2)=rownames(cv5_3)=rownames(cv5_4)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
cv5_1$N=cv5_2$N=cv5_3$N=cv5_4$N="n_A=500;n_B=200"

cv4_1=data.frame(rbind(apply(bootSE1_1_1_na200_nb500[cv_index,], 1, mean),apply(bootSE0_1_1_na200_nb500[cv_index,], 1, mean),apply(bootSE1_0_1_na200_nb500[cv_index,], 1, mean),apply(bootSE0_0_1_na200_nb500[cv_index,], 1, mean)))
cv4_2=data.frame(rbind(apply(bootSE1_1_2_na200_nb500[cv_index,], 1, mean),apply(bootSE0_1_2_na200_nb500[cv_index,], 1, mean),apply(bootSE1_0_2_na200_nb500[cv_index,], 1, mean),apply(bootSE0_0_2_na200_nb500[cv_index,], 1, mean)))
cv4_3=data.frame(rbind(apply(bootSE1_1_3_na200_nb500[cv_index,], 1, mean),apply(bootSE0_1_3_na200_nb500[cv_index,], 1, mean),apply(bootSE1_0_3_na200_nb500[cv_index,], 1, mean),apply(bootSE0_0_3_na200_nb500[cv_index,], 1, mean)))
cv4_4=data.frame(rbind(apply(bootSE1_1_4_na200_nb500[cv_index,], 1, mean),apply(bootSE0_1_4_na200_nb500[cv_index,], 1, mean),apply(bootSE1_0_4_na200_nb500[cv_index,], 1, mean),apply(bootSE0_0_4_na200_nb500[cv_index,], 1, mean)))
colnames(cv4_1)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)")
colnames(cv4_2)=colnames(cv4_3)=colnames(cv4_4)=colnames(cv4_1)
cv4_1$Misspecification=cv4_2$Misspecification=cv4_3$Misspecification=cv4_4$Misspecification=c("Neither","OR","PS","Both")
rownames(cv4_1)=rownames(cv4_2)=rownames(cv4_3)=rownames(cv4_4)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
cv4_1$N=cv4_2$N=cv4_3$N=cv4_4$N="n_A=200;n_B=500"

cv3_1=data.frame(rbind(apply(bootSE1_1_1_na1000_nb200[cv_index,], 1, mean),apply(bootSE0_1_1_na1000_nb200[cv_index,], 1, mean),apply(bootSE1_0_1_na1000_nb200[cv_index,], 1, mean),apply(bootSE0_0_1_na1000_nb200[cv_index,], 1, mean)))
cv3_2=data.frame(rbind(apply(bootSE1_1_2_na1000_nb200[cv_index,], 1, mean),apply(bootSE0_1_2_na1000_nb200[cv_index,], 1, mean),apply(bootSE1_0_2_na1000_nb200[cv_index,], 1, mean),apply(bootSE0_0_2_na1000_nb200[cv_index,], 1, mean)))
cv3_3=data.frame(rbind(apply(bootSE1_1_3_na1000_nb200[cv_index,], 1, mean),apply(bootSE0_1_3_na1000_nb200[cv_index,], 1, mean),apply(bootSE1_0_3_na1000_nb200[cv_index,], 1, mean),apply(bootSE0_0_3_na1000_nb200[cv_index,], 1, mean)))
cv3_4=data.frame(rbind(apply(bootSE1_1_4_na1000_nb200[cv_index,], 1, mean),apply(bootSE0_1_4_na1000_nb200[cv_index,], 1, mean),apply(bootSE1_0_4_na1000_nb200[cv_index,], 1, mean),apply(bootSE0_0_4_na1000_nb200[cv_index,], 1, mean)))
colnames(cv3_1)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)")
colnames(cv3_2)=colnames(cv3_3)=colnames(cv3_4)=colnames(cv3_1)
cv3_1$Misspecification=cv3_2$Misspecification=cv3_3$Misspecification=cv3_4$Misspecification=c("Neither","OR","PS","Both")
rownames(cv3_1)=rownames(cv3_2)=rownames(cv3_3)=rownames(cv3_4)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
cv3_1$N=cv3_2$N=cv3_3$N=cv3_4$N="n_A=1000;n_B=200"
#n=500
cv2_1=data.frame(rbind(apply(bootSE1_1_1_na1000_nb500[cv_index,], 1, mean),apply(bootSE0_1_1_na1000_nb500[cv_index,], 1, mean),apply(bootSE1_0_1_na1000_nb500[cv_index,], 1, mean),apply(bootSE0_0_1_na1000_nb500[cv_index,], 1, mean)))
cv2_2=data.frame(rbind(apply(bootSE1_1_2_na1000_nb500[cv_index,], 1, mean),apply(bootSE0_1_2_na1000_nb500[cv_index,], 1, mean),apply(bootSE1_0_2_na1000_nb500[cv_index,], 1, mean),apply(bootSE0_0_2_na1000_nb500[cv_index,], 1, mean)))
cv2_3=data.frame(rbind(apply(bootSE1_1_3_na1000_nb500[cv_index,], 1, mean),apply(bootSE0_1_3_na1000_nb500[cv_index,], 1, mean),apply(bootSE1_0_3_na1000_nb500[cv_index,], 1, mean),apply(bootSE0_0_3_na1000_nb500[cv_index,], 1, mean)))
cv2_4=data.frame(rbind(apply(bootSE1_1_4_na1000_nb500[cv_index,], 1, mean),apply(bootSE0_1_4_na1000_nb500[cv_index,], 1, mean),apply(bootSE1_0_4_na1000_nb500[cv_index,], 1, mean),apply(bootSE0_0_4_na1000_nb500[cv_index,], 1, mean)))
colnames(cv2_1)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)")
colnames(cv2_2)=colnames(cv2_3)=colnames(cv2_4)=colnames(cv2_1)
cv2_1$Misspecification=cv2_2$Misspecification=cv2_3$Misspecification=cv2_4$Misspecification=c("Neither","OR","PS","Both")
rownames(cv2_1)=rownames(cv2_2)=rownames(cv2_3)=rownames(cv2_4)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
cv2_1$N=cv2_2$N=cv2_3$N=cv2_4$N="n_A=1000;n_B=500"
#N=200
cv1_1=data.frame(rbind(apply(bootSE1_1_1_na200_nb1000[cv_index,], 1, mean),apply(bootSE0_1_1_na200_nb1000[cv_index,], 1, mean),apply(bootSE1_0_1_na200_nb1000[cv_index,], 1, mean),apply(bootSE0_0_1_na200_nb1000[cv_index,], 1, mean)))
cv1_2=data.frame(rbind(apply(bootSE1_1_2_na200_nb1000[cv_index,], 1, mean),apply(bootSE0_1_2_na200_nb1000[cv_index,], 1, mean),apply(bootSE1_0_2_na200_nb1000[cv_index,], 1, mean),apply(bootSE0_0_2_na200_nb1000[cv_index,], 1, mean)))
cv1_3=data.frame(rbind(apply(bootSE1_1_3_na200_nb1000[cv_index,], 1, mean),apply(bootSE0_1_3_na200_nb1000[cv_index,], 1, mean),apply(bootSE1_0_3_na200_nb1000[cv_index,], 1, mean),apply(bootSE0_0_3_na200_nb1000[cv_index,], 1, mean)))
cv1_4=data.frame(rbind(apply(bootSE1_1_4_na200_nb1000[cv_index,], 1, mean),apply(bootSE0_1_4_na200_nb1000[cv_index,], 1, mean),apply(bootSE1_0_4_na200_nb1000[cv_index,], 1, mean),apply(bootSE0_0_4_na200_nb1000[cv_index,], 1, mean)))
colnames(cv1_1)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)")
colnames(cv1_2)=colnames(cv1_3)=colnames(cv1_4)=colnames(cv1_1)
cv1_1$Misspecification=cv1_2$Misspecification=cv1_3$Misspecification=cv1_4$Misspecification=c("Neither","OR","PS","Both")
rownames(cv1_1)=rownames(cv1_2)=rownames(cv1_3)=rownames(cv1_4)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
cv1_1$N=cv1_2$N=cv1_3$N=cv1_4$N="n_A=200;n_B=1000"

cv0_1=data.frame(rbind(apply(bootSE1_1_1_na500_nb1000[cv_index,], 1, mean),apply(bootSE0_1_1_na500_nb1000[cv_index,], 1, mean),apply(bootSE1_0_1_na500_nb1000[cv_index,], 1, mean),apply(bootSE0_0_1_na500_nb1000[cv_index,], 1, mean)))
cv0_2=data.frame(rbind(apply(bootSE1_1_2_na500_nb1000[cv_index,], 1, mean),apply(bootSE0_1_2_na500_nb1000[cv_index,], 1, mean),apply(bootSE1_0_2_na500_nb1000[cv_index,], 1, mean),apply(bootSE0_0_2_na500_nb1000[cv_index,], 1, mean)))
cv0_3=data.frame(rbind(apply(bootSE1_1_3_na500_nb1000[cv_index,], 1, mean),apply(bootSE0_1_3_na500_nb1000[cv_index,], 1, mean),apply(bootSE1_0_3_na500_nb1000[cv_index,], 1, mean),apply(bootSE0_0_3_na500_nb1000[cv_index,], 1, mean)))
cv0_4=data.frame(rbind(apply(bootSE1_1_4_na500_nb1000[cv_index,], 1, mean),apply(bootSE0_1_4_na500_nb1000[cv_index,], 1, mean),apply(bootSE1_0_4_na500_nb1000[cv_index,], 1, mean),apply(bootSE0_0_4_na500_nb1000[cv_index,], 1, mean)))
colnames(cv0_1)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)")
colnames(cv0_2)=colnames(cv0_3)=colnames(cv0_4)=colnames(cv0_1)
cv0_1$Misspecification=cv0_2$Misspecification=cv0_3$Misspecification=cv0_4$Misspecification=c("Neither","OR","PS","Both")
rownames(cv0_1)=rownames(cv0_2)=rownames(cv0_3)=rownames(cv0_4)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
cv0_1$N=cv0_2$N=cv0_3$N=cv0_4$N="n_A=500;n_B=1000"
cv_1=rbind(cv0_1,cv1_1,cv2_1,cv3_1,cv4_1,cv5_1)
cv_1$outcome="Linear"
cv_2=rbind(cv0_2,cv1_2,cv2_2,cv3_2,cv4_2,cv5_2)
cv_2$outcome="Binary (linear)"
cv_3=rbind(cv0_3,cv1_3,cv2_3,cv3_3,cv4_3,cv5_3)
cv_3$outcome="Nonlinear"
cv_4=rbind(cv0_4,cv1_4,cv2_4,cv3_4,cv4_4,cv5_4)
cv_4$outcome="Binary (nonlinear)"
cv=rbind(cv_1,cv_2,cv_3,cv_4)
# Melt to long format
cv_long <- melt(cv,
                id.vars = c("Misspecification", "N", "outcome"),
                variable.name = "Method",
                value.name = "ciLength")#Coverage

cv_long$N=factor(cv_long$N,levels = c("n_A=200;n_B=500","n_A=500;n_B=200",
                                      "n_A=200;n_B=1000","n_A=1000;n_B=200",
                                      "n_A=500;n_B=1000","n_A=1000;n_B=500"))
filtered_cv_long <- cv_long %>%
  filter(!(Method %in% c("IPW1", "IPW2", "MAIC1", "MAIC2", "RA", "STC") &
             Misspecification %in% c("Neither", "Both")))
ggplot(filtered_cv_long, aes(x = Method, y = ciLength, fill = Misspecification)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  #geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
  #geom_bar(stat = "identity", position = "dodge") +
  #geom_point(position = position_dodge(width=0.9),size=2)+
  #geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
  #              width = 0.2, position = position_dodge(width = 0.9)) +
  facet_grid(outcome ~ N,scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  theme_minimal(base_size = 13) +
  labs(x = "Method", y = "C.I. Length") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))

#covariate shift+instrumental strength
seeds <- 1:1000
count <- 0
result1_2_x1LS=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=1,PSgood=1,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=3,shift_param1=2)})
seeds <- 1:1000
count <- 0
result2_2_x1LS=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=0,PSgood=1,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=3,shift_param1=2)})
seeds <- 1:1000
count <- 0
result3_2_x1LS=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=1,PSgood=0,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=3,shift_param1=2)})
seeds <- 1:1000
count <- 0
result4_2_x1LS=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=0,PSgood=0,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=3,shift_param1=2)})

x1LS=data.frame(rbind(apply(result1_2_x1LS, 1, mean),apply(result2_2_x1LS, 1, mean),apply(result3_2_x1LS, 1, mean),apply(result4_2_x1LS, 1, mean)))
x1LS_sd=data.frame(rbind(apply(result1_2_x1LS, 1, sd),apply(result2_2_x1LS, 1, sd),apply(result3_2_x1LS, 1, sd),apply(result4_2_x1LS, 1, sd)))
x1LS$Misspecification=c("Neither","OR","PS","Both")
colnames(x1LS)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","Misspecification")
x1LS_sd$Misspecification=c("Neither","OR","PS","Both")
colnames(x1LS_sd)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","Misspecification")
rownames(x1LS)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
x1LS$shift="Normal(3,2)"

seeds <- 1:1000
count <- 0
result1_2_x1LS1=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=1,PSgood=1,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=3,shift_param1=3)})
seeds <- 1:1000
count <- 0
result2_2_x1LS1=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=0,PSgood=1,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=3,shift_param1=3)})
seeds <- 1:1000
count <- 0
result3_2_x1LS1=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=1,PSgood=0,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=3,shift_param1=3)})
seeds <- 1:1000
count <- 0
result4_2_x1LS1=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=0,PSgood=0,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=3,shift_param1=3)})



seeds <- 1:1000
count <- 0
result1_2_x1RS=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=1,PSgood=1,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=10,shift_param1=2)})
seeds <- 1:1000
count <- 0
result2_2_x1RS=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=0,PSgood=1,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=10,shift_param1=2)})
seeds <- 1:1000
count <- 0
result3_2_x1RS=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=1,PSgood=0,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=10,shift_param1=2)})
seeds <- 1:1000
count <- 0
result4_2_x1RS=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=0,PSgood=0,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=10,shift_param1=2)})
x1RS=data.frame(rbind(apply(result1_2_x1RS, 1, mean),apply(result2_2_x1RS, 1, mean),apply(result3_2_x1RS, 1, mean),apply(result4_2_x1RS, 1, mean)))
x1RS_sd=data.frame(rbind(apply(result1_2_x1RS, 1, sd),apply(result2_2_x1RS, 1, sd),apply(result3_2_x1RS, 1, sd),apply(result4_2_x1RS, 1, sd)))
x1RS$Misspecification=c("Neither","OR","PS","Both")
colnames(x1RS)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","Misspecification")
x1RS_sd$Misspecification=c("Neither","OR","PS","Both")
colnames(x1RS_sd)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","Misspecification")
rownames(x1RS)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
x1RS$shift="Normal(10,2)"

seeds <- 1:1000
count <- 0
result1_2_x1RS1=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=1,PSgood=1,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=10,shift_param1=3)})
seeds <- 1:1000
count <- 0
result2_2_x1RS1=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=0,PSgood=1,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=10,shift_param1=3)})
seeds <- 1:1000
count <- 0
result3_2_x1RS1=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=1,PSgood=0,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=10,shift_param1=3)})
seeds <- 1:1000
count <- 0
result4_2_x1RS1=replicate(1000, { count <<- count + 1; set.seed(seeds[count]); sim2_norm(ORgood=0,PSgood=0,n=10000,n_A=500,n_B=500,linear=TRUE,outcome="binary",shift=TRUE,shift1=TRUE,shift_param=10,shift_param1=3)})
x1RS1=data.frame(rbind(apply(result1_2_x1RS1, 1, mean),apply(result2_2_x1RS1, 1, mean),apply(result3_2_x1RS1, 1, mean),apply(result4_2_x1RS1, 1, mean)))
x1RS1_sd=data.frame(rbind(apply(result1_2_x1RS1, 1, sd),apply(result2_2_x1RS1, 1, sd),apply(result3_2_x1RS1, 1, sd),apply(result4_2_x1RS1, 1, sd)))
x1RS1$Misspecification=c("Neither","OR","PS","Both")
colnames(x1RS1)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","Misspecification")
x1RS1_sd$Misspecification=c("Neither","OR","PS","Both")
colnames(x1RS1_sd)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","Misspecification")
rownames(x1RS1)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
x1RS1$shift="Normal(10,3)"

df_cs=rbind(x1LS,x1LS1,x1RS,x1RS1)
df_cs_sd=rbind(x1LS_sd,x1LS1_sd,x1RS_sd,x1RS1_sd)
save(df_cs,file="bias_cs.RData")
save(df_cs,file="bias_sd_cs.RData")
# Melt to long format
df_cs_long <- melt(df_cs,
                   id.vars = c("Misspecification", "shift"),
                   variable.name = "Method",
                   value.name = "Estimate")
df_long_cs_sd <- melt(df_cs_sd,
                      id.vars = c("Misspecification"),
                      variable.name = "Method",
                      value.name = "SE")
df_cs_long$CI_lower <- df_cs_long$Estimate - 1.96 * df_long_cs_sd$SE
df_cs_long$CI_upper <- df_cs_long$Estimate + 1.96 * df_long_cs_sd$SE

#df_long$N=ifelse(df_long$N=="n=200","n=200 per study",ifelse(df_long$N=="n=500","n=500 per study","n=1000 per study"))
df_cs_long$shift=factor(df_cs_long$shift,levels = c("Normal(3,2)","Normal(3,3)",
                                                    "Normal(10,2)","Normal(10,3)"))
library(dplyr)
filtered_df_cs_long <- df_cs_long %>%
  filter(!(Method %in% c("IPW1", "IPW2", "MAIC1", "MAIC2", "RA", "STC") &
             Misspecification %in% c("Neither", "Both")))



ggplot(filtered_df_cs_long, aes(x = Method, y = Estimate, color = Misspecification)) +
  #geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  #geom_bar(stat = "identity", position = "dodge") +
  geom_point(position = position_dodge(width=0.9),size=2)+
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2, position = position_dodge(width = 0.9)) +
  facet_grid(~ shift,scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  theme_minimal(base_size = 13) +
  labs(x = "Method", y = "Bias") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))





ivShift1=data.frame(rbind(apply(result1_2, 1, mean),apply(result2_2, 1, mean),apply(result3_2, 1, mean),apply(result4_2, 1, mean)))
ivShift1_sd=data.frame(rbind(apply(result1_2, 1, sd),apply(result2_2, 1, sd),apply(result3_2, 1, sd),apply(result4_2, 1, sd)))
ivShift1$Misspecification=c("Neither","OR","PS","Both")
colnames(ivShift1)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","Misspecification")
ivShift1$iv_strength=0.5
ivShift1_sd$Misspecification=c("Neither","OR","PS","Both")
colnames(ivShift1_sd)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","Misspecification")
rownames(ivShift1)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
ivShift2$iv_strength=1
ivShift2_sd$Misspecification=c("Neither","OR","PS","Both")
colnames(ivShift2_sd)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","Misspecification")
rownames(ivShift2)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
ivShift3$iv_strength=0.75
ivShift3_sd$Misspecification=c("Neither","OR","PS","Both")
colnames(ivShift3_sd)=c("IPW1","IPW2","MAIC1","MAIC2","RA","STC","AIPW","DR","DR (park)","Misspecification")
rownames(ivShift3)=c("OR1PS1","OR0PS1","OR1PS0","OR0PS0")
df_iv=rbind(iv_baseline,ivShift1,ivShift2,ivShift3)
df_iv_sd=rbind(df2_2_sd,ivShift1_sd,ivShift2_sd,ivShift3_sd)
save(df_iv,file="bias_iv.RData")
save(df_iv_sd,file="bias_sd_iv.RData")

# Melt to long format
df_iv_long <- melt(df_iv,
                   id.vars = c("Misspecification", "iv_strength"),
                   variable.name = "Method",
                   value.name = "Estimate")
df_long_iv_sd <- melt(df_iv_sd,
                      id.vars = c("Misspecification"),
                      variable.name = "Method",
                      value.name = "SE")
df_iv_long$CI_lower <- df_iv_long$Estimate - 1.96 * df_long_iv_sd$SE
df_iv_long$CI_upper <- df_iv_long$Estimate + 1.96 * df_long_iv_sd$SE

#df_long$N=ifelse(df_long$N=="n=200","n=200 per study",ifelse(df_long$N=="n=500","n=500 per study","n=1000 per study"))

df_iv_long$iv_strength=factor(df_iv_long$iv_strength,labels = c("IV strength = 0.1","IV strength = 0.5","IV strength = 0.75","IV strength = 1"))
filtered_df_iv_long <- df_iv_long %>%
  filter(!(Method %in% c("IPW1", "IPW2", "MAIC1", "MAIC2", "RA", "STC") &
             Misspecification %in% c("Neither", "Both")))
ggplot(filtered_df_iv_long, aes(x = Method, y = Estimate, color = Misspecification)) +
  #geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  #geom_bar(stat = "identity", position = "dodge") +
  geom_point(position = position_dodge(width=0.9),size=2)+
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2, position = position_dodge(width = 0.9)) +
  facet_grid(~ iv_strength,scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  theme_minimal(base_size = 13) +
  labs(x = "Method", y = "Bias") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))