#' A function of power calculation for Stepped Wedge Design Studies
#' @description This function performs power calculations for stepped wedge cluster randomized trials under different scenarios.
#' @param K number of participants at each time period in a cluster
#' @param design I*J dimensional data set that describes the study design (control 0, intervention 1), I is the number of clusters, J is the number of time periods
#' @param family family of responses, specify family="gaussian" for continuous outcome and family="binomial" for binary outcome, with default value of "binomial"
#' @param model choose from conditional model (model="conditional") and marginal model (model="marginal"), with default value of applying conditional model
#' @param link choose link function from link="identity", link="log" and link="logit", with default value of identity link
#' @param type choose the study type, specify type="cohort" for closed cohort study and type="cross-sectional" for cross-sectional study, with default value of cross-sectional study
#' @param meanresponse_start the anticipated mean response rate in the control group at the start of the study
#' @param meanresponse_end0 the anticipated mean response rate in the control group at the end of the study, with default value equals to meanresponse_start (no time effects)
#' @param meanresponse_end1 the anticipated mean response rate in the intervention group at the end of the study
#' @param effectsize_beta the anticipated effect size, just omit this parameter if you don't need to specify it. In all scenarios, you can choose to specify the three parameters about mean responses without specifying this effect size, or alternatively specify meanresponse_start, meanresponse_end0 and this effect size. For continuous outcomes, users can conduct power calculations by only specifying this parameter without the above three parameters about mean responses (as the power is dependent just on it), then calculation will be implemented assuming scenarios without time effects. If you would consider scenarios with time effects and continuous outcomes, please specify meanresponse_start, meanresponse_end0 (donot require accurate information, just make sure they are not equal) and this effectsize_beta.
#' @param sigma2 marginal variance of the outcome (only needed for continuous outcomes and should not be an input for binary outcomes), with default value of 0.
#' @param typeIerror two-sided type I error, with default value of 0.05
#' @param alpha0 within-period correlation, with default value of 0.1
#' @param alpha1 between-period correlation, with default value of alpha0/2
#' @param alpha2 within-individual correlation, should not be an input under cross-sectional designs although it is numerically identical to alpha1 in this scenario by definition
#' @return The object returned is a list, which includes the design matrix and a summary table of this design (including the power)
#' @examples
#' library(swdpwr)
#' #a cross-sectional design with 12 clusters, 3 periods and binary outcomes applying conditional model
#' #alpha2 should not be specified, as the current version does not support power calculation using
#' #conditional models with binary outcomes in a cohort design
#' #create a 12*3 matrix which describes the study design,
#' #0 means control status, 1 means intervention status
#' dataset = matrix(c(rep(c(0,1,1),6),rep(c(0,0,1),6)),12,3,byrow=TRUE)
#'
#' #specify meanresponse_start, meanresponse_end0 and meanresponse_end1
#' swdpower(K = 30, design = dataset, family = "binomial", model = "conditional", link = "logit",
#' type = "cross-sectional", meanresponse_start = 0.2, meanresponse_end0 = 0.3,
#' meanresponse_end1 = 0.4, typeIerror = 0.05, alpha0 = 0.01, alpha1 = 0.01)
#'
#' #specify meanresponse_start, meanresponse_end0 and effectsize_beta
#' swdpower(K = 30, design = dataset, family = "binomial", model = "conditional", link = "logit",
#' type = "cross-sectional", meanresponse_start = 0.2, meanresponse_end0 = 0.3, effectsize_beta = 0.6,
#' typeIerror = 0.05, alpha0 = 0.01, alpha1 = 0.01)
#'
#' #a cohort design with 8 clusters, 3 periods and continuous outcomes applying marginal model
#' #sigma2 should be specified, as continuous outcomes require marginal variance in calculation
#' #create a 8*3 matrix which describes the study design,
#' #0 means control status, 1 means intervention status
#' dataset = matrix(c(rep(c(0,1,1),4),rep(c(0,0,1),4)),8,3, byrow=TRUE)
#'
#' #specify meanresponse_start, meanresponse_end0 and meanresponse_end1
#' swdpower(K = 24, design = dataset, family = "gaussian", model = "marginal", link = "identity",
#' type = "cohort", meanresponse_start = 0.1, meanresponse_end0 = 0.2,  meanresponse_end1 = 0.4,
#' sigma2 = 0.095, typeIerror = 0.05, alpha0 = 0.03, alpha1 = 0.015, alpha2 = 0.2)
#'
#' #specify effectsize_beta only, then the program runs assuming no time effects
#' swdpower(K = 24, design = dataset, family = "gaussian", model = "marginal", link = "identity",
#' type = "cohort",effectsize_beta=0.3, sigma2 = 0.095, typeIerror = 0.05, alpha0 = 0.03,
#' alpha1 = 0.015, alpha2 = 0.2)
#' @export

swdpower<-function(K, design, family="binomial", model="conditional", link="identity", type="cross-sectional", meanresponse_start=NA, meanresponse_end0=meanresponse_start, meanresponse_end1=NA, effectsize_beta=NA, sigma2=0, typeIerror=0.05, alpha0=0.1, alpha1=alpha0/2, alpha2=NA)
{
  ICC0=alpha0
  ICC1=alpha1
  ICC2=alpha2
  if(type!="cohort"&&type!="cross-sectional") stop("Wrong input for study type.")
  if(model=="conditional"&&type=="cohort"&&family=="binomial") {
    warning("For binary outcomes, conditional model only allows for cross-sectional settings.'type='cross-sectional'' is forced.")
    type="cross-sectional"}
  if(type=="cross-sectional"&is.na(ICC2)==0) warning("For a cross-sectional study, alpha2 is undefined. So please do not specify it although it equals to alpha1 by definition.")
  if(type=="cross-sectional")  ICC2=ICC1
  if(family=="binomial"&&sigma2!=0)  warning("Argument 'sigma2=' is not used because 'family='binomial'' is specified. For binary outcomes, the marginal variance of outcome sigma2 should not be an input.")
  if(family=="gaussian") {response=1} else
    if(family=="binomial") {response=2} else
    {stop("Wrong input of class for family of responses.")}
  if(family=="gaussian"&&link!="identity") {
    warning("Continuous outcomes only allow for identity link function and the link function is forced to be identity.")
    link="identity"
  }
  lf=link
  if(link=="identity") {link=1} else
    if(link=="log") {link=2} else
      if(link=="logit") {link=3} else
      {stop("Wrong input of class for link functions.")}

  if(model=="conditional"&&type=="cross-sectional"&&family=="binomial"&&(ICC0!=ICC1)) {
    warning("For cross-sectional setting and binary outcome, conditional model requires alpha0=alpha1 as time effects are fixed. Value of 'alpha0' is set to 'alpha1'.")
    ICC0=ICC1}

  md=model
  if(model=="conditional") {model=1} else
    if(model=="marginal") {model=2} else
    {stop("Wrong input of class for models.")}
  if(type=="cohort"&&is.na(ICC2)==1) stop("ICC2 should be specified for cohort study")

  if(!is.na(meanresponse_start)&&!is.na(meanresponse_end0)&&!is.na(meanresponse_end1)&&!is.na(effectsize_beta)) stop("Don't specify all of meanresponse_start, meanresponse_end0, meanresponse_end1, effectsize_beta at the same time. Only two ways of parameter specification are allowed: 1)specify: meanresponse_start, meanresponse_end0, meanresponse_end1; 2)specify: meanresponse_start, meanresponse_end0, effectsize_beta.")
  if(family=="gaussian"&&!is.na(effectsize_beta)&&!is.na(meanresponse_end1-meanresponse_end0)&&((meanresponse_end1-meanresponse_end0)!=effectsize_beta)) {stop("The effect size has contradictive input as meanresponse_end1-meanresponse_end0 is not equal to effectsize_beta.")}
  if(family=="gaussian"&!is.na(effectsize_beta)) {
    if(!is.na(meanresponse_end0-meanresponse_start)) {meanresponse_end1=meanresponse_end0+effectsize_beta}
    else {
      meanresponse_start=0
      meanresponse_end0=0
      meanresponse_end1=meanresponse_end0+effectsize_beta}
  }
  if(family=="binomial"&&is.na(meanresponse_end1)) {
    if(link==1) meanresponse_end1=meanresponse_end0+effectsize_beta
    if(link==2) meanresponse_end1=(meanresponse_end0)*exp(effectsize_beta)
    if(link==3) meanresponse_end1=((1/(1-meanresponse_end0)-1)*exp(effectsize_beta))/(1+((1/(1-meanresponse_end0)-1)*exp(effectsize_beta)))
  }

  if(!is.na(meanresponse_start)*!is.na(meanresponse_end0)*!is.na(meanresponse_end1)) {}
  else stop("Missing input parameter for effect size calculation.")
  m1=meanresponse_start
  m2=meanresponse_end0
  m3=meanresponse_end1
  tem1=min(m1,m2,m3)
  tem2=max(m1,m2,m3)
  if(family=="binomial"&&tem1<0) stop("Violation of valid probability, given input parameters: min(meanresponse_start,meanresponse_end0, meanresponse_end1)<0. Please check whether any of these values are out of range and revise one or more of them.")
  if(family=="binomial"&&tem2>1) stop("Violation of valid probability, given input parameters: max(meanresponse_start,meanresponse_end0, meanresponse_end1)>1. Please check whether any of these values are out of range and revise one or more of them.")

  if(family=="binomial"&&model==1&&(K>150)&&(m1!=m2)) stop("K should be at least smaller than 150 for this scenario as the running time is too long with this K for the power calculation of binary outcomes under conditional model with time effects. Please reduce K or use the model without time effects or use marginal models.")

  dataset=design
  I=dim(dataset)[1]
  J=dim(dataset)[2]
  II=I
  JJ=J
  KK=K
  X_in=dataset
  res=response
  opt=model
  if(link==1) {
    mu=m1
    gammaJ=m2-m1
    beta=m3-m2
  }
  if(link==2) {
    mu=log(m1)
    gammaJ=log(m2/m1)
    beta=log(m3/m2)
  }
  if(link==3) {
    mu=log(1/(1-m1)-1)
    gammaJ=log(1/(1-m2)-1)-log(1/(1-m1)-1)
    beta=log(1/(1-m3)-1)-log(1/(1-m2)-1)
  }

  alpha=typeIerror
  mu=mu
  beta=beta
  p0totalchange=gammaJ
  sigma2=sigma2
  rho0=ICC1
  typeone=alpha
  link=link-1


  TOLERANCE = 1e-5
  convergence = 0
  p0 = rep(0,JJ)
  gamma = rep(0,JJ)
  mincomp = rep(0,JJ+2)
  maxcomp = rep(0,JJ+2)
  #initialization
  tau2=1
  power=0



  temp1=c(ICC0,ICC1,ICC2)
  if(max(temp1)>1) stop('Violate range of intraclass correlations: max(alpha0,alpha1,alpha2)>1.Please correct the values of correlation parameters, they must be between 0 and 1.')

  temp1=c(ICC0,ICC1,ICC2)
  if(min(temp1)<0) stop('Violate range of intraclass correlations: min(alpha0,alpha1,alpha2)<0.Please correct the values of correlation parameters, they must be between 0 and 1.')

  if(typeone>1) stop("Type I error provided is larger than 1, it must be between 0 and 1, and is usually 0.05.")
  if(typeone<0) stop("Type I error provided is less than 0, it must be between 0 and 1, and is usually 0.05.")
  alpha0= ICC0
  alpha1= rho0
  alpha2= ICC2
  if(res>1.5){

    if(opt<1.5)
    {
      if (p0totalchange > TOLERANCE | p0totalchange < -TOLERANCE)
      {
        if(link<1)
        {
          GQ=100
          GQX = rep(0,GQ)
          GQW = rep(0,GQ)
          p01=mu
          p11=mu+beta
          p0[1] = p01
          p11 = p11
          p0stepchange = p0totalchange/(JJ-1)
          for(j in 2:JJ)
          {
            p0[j]= p0[j-1] + p0stepchange
          }
          #  call computeparameter(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
          #dyn.load("computeparameter.so")
          mu=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$mu)
          beta=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$beta)
          gamma=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$gamma)
          tau2=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$tau2)

          period=NULL
          period=gamma+mu
          period[JJ+1]=beta
          X=matrix(0,JJ,JJ+1)
          X[,1:JJ]=diag(JJ)
          for(i in 1:II)
          {
            X[,JJ+1]=X_in[i,]
            gmu=X%*%as.matrix(period)
            mmu=(gmu*(1-gmu))^0.5
            for(m in 1:JJ)
            {
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))>gmu[m]) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))>gmu(m)')
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))<max(2*gmu[m]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")# stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))<max(0,2*gmu(m)-1)')
            }
            for(m in 1:(JJ-1))
              for(k in (m+1):JJ)
              {
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")# stop('((gmu(m)*gmu(k)+alpha1*(1/mmu(m,m))*(1/mmu(k,k))))>min(gmu(m),gmu(k))')
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu(m)*gmu(k)+alpha1*(1/mmu(m,m))*(1/mmu(k,k))))<max(0,gmu(m)+gmu(k)-1)')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu(m)*gmu(k)+alpha2*(1/mmu(m,m))*(1/mmu(k,k))))>min(gmu(m),gmu(k))')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu(m)*gmu(k)+alpha2*(1/mmu(m,m))*(1/mmu(k,k))))<max(0,gmu(m)+gmu(k)-1)')
              }
          }



          temp1=c(p01+p0totalchange,p11+p0totalchange,p01,p11)
          if(min(temp1)<0) stop("Violate theory of probability under identity link: min(mu+gammaJ,mu+beta+gammaJ,mu,mu+beta)<0. Please check whether any of these four values are out of range.")

          temp1=c(p01+p0totalchange,p11+p0totalchange,p01,p11)
          if(max(temp1)>1) stop('Violate theory of probability under identity link: max(mu+gammaJ,mu+beta+gammaJ,mu,mu+beta)>1. Please check whether any of these four values are out of range.')

          a = 100
          b = -100
          for(j in 1:JJ)
          {
            temp = mu+gamma[j]
            if (temp < a){
              a = temp
              mincomp=rep(0,JJ+2)
              mincomp[JJ+1] = 1
              mincomp[j] = 1
            }

            if (temp > b){
              b = temp
              maxcomp=rep(0,JJ+2)
              maxcomp[JJ+1] = 1
              maxcomp[j] = 1
            }
            temp = mu+beta+gamma[j]
            if (temp < a){
              a = temp
              mincomp= rep(0,JJ+2)
              mincomp[JJ+1] = 1
              mincomp[JJ+2] = 1
              mincomp[j] = 1
            }

            if (temp > b){
              b = temp
              maxcomp= rep(0,JJ+2)
              maxcomp[JJ+1] = 1
              maxcomp[JJ+2] = 1
              maxcomp[j] = 1
            }
          }
          a = -a
          b = 1-b
          #call legendre_handle (GQ, a, b, GQX, GQW)
          ##dyn.load("legendre_rule.so")
          GQX=suppressWarnings(.Fortran("legendrehandle",GQ=as.integer(GQ), a=as.numeric(a), b=as.numeric(b), GQX=as.numeric(GQX), GQW=as.numeric(GQW) )$GQX)
          GQW=suppressWarnings(.Fortran("legendrehandle",GQ=as.integer(GQ), a=as.numeric(a), b=as.numeric(b), GQX=as.numeric(GQX), GQW=as.numeric(GQW) )$GQW)
          #Gaussian Legendre will not take two limits, a and b
          #power = LinearPower_time(mu, beta, gamma, tau2, II, JJ, KK, a, b, mincomp, maxcomp, GQ, GQX, GQW, X_in)
          #dyn.load("power_cal_wrapper.so")
          power=suppressWarnings(.Fortran("LinearPowertimewrapper",mu=mu, beta, gamma=gamma, tau2=tau2, II=as.integer(II), JJ=as.integer(JJ), KK=as.integer(KK), a=as.numeric(a), b=as.numeric(b), mincomp=as.integer(mincomp), maxcomp=as.integer(maxcomp), GQ=as.integer(GQ), GQX=as.numeric(GQX), GQW=as.numeric(GQW), X_in=as.integer(X_in),typeone=as.numeric(typeone),power=power )$power)

        }
        else if(link<2 & link>0)
        {
          GQ=500
          GQX = rep(0,GQ)
          GQW = rep(0,GQ)
          gamma = rep(0,JJ)
          p01= exp(mu)
          p0[1] = p01
          p11 = exp(mu+beta)
          p11 = p11
          p0stepchange = p0totalchange/(JJ-1)
          for(j in 2:JJ)
          {
            gamma[j]= gamma[j-1] + p0stepchange
          }
          for(j in 2:JJ)
          {
            p0[j]= exp(gamma[j]+mu)
          }
          #  call computeparameter(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
          #dyn.load("computeparameterlog.so")
          mu=suppressWarnings(.Fortran("computeparameterlog",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$mu)
          beta=suppressWarnings(.Fortran("computeparameterlog",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$beta)
          gamma=suppressWarnings(.Fortran("computeparameterlog",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$gamma)
          tau2=suppressWarnings(.Fortran("computeparameterlog",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$tau2)

          period=NULL
          period=gamma+mu
          period[JJ+1]=beta
          X=matrix(0,JJ,JJ+1)
          X[,1:JJ]=diag(JJ)
          for(i in 1:II)
          {
            X[,JJ+1]=X_in[i,]
            gmu=X%*%as.matrix(period)
            gmu=exp(gmu)
            mmu=(gmu*(1-gmu))^0.5
            for(m in 1:JJ)
            {
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))>gmu[m]) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))>gmu(m)')
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))<max(2*gmu[m]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))<max(0,2*gmu(m)-1)')
            }
            for(m in 1:(JJ-1))
              for(k in (m+1):JJ)
              {
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu(m)*gmu(k)+alpha1*(1/mmu(m,m))*(1/mmu(k,k))))>min(gmu(m),gmu(k))')
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0))  stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu(m)*gmu(k)+alpha1*(1/mmu(m,m))*(1/mmu(k,k))))<max(0,gmu(m)+gmu(k)-1)')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu(m)*gmu(k)+alpha2*(1/mmu(m,m))*(1/mmu(k,k))))>min(gmu(m),gmu(k))')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu(m)*gmu(k)+alpha2*(1/mmu(m,m))*(1/mmu(k,k))))<max(0,gmu(m)+gmu(k)-1)')
              }
          }


          temp1=c(mu+p0totalchange,mu+beta+p0totalchange,mu,mu+beta)
          if(max(temp1)>0) stop('Violate theory of probability under log link: max(mu+gammaJ,mu+beta+gammaJ,mu,mu+beta)>0.Please check whether any of these four values are out of range.')




          a = 100
          b = -100
          for(j in 1:JJ)
          {
            temp = mu+gamma[j]
            if (temp < a){
              a = temp
              mincomp= rep(0,JJ+2)
              mincomp[JJ+1] = 1
              mincomp[j] = 1
            }

            if (temp > b){
              b = temp
              maxcomp= rep(0,JJ+2)
              maxcomp[JJ+1] = 1
              maxcomp[j] = 1
            }
            temp = mu+beta+gamma[j]
            if (temp < a){
              a = temp
              mincomp= rep(0,JJ+2)
              mincomp[JJ+1] = 1
              mincomp[JJ+2] = 1
              mincomp[j] = 1
            }

            if (temp > b){
              b = temp
              maxcomp= rep(0,JJ+2)
              maxcomp[JJ+1] = 1
              maxcomp[JJ+2] = 1
              maxcomp[j] = 1
            }
          }
          if (beta>0)
          {
            a= 5000
          }else
          { a= 5000 }

          b = b
          #call legendre_handle (GQ, b, a, GQX, GQW)
          #!call legendre_handle2 (GQ, b, a, GQX, GQW)
          #dyn.load("legendre_rule.so")
          GQX=suppressWarnings(.Fortran("legendrehandle",GQ=as.integer(GQ), a=b, b=a, GQX=as.numeric(GQX), GQW=as.numeric(GQW) )$GQX)
          GQW=suppressWarnings(.Fortran("legendrehandle",GQ=as.integer(GQ), a=b, b=a, GQX=as.numeric(GQX), GQW=as.numeric(GQW) )$GQW)

          #Gaussian Legendre will not take two limits, a and b
          #power = LogPower_time(mu, beta, gamma, tau2, II, JJ, KK, b, a, mincomp, maxcomp, GQ, GQX, GQW, X_in)
          #dyn.load("power_cal_wrapper.so")
          power=suppressWarnings(.Fortran("LogPowertimewrapper",mu=mu, beta=beta, gamma=gamma, tau2=tau2, II=as.integer(II), JJ=as.integer(JJ), KK=as.integer(KK), a=b, b=a, mincomp=as.integer(mincomp), maxcomp= as.integer(maxcomp), GQ=as.integer(GQ),GQX=as.numeric(GQX), GQW=as.numeric(GQW), X_in=as.integer(X_in),typeone=as.numeric(typeone),power=power )$power)

        }else
        {
          GQ = 40
          GQX = rep(0,GQ)
          GQW = rep(0,GQ)
          #call HERZO(GQ,GQX,GQW)
          #dyn.load("herzo.so")
          GQW=suppressWarnings(.Fortran("HERZO",GQ=as.integer(GQ),GQX=as.numeric(GQX),GQW=as.numeric(GQX) )$GQW)
          GQX=suppressWarnings(.Fortran("HERZO",GQ=as.integer(GQ),GQX=as.numeric(GQX),GQW=as.numeric(GQX) )$GQX)
          p01=exp(mu)/(1+exp(mu))
          p11=exp(mu+beta)/(1+exp(mu+beta))
          p0[1] = p01
          p11 = p11
          p0stepchange = p0totalchange/(JJ-1)
          for(j in 2:JJ)
          {
            gamma[j]= gamma[j-1] + p0stepchange
          }
          for(j in 2:JJ)
          {
            p0[j]= exp(gamma[j]+mu)/(1+exp(mu+gamma[j]))
          }
          #initial values
          mual=mu
          betaal=beta
          gammaal=gamma
          tau2 = 1.0
          #gamma = 0.0
          #call computeparameterlogit(JJ, mu, beta, gamma, tau2, p0, p11, rho0, GQ, GQX, GQW, convergence)
          #dyn.load("computeparameterlogit.so")
          mu=suppressWarnings(.Fortran("computeparameterlogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0, GQ=as.integer(GQ), GQX=as.numeric(GQX), GQW=as.numeric(GQW), convergence=convergence)$mu)
          beta=suppressWarnings(.Fortran("computeparameterlogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0, GQ=as.integer(GQ), GQX=as.numeric(GQX), GQW=as.numeric(GQW), convergence=convergence)$beta)
          gamma=suppressWarnings(.Fortran("computeparameterlogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0, GQ=as.integer(GQ), GQX=as.numeric(GQX), GQW=as.numeric(GQW), convergence=convergence)$gamma)
          tau2=suppressWarnings(.Fortran("computeparameterlogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0, GQ=as.integer(GQ), GQX=as.numeric(GQX), GQW=as.numeric(GQW), convergence=convergence)$tau2)

          period=NULL
          period=gamma+mu
          period[JJ+1]=beta
          X=matrix(0,JJ,JJ+1)
          X[,1:JJ]=diag(JJ)
          for(i in 1:II)
          {
            X[,JJ+1]=X_in[i,]
            gmu=X%*%as.matrix(period)
            gmu=exp(gmu)/(1+exp(gmu))
            mmu=(gmu*(1-gmu))^0.5
            for(m in 1:JJ)
            {
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))>gmu[m]) stop('Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.')
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))<max(2*gmu[m]-1,0)) stop('Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.')
            }
            for(m in 1:(JJ-1))
              for(k in (m+1):JJ)
              {
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")# stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")# stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
              }
          }

          #mu=mual
          #beta=betaal
          #gamma=gammaal
          #power = LogitPower_time(mu, beta, gamma, tau2, II, JJ, KK, GQ, GQX, GQW,X_in)
          #dyn.load("power_cal_wrapper.so")
          power=suppressWarnings(.Fortran("LogitPowertimewrapper",mu=mu, beta=beta, gamma=gamma, tau2=tau2, II=as.integer(II), JJ=as.integer(JJ), KK=as.integer(KK), GQ=as.integer(GQ), GQX=as.numeric(GQX), GQW=as.numeric(GQW), X_in=as.integer(X_in),typeone=as.numeric(typeone), power=power )$power)
          if(!is.na(effectsize_beta)) beta=effectsize_beta
        }
      }
      else
      {
        if(link<1)
        {

          GQ=100
          GQX = rep(0,GQ)
          GQW = rep(0,GQ)
          p01=mu
          p11=beta+mu
          p0[1] = p01
          p11 = p11
          p0stepchange = p0totalchange/(JJ-1)
          for(j in 2:JJ)
          {
            p0[j]= p0[j-1] + p0stepchange
          }
          #call computeparameter(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
          #dyn.load("computeparameter.so")
          mu=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$mu)
          beta=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$beta)
          tau2=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$tau2)

          period=NULL
          period[1]=mu
          period[2]=beta
          X=matrix(1,JJ,2)

          for(i in 1:II)
          {
            X[,2]=X_in[i,]
            gmu=X%*%as.matrix(period)
            mmu=(gmu*(1-gmu))^0.5
            for(m in 1:JJ)
            {
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))>gmu[m]) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")# stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))>gmu(m)')
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))<max(2*gmu[m]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))<max(0,2*gmu(m)-1)')
            }
            for(m in 1:(JJ-1))
              for(k in (m+1):JJ)
              {
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")# stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
              }
          }

          temp1=c(p01+p0totalchange,p11+p0totalchange,p01,p11)
          if(min(temp1)<0) stop("Violate theory of probability under identity link: min(mu+gammaJ,mu+beta+gammaJ,mu,mu+beta)<0.Please check whether any of these four values are out of range.")

          temp1=c(p01+p0totalchange,p11+p0totalchange,p01,p11)
          if(max(temp1)>1) stop('Violate theory of probability under identity link: max(mu+gammaJ,mu+beta+gammaJ,mu,mu+beta)>1.Please check whether any of these four values are out of range.')

          if (beta>0){
            a = - mu
            b = 1 - mu - beta
          }else
          {
            a = - mu - beta
            b = 1 - mu
          }
          #! a = a / dsqrt(2.0d0*tau2)
          #! b = b / dsqrt(2.0d0*tau2)
          #call legendre_handle (GQ, a, b, GQX, GQW)
          #! Gaussian Legendre will not take two limits, a and b
          #power = LinearPower_notime(mu, beta, tau2, II, JJ, KK, a, b, GQ, GQX, GQW, X_in)
          #dyn.load("legendre_rule.so")
          GQX=suppressWarnings(.Fortran("legendrehandle",GQ=as.integer(GQ), a=a, b=b, GQX=as.numeric(GQX), GQW=as.numeric(GQW) )$GQX)
          GQW=suppressWarnings(.Fortran("legendrehandle",GQ=as.integer(GQ), a=a, b=b, GQX=as.numeric(GQX), GQW=as.numeric(GQW) )$GQW)
          #dyn.load("power_cal_wrapper.so")
          power=suppressWarnings(.Fortran("LinearPowernotimewrapper",mu=mu, beta=beta,  tau2=tau2, II=as.integer(II), JJ=as.integer(JJ), KK=as.integer(KK), a=a, b=b, GQ=as.integer(GQ), GQX=as.numeric(GQX), GQW=as.numeric(GQW), X_in=as.integer(X_in) ,typeone=as.numeric(typeone), power=power )$power)
        }
        else if(link<2 & link>0)
        {

          GQ=500
          GQX = rep(0,GQ)
          GQW = rep(0,GQ)
          p01= exp(mu)
          p0[1] = p01
          p11 = exp(mu+beta)
          p11 = p11
          p0stepchange = p0totalchange/(JJ-1)
          for(j in 2:JJ)
          {
            gamma[j]= gamma[j-1] + p0stepchange
          }
          for(j in 2:JJ)
          {
            p0[j] = exp(gamma[j]+mu)
          }

          temp1=c(mu+p0totalchange,mu+beta+p0totalchange,mu,mu+beta)
          if(max(temp1)>0) stop('Violate theory of probability under log link: max(mu+gammaJ,mu+beta+gammaJ,mu,mu+beta)>0.Please check whether any of these four values are out of range.')

          #call computeparameterlog(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
          #dyn.load("computeparameterlog.so")
          mu=suppressWarnings(.Fortran("computeparameterlog",JJ=as.integer(JJ),mu=mu,beta=beta,gamma=gamma,tau2=tau2,p0=p0,p11=p11, rho0=rho0 )$mu)
          beta=suppressWarnings(.Fortran("computeparameterlog",JJ=as.integer(JJ), mu=mu, beta=beta,gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$beta)
          tau2=suppressWarnings(.Fortran("computeparameterlog",JJ=as.integer(JJ), mu=mu, beta=beta,gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$tau2)

          period=NULL
          period[1]=mu
          period[2]=beta
          X=matrix(1,JJ,2)
          for(i in 1:II)
          {
            X[,2]=X_in[i,]
            gmu=X%*%as.matrix(period)
            gmu=exp(gmu)
            mmu=(gmu*(1-gmu))^0.5
            for(m in 1:JJ)
            {
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))>gmu[m]) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")# stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))>gmu(m)')
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))<max(2*gmu[m]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))<max(0,2*gmu(m)-1)')
            }
            for(m in 1:(JJ-1))
              for(k in (m+1):JJ)
              {
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0))  stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
              }
          }



          if (beta>0) {
            #!inf= 0
            #!a = 0.5d0/tau2
            a= 5000.0
            #!a=0.05
            b = - mu - beta

          }else{
            # !inf= 0
            #!a =0.5d0/tau2
            a= 5000.0
            #!a=0.05
            b =  - mu
          }
          #! a = a / dsqrt(2.0d0*tau2)
          #! b = b / dsqrt(2.0d0*tau2)
          b=-b
          # !b=0
          # call legendre_handle (GQ, b, a, GQX, GQW)
          # !  call legendre_handle2 (GQ, b, a, GQX, GQW)
          # ! Gaussian Legendre will not take two limits, a and b
          # power = LogPower_notime(mu, beta, tau2, II, JJ, KK, b, a, GQ, GQX, GQW, X_in)
          #dyn.load("legendre_rule.so")
          GQX=suppressWarnings(.Fortran("legendrehandle",GQ=as.integer(GQ), a=b, b=a, GQX=as.numeric(GQX), GQW=as.numeric(GQW) )$GQX)
          GQW=suppressWarnings(.Fortran("legendrehandle",GQ=as.integer(GQ), a=b, b=a, GQX=as.numeric(GQX), GQW=as.numeric(GQW) )$GQW)
          #Gaussian Legendre will not take two limits, a and b
          #dyn.load("power_cal_wrapper.so")
          power=suppressWarnings(.Fortran("LogPowernotimewrapper",mu=mu, beta=beta, tau2=tau2, II=as.integer(II), JJ=as.integer(JJ), KK=as.integer(KK), a=b,b=a,GQ=as.integer(GQ),GQX=as.numeric(GQX), GQW=as.numeric(GQW), X_in=as.integer(X_in) ,typeone=as.numeric(typeone), power=power )$power)

        }
        else
        {
          GQ = 40
          GQX = rep(0,GQ)
          GQW = rep(0,GQ)
          #call HERZO(GQ,GQX,GQW)
          #dyn.load("herzo.so")
          GQX=suppressWarnings(.Fortran("HERZO",GQ=as.integer(GQ),GQX=as.numeric(GQX),GQW=as.numeric(GQW) )$GQX)
          GQW=suppressWarnings(.Fortran("HERZO",GQ=as.integer(GQ),GQX=as.numeric(GQX),GQW=as.numeric(GQW) )$GQW)
          p0totalchange=0
          p01=exp(mu)/(1+exp(mu))
          p11=exp(mu+beta)/(1+exp(mu+beta))
          p0[1] = p01
          p11 = p11
          p0stepchange = p0totalchange/(JJ-1)
          for(j in 2:JJ)
          {
            gamma[j]= gamma[j-1] + p0stepchange
          }
          for(j in 2:JJ)
          {
            p0[j]= exp(gamma[j]+mu)/(1+exp(mu+gamma[j]))
          }

          mual=mu
          betaal=beta
          gammaal=gamma
          tau2 = 1.0
          #gamma = 0.0
          #call computeparameterlogit(JJ, mu, beta, gamma, tau2, p0, p11, rho0, GQ, GQX, GQW, convergence)
          #power = LogitPower_notime(mu, beta, tau2, II, JJ, KK, GQ, GQX, GQW,X_in)
          #dyn.load("computeparameterlogit.so")
          mu=suppressWarnings(.Fortran("computeparameterlogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0,GQ=as.integer(GQ), GQX=as.numeric(GQX), GQW=as.numeric(GQW), convergence=convergence )$mu)
          beta=suppressWarnings(.Fortran("computeparameterlogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0,GQ=as.integer(GQ), GQX=as.numeric(GQX), GQW=as.numeric(GQW), convergence=convergence )$beta)
          tau2=suppressWarnings(.Fortran("computeparameterlogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0,GQ=as.integer(GQ), GQX=as.numeric(GQX), GQW=as.numeric(GQW), convergence=convergence )$tau2)

          period=NULL
          period[1]=mu
          period[2]=beta
          X=matrix(1,JJ,2)
          for(i in 1:II)
          {
            X[,2]=X_in[i,]
            gmu=X%*%as.matrix(period)
            gmu=exp(gmu)/(1+exp(gmu))
            mmu=(gmu*(1-gmu))^0.5
            for(m in 1:JJ)
            {
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))>gmu[m]) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))>gmu(m)')
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))<max(2*gmu[m]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))<max(0,2*gmu(m)-1)')
            }
            for(m in 1:(JJ-1))
              for(k in (m+1):JJ)
              {
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0))  stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
              }
          }

          #dyn.load("power_cal_wrapper.so")
          #mu=mual
          #beta=betaal
          #gamma=gammaal
          power=suppressWarnings(.Fortran("LogitPowernotimewrapper",mu=mu, beta=beta, tau2=tau2, II=as.integer(II), JJ=as.integer(JJ), KK=as.integer(KK), GQ=as.integer(GQ), GQX=as.numeric(GQX), GQW=as.numeric(GQW), X_in=as.integer(X_in), typeone=as.numeric(typeone),power=power )$power)
          if(!is.na(effectsize_beta)) beta=effectsize_beta
        }
      }
    }
    else
    {
      alpha0= ICC0
      alpha1= rho0
      alpha2= ICC2
      lambda1 = 1- alpha0 + alpha1 - alpha2
      lambda2 = 1- alpha0 - (JJ-1)*(alpha1 - alpha2)
      lambda3 = 1 + (KK-1)*(alpha0 - alpha1) - alpha2
      lambda4 = 1 + (KK-1)*alpha0 + (JJ-1)*(KK-1)*alpha1 + (JJ-1)*alpha2
      temp2=c(lambda1,lambda2,lambda3,lambda4)
      if(min(temp2)<0) stop('Correlation matrix R is no longer positive definite. Please check whether the between-period correlation is unrealistically larger than the within-period correlation or the within-individual correlation.')

      #if(ICC1>ICC0)  stop('ICC1>ICC0')


      #if(ICC1>ICC2)  stop('ICC1>ICC2')

      if (p0totalchange > TOLERANCE | p0totalchange < -TOLERANCE)
      {
        if(link<1)
        {
          p01=mu
          p11=mu+beta
          p0[1] = p01
          p11 = p11
          p0stepchange = p0totalchange/(JJ-1)
          for(j in 2:JJ)
          {
            p0[j]= p0[j-1] + p0stepchange
          }

          temp1=c(p01+p0totalchange,p11+p0totalchange,p01,p11)
          if(min(temp1)<0) stop("Violate theory of probability under identity link: min(mu+gammaJ,mu+beta+gammaJ,mu,mu+beta)<0.Please check whether any of these four values are out of range.")

          temp1=c(p01+p0totalchange,p11+p0totalchange,p01,p11)
          if(max(temp1)>1) stop('Violate theory of probability under identity link: max(mu+gammaJ,mu+beta+gammaJ,mu,mu+beta)>1.Please check whether any of these four values are out of range.')

          #call computeparameter(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
          #dyn.load("computeparameter.so")
          mu=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$mu)
          beta=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$beta)
          gamma=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$gamma)
          rho0=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$rho0)
          #power=LinearPower_GEE(mu, beta, gamma, rho0, II, JJ, KK, X_in)
          #dyn.load("power_cal_wrapper.so")
          period=NULL
          period=gamma+mu
          period[JJ+1]=beta
          X=matrix(0,JJ,JJ+1)
          X[,1:JJ]=diag(JJ)
          for(i in 1:II)
          {
            X[,JJ+1]=X_in[i,]
            gmu=X%*%as.matrix(period)
            mmu=(gmu*(1-gmu))^0.5
            for(m in 1:JJ)
            {
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))>gmu[m]) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))>gmu(m)')
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))<max(2*gmu[m]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")# stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))<max(0,2*gmu(m)-1)')
            }
            for(m in 1:(JJ-1))
              for(k in (m+1):JJ)
              {
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")# stop('((gmu(m)*gmu(k)+alpha1*(1/mmu(m,m))*(1/mmu(k,k))))>min(gmu(m),gmu(k))')
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu(m)*gmu(k)+alpha1*(1/mmu(m,m))*(1/mmu(k,k))))<max(0,gmu(m)+gmu(k)-1)')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu(m)*gmu(k)+alpha2*(1/mmu(m,m))*(1/mmu(k,k))))>min(gmu(m),gmu(k))')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu(m)*gmu(k)+alpha2*(1/mmu(m,m))*(1/mmu(k,k))))<max(0,gmu(m)+gmu(k)-1)')
              }
          }



          power=suppressWarnings(.Fortran("LinearPowerGEEwrapper",mu=mu, beta=beta, gamma=gamma, rho0=rho0, II=as.integer(II), JJ=as.integer(JJ), KK=as.integer(KK), X_in=as.integer(X_in),ICC0=ICC0,ICC2=ICC2,typeone=as.numeric(typeone),power=power )$power)

        }
        else if(link<2 & link>0)
        {
          p01=exp(mu)
          p11=exp(mu+beta)
          p0[1] = p01
          p11 = p11
          p0stepchange = p0totalchange/(JJ-1)
          for(j in 2:JJ)
          {
            gamma[j]= gamma[j-1] + p0stepchange
          }
          for(j in 2:JJ)
          {
            p0[j] = exp(gamma[j]+mu)
          }
          temp1=c(mu+p0totalchange,mu+beta+p0totalchange,mu,mu+beta)
          if(max(temp1)>0) stop('Violate theory of probability under log link: max(mu+gammaJ,mu+beta+gammaJ,mu,mu+beta)>0.Please check whether any of these four values are out of range.')

          #call computeparameterGEElog(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
          #dyn.load("computeparameterGEElog.so")
          mu=suppressWarnings(.Fortran("computeparameterGEElog",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$mu)
          beta=suppressWarnings(.Fortran("computeparameterGEElog",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$beta)
          gamma=suppressWarnings(.Fortran("computeparameterGEElog",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$gamma)
          rho0=suppressWarnings(.Fortran("computeparameterGEElog",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$rho0)
          #power=LogPower_GEE(mu, beta, gamma, rho0, II, JJ, KK, X_in)
          #dyn.load("power_cal_wrapper.so")
          period=NULL
          period=gamma+mu
          period[JJ+1]=beta
          X=matrix(0,JJ,JJ+1)
          X[,1:JJ]=diag(JJ)
          for(i in 1:II)
          {
            X[,JJ+1]=X_in[i,]
            gmu=X%*%as.matrix(period)
            gmu=exp(gmu)
            mmu=(gmu*(1-gmu))^0.5
            for(m in 1:JJ)
            {
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))>gmu[m]) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))>gmu(m)')
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))<max(2*gmu[m]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))<max(0,2*gmu(m)-1)')
            }
            for(m in 1:(JJ-1))
              for(k in (m+1):JJ)
              {
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu(m)*gmu(k)+alpha1*(1/mmu(m,m))*(1/mmu(k,k))))>min(gmu(m),gmu(k))')
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0))  stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu(m)*gmu(k)+alpha1*(1/mmu(m,m))*(1/mmu(k,k))))<max(0,gmu(m)+gmu(k)-1)')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu(m)*gmu(k)+alpha2*(1/mmu(m,m))*(1/mmu(k,k))))>min(gmu(m),gmu(k))')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu(m)*gmu(k)+alpha2*(1/mmu(m,m))*(1/mmu(k,k))))<max(0,gmu(m)+gmu(k)-1)')
              }
          }


          power=suppressWarnings(.Fortran("LogPowerGEEwrapper",mu=mu, beta=beta, gamma=gamma, rho0=rho0, II=as.integer(II), JJ=as.integer(JJ), KK=as.integer(KK), X_in=as.integer(X_in),ICC0=ICC0,ICC2=ICC2,typeone=as.numeric(typeone),power=power )$power)
        }
        else
        {
          p01=exp(mu)/(1+exp(mu))
          p11=exp(mu+beta)/(1+exp(mu+beta))
          p0[1] = p01
          p11 = p11
          p0stepchange = p0totalchange/(JJ-1)
          for(j in 2:JJ)
          {
            gamma[j]= gamma[j-1] + p0stepchange
          }
          for(j in 2:JJ)
          {
            p0[j]= exp(gamma[j]+mu)/(1+exp(mu+gamma[j]))
          }
          #call computeparameterGEElogit(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
          #power=LogitPower_GEE(mu, beta, gamma, rho0, II, JJ, KK, X_in)
          #dyn.load("computeparameterGEElogit.so")
          mu=suppressWarnings(.Fortran("computeparameterGEElogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$mu)
          beta=suppressWarnings(.Fortran("computeparameterGEElogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$beta)
          gamma=suppressWarnings(.Fortran("computeparameterGEElogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$gamma)
          rho=suppressWarnings(.Fortran("computeparameterGEElogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$rho0)
          #dyn.load("power_cal_wrapper.so")
          period=NULL
          period=gamma+mu
          period[JJ+1]=beta
          X=matrix(0,JJ,JJ+1)
          X[,1:JJ]=diag(JJ)
          for(i in 1:II)
          {
            X[,JJ+1]=X_in[i,]
            gmu=X%*%as.matrix(period)
            gmu=exp(gmu)/(1+exp(gmu))
            mmu=(gmu*(1-gmu))^0.5
            for(m in 1:JJ)
            {
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))>gmu[m]) stop('Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.')
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))<max(2*gmu[m]-1,0)) stop('Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.')
            }
            for(m in 1:(JJ-1))
              for(k in (m+1):JJ)
              {
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")# stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")# stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
              }
          }

          power=suppressWarnings(.Fortran("LogitPowerGEEwrapper",mu=mu, beta=beta, gamma=gamma, rho0=rho0, II=as.integer(II), JJ=as.integer(JJ), KK=as.integer(KK), X_in=as.integer(X_in),ICC0=ICC0,ICC2=ICC2,typeone=as.numeric(typeone),power= power )$power)
        }
      }
      else
      {
        if(link<1)
        {
          p01=mu
          p11=beta+mu
          p0[1] = p01
          p11 = p11
          p0stepchange = p0totalchange/(JJ-1)
          for(j in 2:JJ)
          {
            p0[j]= p0[j-1] + p0stepchange
          }
          temp1=c(p01+p0totalchange,p11+p0totalchange,p01,p11)
          if(min(temp1)<0) stop("Violate theory of probability under identity link: min(mu+gammaJ,mu+beta+gammaJ,mu,mu+beta)<0.Please check whether any of these four values are out of range.")

          temp1=c(p01+p0totalchange,p11+p0totalchange,p01,p11)
          if(max(temp1)>1) stop('Violate theory of probability under identity link: max(mu+gammaJ,mu+beta+gammaJ,mu,mu+beta)>1.Please check whether any of these four values are out of range.')

          #call computeparameter(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
          #power=LinearPower_GEE_notime(mu, beta, gamma, rho0, II, JJ, KK, X_in)
          #dyn.load("computeparameter.so")
          mu=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$mu)
          beta=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$beta)
          gamma=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$gamma)
          rho0=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$rho0)
          #power=LinearPower_GEE(mu, beta, gamma, rho0, II, JJ, KK, X_in)
          #dyn.load("power_cal_wrapper.so")
          period=NULL
          period[1]=mu
          period[2]=beta
          X=matrix(1,JJ,2)

          for(i in 1:II)
          {
            X[,2]=X_in[i,]
            gmu=X%*%as.matrix(period)
            mmu=(gmu*(1-gmu))^0.5
            for(m in 1:JJ)
            {
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))>gmu[m]) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")# stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))>gmu(m)')
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))<max(2*gmu[m]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))<max(0,2*gmu(m)-1)')
            }
            for(m in 1:(JJ-1))
              for(k in (m+1):JJ)
              {
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")# stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
              }
          }
          power=suppressWarnings(.Fortran("LinearPowerGEEnotimewrapper",mu=mu, beta=beta, gamma=gamma, rho0=rho0, II=as.integer(II), JJ=as.integer(JJ), KK=as.integer(KK), X_in=as.integer(X_in),ICC0=ICC0,ICC2=ICC2,typeone=as.numeric(typeone), power=power )$power)
        }
        else if(link<2 & link>0)
        {
          p01=exp(mu)
          p11=exp(mu+beta)
          p0[1] = p01
          p11 = p11
          p0stepchange = p0totalchange/(JJ-1)
          for(j in 2:JJ)
          {
            gamma[j]= gamma[j-1] + p0stepchange
          }
          for(j in 2:JJ)
          {
            p0[j] = exp(gamma[j]+mu)
          }
          temp1=c(mu+p0totalchange,mu+beta+p0totalchange,mu,mu+beta)
          if(max(temp1)>0) stop('Violate theory of probability under log link: max(mu+gammaJ,mu+beta+gammaJ,mu,mu+beta)>0.Please check whether any of these four values are out of range.')

          #call computeparameterGEElog(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
          #power=LogPower_GEE_notime(mu, beta, gamma, rho0, II, JJ, KK, X_in)
          #dyn.load("computeparameterGEElog.so")
          mu=suppressWarnings(.Fortran("computeparameterGEElog",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$mu)
          beta=suppressWarnings(.Fortran("computeparameterGEElog",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$beta)
          gamma=suppressWarnings(.Fortran("computeparameterGEElog",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$gamma)
          rho0=suppressWarnings(.Fortran("computeparameterGEElog",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$rho0)
          #power=LogPower_GEE(mu, beta, gamma, rho0, II, JJ, KK, X_in)
          #dyn.load("power_cal_wrapper.so")
          period=NULL
          period[1]=mu
          period[2]=beta
          X=matrix(1,JJ,2)
          for(i in 1:II)
          {
            X[,2]=X_in[i,]
            gmu=X%*%as.matrix(period)
            gmu=exp(gmu)
            mmu=(gmu*(1-gmu))^0.5
            for(m in 1:JJ)
            {
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))>gmu[m]) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")# stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))>gmu(m)')
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))<max(2*gmu[m]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))<max(0,2*gmu(m)-1)')
            }
            for(m in 1:(JJ-1))
              for(k in (m+1):JJ)
              {
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0))  stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
              }
          }

          power=suppressWarnings(.Fortran("LogPowerGEEnotimewrapper",mu=mu, beta=beta, gamma=gamma, rho0=rho0, II=as.integer(II), JJ=as.integer(JJ), KK=as.integer(KK), X_in=as.integer(X_in),ICC0=ICC0,ICC2=ICC2,typeone=as.numeric(typeone),power=power )$power)
        }
        else
        {
          p01=exp(mu)/(1+exp(mu))
          p11=exp(mu+beta)/(1+exp(mu+beta))
          p0[1] = p01
          p11 = p11
          p0stepchange = p0totalchange/(JJ-1)
          for(j in 2:JJ)
          {
            gamma[j]= gamma[j-1] + p0stepchange
          }
          for(j in 2:JJ)
          {
            p0[j]= exp(gamma[j]+mu)/(1+exp(mu+gamma[j]))
          }
          #call computeparameterGEElogit(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
          #power=LogitPower_GEE_notime(mu, beta, gamma, rho0, II, JJ, KK, X_in)
          #dyn.load("computeparameterGEElogit.so")
          mu=suppressWarnings(.Fortran("computeparameterGEElogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$mu)
          beta=suppressWarnings(.Fortran("computeparameterGEElogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$beta)
          gamma=suppressWarnings(.Fortran("computeparameterGEElogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$gamma)
          rho0=suppressWarnings(.Fortran("computeparameterGEElogit",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$rho0)
          #dyn.load("power_cal_wrapper.so")
          period=NULL
          period[1]=mu
          period[2]=beta
          X=matrix(1,JJ,2)
          for(i in 1:II)
          {
            X[,2]=X_in[i,]
            gmu=X%*%as.matrix(period)
            gmu=exp(gmu)/(1+exp(gmu))
            mmu=(gmu*(1-gmu))^0.5
            for(m in 1:JJ)
            {
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))>gmu[m]) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))>gmu(m)')
              if((gmu[m]*gmu[m]+alpha0*gmu[m]*(1-gmu[m]))<max(2*gmu[m]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('(gmu(m)*gmu(m)+alpha0*gmu(m)*(1-gmu(m)))<max(0,2*gmu(m)-1)')
            }
            for(m in 1:(JJ-1))
              for(k in (m+1):JJ)
              {
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0))  stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha1*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))>min(c(gmu[m],gmu[k]))')
                if(((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)) stop("Correlation parameters do not satisfy the restrictions of Qaqish (2003). Please  check whether it is possible to reduce the effect size, or make adjustments to  the intraclass correlations.")#stop('((gmu[m]*gmu[k]+alpha2*(mmu[m])*(mmu[k])))<max(gmu[m]+gmu[k]-1,0)')
              }
          }
          power=suppressWarnings(.Fortran("LogitPowerGEEnotimewrapper",mu=mu, beta=beta, gamma=gamma, rho0=rho0, II=as.integer(II), JJ=as.integer(JJ), KK=as.integer(KK), X_in=as.integer(X_in),ICC0=ICC0,ICC2=ICC2,typeone=as.numeric(typeone),power=power )$power)
        }
      }
    }
  }
  else{
    alpha0= ICC0
    alpha1= rho0
    alpha2= ICC2
    lambda1 = 1- alpha0 + alpha1 - alpha2
    lambda2 = 1- alpha0 - (JJ-1)*(alpha1 - alpha2)
    lambda3 = 1 + (KK-1)*(alpha0 - alpha1) - alpha2
    lambda4 = 1 + (KK-1)*alpha0 + (JJ-1)*(KK-1)*alpha1 + (JJ-1)*alpha2
    temp2=c(lambda1,lambda2,lambda3,lambda4)
    if(min(temp2)<0) stop('Correlation matrix R is no longer positive definite. Please check whether the inter-period correlation is unrealistically larger than the within-period correlation or the within-individual correlation.')
    # if(ICC1>ICC0)  stop('ICC1>ICC0')


    #if(ICC1>ICC2)  stop('ICC1>ICC2')
    if (p0totalchange > TOLERANCE | p0totalchange < -TOLERANCE)
    {   p01=mu
    p11=mu+beta
    p0[1] = p01
    p11 = p11
    p0stepchange = p0totalchange/(JJ-1)
    for(j in 2:JJ)
    {
      p0[j]= p0[j-1] + p0stepchange
    }
    #  call computeparameter(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
    #dyn.load("computeparameter.so")
    mu=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$mu)
    beta=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$beta)
    gamma=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$gamma)
    tau2=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$tau2)
    power=suppressWarnings(.Fortran("ContinuousPowerGEEtimewrapper",mu=mu, beta=beta, gamma=gamma,rho0=rho0, II=as.integer(II), JJ=as.integer(JJ), KK=as.integer(KK), X_in=as.integer(X_in),sigma2=sigma2,ICC0=ICC0,ICC2=ICC2,typeone=as.numeric(typeone),power=power )$power)
    }
    else
    { p01=mu
    p11=mu+beta
    p0[1] = p01
    p11 = p11
    p0stepchange = p0totalchange/(JJ-1)
    for(j in 2:JJ)
    {
      p0[j]= p0[j-1] + p0stepchange
    }
    #  call computeparameter(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
    #dyn.load("computeparameter.so")
    mu=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$mu)
    beta=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$beta)
    gamma=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$gamma)
    tau2=suppressWarnings(.Fortran("computeparameter",JJ=as.integer(JJ), mu=mu, beta=beta, gamma=gamma, tau2=tau2, p0=p0, p11=p11, rho0=rho0 )$tau2)
    power=suppressWarnings(.Fortran("ContinuousPowerGEEnotimewrapper",mu=mu, beta=beta, gamma=gamma,rho0=rho0, II=as.integer(II), JJ=as.integer(JJ), KK=as.integer(KK), X_in=as.integer(X_in),sigma2=sigma2,ICC0=ICC0,ICC2=ICC2,typeone=as.numeric(typeone),power=power )$power)
    }
  }
  #summary

  if(type=="cohort") SS=II*KK
  else SS=II*KK*JJ

  if(response<1.5) {
    summary=matrix(data=c(I,J,K,SS,type,family,md,lf,round(beta,3),round(gamma[J],3),sigma2,ICC0,ICC1,ICC2,alpha,round(power,3)),ncol=1)
    rownames(summary)=c("I","J","K","total sample size","study type","family of outcomes","model","link","treatment effect beta","time effect gamma_J","marginal variance","alpha0","alpha1","alpha2","Type I error","Power")
    design_matrix=dataset
    list_data =list(design_matrix,summary)
    names(list_data)=c("design matrix dataset, row-cluster col-time","Summary")
    return(list_data)
  }

  if(response>1.5) {
    summary=matrix(data=c(I,J,K,SS,type,family,md,lf,round(mu,3),round(beta,3),round(gamma[J],3),ICC0,ICC1,ICC2,alpha,round(power,3)),ncol=1)
    rownames(summary)=c("I","J","K","total sample size","study type","family of outcomes","model","link","baseline effect mu","treatment effect beta","time effect gamma_J","alpha0","alpha1","alpha2","Type I error","Power")
    design_matrix=dataset
    list_data =list(design_matrix,summary)
    names(list_data)=c("design matrix dataset, row-cluster col-time","Summary")
    return(list_data)
  }

}

