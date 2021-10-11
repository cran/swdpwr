#' A function for calculating the intracluster correlation coefficient (ICC) for binary outcomes given the cluster level random effects variance
#'
#' @description This function calculates the ICC (intracluster correlation coeffcient which measures the correlation between individuals in the same cluster) under different link funcitions in a cross-sectional stepped wedge CRT with binary outcomes. This model considers only the fixed time effects and does not include cluster by time interaction random effect.
#'
#' @param link choose link function from link="identity", link="log" and link="logit", with default value of identity link
#' @param meanresponse_start the anticipated mean response in the control group at the start of the study
#' @param tau2 also denoted as sigma_b: variance of the between-cluster random effect, default is 0
#' @return The object returned includes the link function and value for the ICC in this study
#' @examples
#' BinICC(link="identity",meanresponse_start=0.2,tau2=0.05)
#'
#' @import spatstat.core
#'
#' @export
#'
BinICC<-function(link="identity", meanresponse_start, tau2=0){
  if(meanresponse_start>1| meanresponse_start<0) stop("meanresponse_start for binary outcomes should be between 0 and 1")
  if(link!="identity"&&link!="log"&&link!="logit") stop("Wrong input for link function names.")
  if(link=="identity") mu = meanresponse_start
  if(link=="log") mu = log(meanresponse_start)
  if(link=="logit") mu = log(1/(1-meanresponse_start)-1)

  if(link=="identity") ICC=tau2/(mu*(1-mu))
  if(link=="log") ICC=(exp(2*mu+2*tau2)-exp(2*mu+tau2))/(exp(mu+tau2/2)-exp(2*mu+tau2))
  if(link=="logit")
    {temp = gauss.hermite(function(x) (exp(mu+x)/(1+exp(mu+x)))^2, 0, tau2)
   ICC=(temp- meanresponse_start^2)/(temp- meanresponse_start^2+meanresponse_start-temp)
  }
  if(ICC>1|ICC<0) stop("ICC should be between 0 and 1, please modify tau2 or meanresponse_start.")
  list_data=list(study_type="cross-sectional",link_function=link ,ICC=ICC)
  return(list_data)
}

