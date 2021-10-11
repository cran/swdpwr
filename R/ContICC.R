#' A function for calculating ICCs for continuous outcomes given random effects variances
#'
#' @description This function calculates the within-period, between-period, and within-individual correlation parameters for continuous outcomes in a stepped wedge CRT
#'
#' @param type choose the study type, specify type="cohort" for closed cohort study and type="cross-sectional" for cross-sectional study, default is "cross-sectional"
#' @param sigma2 marginal variance of the outcom, default is 1
#' @param sigma_b variance of the between-cluster random effect, default is 0
#' @param sigma_c variance of the cluster-by-time interaction random effect, default is 0
#' @param sigma_pi variance of the random effect for repeated measures of one individual, this parameter should not be specified for cross-sectional studies, default is NA
#' @return The object returned includes the study type and values for the ICCs in this study
#' @examples
#' ContICC(type="cohort",sigma2=1.5,sigma_b=0.5,sigma_c=0.2,sigma_pi=0.3)
#' @export
#'
#'
ContICC<-function(type="cross-sectional", sigma2=1, sigma_b=0, sigma_c=0, sigma_pi=NA){
  if(type!="cohort"&&type!="cross-sectional") stop("Wrong input for study type.")
  if(type=="cross-sectional"&&is.na(sigma_pi)!=1) stop("sigma_pi should not be specified for a cross-sectional study.")
  if(type=="cross-sectional") sigma_pi=0
  if(sum(sigma_b+sigma_c+sigma_pi)>sigma2) stop("The marginal variance should be larger than the sum of its components.")
  alpha0= (sigma_b+sigma_c)/sigma2
  alpha1= (sigma_b)/sigma2
  alpha2= (sigma_b+sigma_pi)/sigma2
  if(type=="cohort") {list_data=list(study_type=type,alpha0=alpha0,alpha1=alpha1,alpha2=alpha2)
  return(list_data)}
  if(type=="cross-sectional") {list_data=list(study_type=type,alpha0=alpha0,alpha1=alpha1)
  return(list_data)}
  }

