#' Print the results of swdpower
#'
#' The \code{\link{print}} method for class "swdpower"
#'
#' @param x an object used to select a method.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The output from \code{\link{print}}
#' @export
#'
#'
print.swdpower<-function(x,...){
study<-paste(x$study.type,collapse = ', ')
sampesize<-paste(x$total.sample.size)
cat('This',study, "study has total sample size of", sampesize, "\n")
power<-paste(x$Power,collapse = ', ')
beta<-paste(x$treatment.effect.beta,collapse = ',')
error<-paste(x$Type.I.error,collapse = ',')
cat('Power for this scenario is',power, "for the alternative hypothesis treatment effect beta =",beta, "(two-sided Type I error = ",error, ") \n")
}
