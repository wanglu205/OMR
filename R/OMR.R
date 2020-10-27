########################################################################################################################
# Package: OMR
# Version: 1.0
# Date   : 2020-10-21
# Title  : Probabilistic Mendelian Randomization under the Omnigenic Genetic Architecture
# Authors: Lu Wang, Boran Gao, Yue Fun, and Xiang Zhou
# Contact: willa0205@yeah.net
#          University of Michigan, Department of Biostatistics
########################################################################################################################

#' Probabilistic Mendelian Randomization under the Omnigenic Genetic Architecture
#' 
#' Two-sampling Mendelian randomization analysis under omnigenic genetic architecture assumption.
#' 
#' OMR makes an omnigenic modeling assumption on the SNP effects on the exposure variable and uses genome-wide SNPs as instrumental variables without pre-selection to enable more powerful MR method.
#'   
#' 
#' @param nx an integer specifying the sample size in exposure study.
#' @param ny an integer specifying the sample size in outcome study.
#' @param n_ref an integer specifying the sample size in reference study.
#' @param l.j a vector representing the LD Scores calculated from reference study.
#' @param Z a matrix containing marginal z-scores, with the first coloum containing z-scores for outcome study and the second coloum containing z-scores for exposure study.
#' @param numCore a positive integer specifying the number of cores for parallel computing (default = 1)
#' @param two_study whether the study is a two-sample setting, the defult value is TRUE.
#' @return A list that contains
#' @return \item{alpha}{Estimate of causal effect. Assuming the summary statistics are standardized, \code{alpha} represents increase in mean value of Y in s.d. unit of Y (for continuous outcomes) or log-OR of Y (for binary outcomes) associated with per s.d. unit increase in values of X (for continuous exposures) or values of X changing from 0 to 1 (for binary exposures).}
#' @return \item{beta}{Proportion of SNP heritability in the exposure variable}
#' @return \item{gamma}{Proportion of SNP horizontal pleiotropy effect in the outcome variable}
#' @return \item{z_stat}{Z-statistic for test of the causal effect estimate.}
#' @return \item{se_beta}{Standard error of causal effect estimate.}
#' @return \item{pvalue}{P value for the causal effect, based on the z test.}
#' @author Lu Wang, Boran Gao, Yue Fun, and Xiang Zhou
#' @examples 
#' data(exampledata)
#' attach(exampledata)
#' res=OMR(ny,nx,num.per,l.j,Z)
#' closeAllConnections()
#' detach(exampledata)


############################################
#
#         OMR Function
#
#############################################
OMR <- function(ny,nx,num.per,l.j,Z,numCore=1,two_study=T){
  registerDoParallel(cores=numCore)
  num.snp<-length(l.j)
  wind_size<-round(num.snp/100)
  seq_1<-c()
  seq_2<-c()
  for(i in 1:100){
    seq_1[i] = wind_size*(i-1)+1
    seq_2[i] = wind_size*i
  }
  seq_2[100] = num.snp
  interval_matrix<-as.matrix(rbind(seq_1,seq_2))
  tmp <- foreach(iVar=1:100, .combine=rbind, .errorhandling = 'remove')%dopar%{
    tmpt <- data.frame()
    jack_lj<-l.j[-c(interval_matrix[1,iVar]:interval_matrix[2,iVar])]
    jack_z<-Z[-c(interval_matrix[1,iVar]:interval_matrix[2,iVar]),]
    iter <- Iteration(n1,n2,num.per,jack_lj,jack_z)
    tmpt=rbind(tmpt,data.frame(alpha_loop=iter[1],sigma_beta_loop=iter[2],sigma_gamma_loop=iter[3],sigma_inverse_loop=iter[5],LogL_new=iter[8]))
    return(tmpt)
  }
  res <- list(alpha = mean(tmp[,1]),
              beta = num.snp*mean(tmp[,2]),
			  gamma  = num.snp*mean(tmp[,3]))
  res$se = sqrt(sum((tmp[,1]-mean(tmp[,1]))^2)*0.99) # Standard error
  res$z_stat = res$alpha/res$se
  res$pvalue = 2*(pnorm(abs(res$z_stat),lower.tail = FALSE))
  return(res)
}

#########################
#
#     Iterations 
#
#########################

Iteration<-function(ny,nx,num.per,l.j,Z,two_study=T){
  num.snp<-length(l.j)
  z_1j <- Z[,1]^2
  z_2j <- Z[,2]^2
  z_12 <- Z[,1]*Z[,2]
  x1<-l.j*nx
  lmod <- lm(z_2j~x1)
  sigma_beta_loop <- coef(lmod)[2]
  x2<-l.j*sqrt(ny)*sqrt(nx)
  Linear_mod_2<-lm(z_12~x2)
  rho_0<-coef(Linear_mod_2)[1]
  if(two_study==T | abs(rho_0) < 0.1){
    rho <- 0
    sigma_inverse_loop <- matrix(c(1,rho,rho,1),ncol = 2,nrow=2)
    alpha_loop<-mean(z_12)/(mean(l.j)*sqrt(ny)*sqrt(nx)+num.snp)/sigma_beta_loop  ###two study
  }else{
    rho<-coef(Linear_mod_2)[1]
    alpha_loop <- (mean(z_12)- rho_0)/(mean(l.j)*sqrt(ny)*sqrt(nx)+num.snp)/sigma_beta_loop
    sigma_inverse_loop <- matrix(c(1,rho,rho,1),ncol = 2,nrow=2)
  }
  x3<-l.j*ny+num.snp
  lmod <- lm(z_1j~x3)
  sigma_gamma_loop <- coef(lmod)[2]-alpha_loop*alpha_loop*sigma_beta_loop
  iter=100
  num_iter=1
  out=EM(ny,nx,num.snp,num.per,sigma_inverse_loop,alpha_loop,sigma_beta_loop,sigma_gamma_loop,l.j,Z)
  alpha_loop = out[1]
  sigma_beta_loop = out[2]
  if(out[3] < 0){sigma_gamma_loop = 0}else{sigma_gamma_loop = out[3]}
  LogL_new = out[4]
  num_iter= num_iter+1 
  out=EM(ny,nx,num.snp,num.per,sigma_inverse_loop,alpha_loop,sigma_beta_loop,sigma_gamma_loop,l.j,Z)
  alpha_loop = out[1]
  sigma_beta_loop = out[2]
  if(out[3] < 0){sigma_gamma_loop = 0}else{sigma_gamma_loop = out[3]}
  LogL_new = out[4]
  num_iter= num_iter+1 
  LogL_old=LogL_new
  while(num_iter<iter){
    out=EM(ny,nx,num.snp,num.per,sigma_inverse_loop,alpha_loop,sigma_beta_loop,sigma_gamma_loop,l.j,Z)
    alpha_loop = out[1]
    sigma_beta_loop = out[2]
    if(out[3] < 0){sigma_gamma_loop = 0}else{sigma_gamma_loop = out[3]}
    LogL_new = out[4]
    flag <- LogL_new-LogL_old
    if(flag <1e-03| is.na(flag)){break}
    LogL_old<- LogL_new
    num_iter= num_iter+1 
  }
  num.iter=0
  iter=10
  flag=0
  LogL_old=LogL_new
  while(num.iter<iter){
    NR_out<-NR(ny,nx,num.snp,num.per,sigma_inverse_loop,alpha_loop,sigma_beta_loop,sigma_gamma_loop,l.j,Z)
    LogL_new = NR_out[4]
    alpha_loop<-NR_out[1]
    sigma_beta_loop<-NR_out[2]
    if(NR_out[3] < 0){sigma_gamma_loop = 0}else{sigma_gamma_loop = NR_out[3]}
    flag <- LogL_new-LogL_old
    if(flag <0 | is.na(flag)){break}
    LogL_old<- LogL_new
    num.iter<-num.iter+1
    if(flag <1e-05){break}
  }
  return(c(alpha_loop,sigma_beta_loop,sigma_gamma_loop,sigma_inverse_loop,LogL_new,flag))
}

# end function
#########################################
#             CODE END                  #
#########################################

