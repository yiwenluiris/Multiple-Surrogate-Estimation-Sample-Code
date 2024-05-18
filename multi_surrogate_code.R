library(matlib)
library(MASS)
library(tidyverse)

expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(p/(1-p))}

####################################################################################
############################# Proposed method ######################################
####################################################################################
# y: True phenotypes (with NA)
# !!!!!! other than validation set, the value of y should be "NA" !!!!!!!!!!!!!!!!

# s: Multiple algorithm-generated phenotypes (surrogates)

# x: Covariates (characteristics of patients)
# !!!!!! categorical predictors have to be transferred into dummy variables !!!!!!
####################################################################################

augmented_est <- function(y, s, x){
  
  id <- which(y != 'NA')
  x <- cbind(1, x)
  
  val_x = x[id,]          # validation set of x
  val_y = y[id]           # validation set of y
  val_s = s[id,]           # validation set of s
  
  nr <- dim(x)[2] #Number of variables
  n <- dim(x)[1] #Total sample size
  m <- length(val_y) #number of validated individuals
  rho <- m/n #Validation set ratio
  
  
  # glm for validation set only with true disease
  glm1 <- glm(val_y~val_x[,-1], family="binomial")
  beta_hat <- glm1$coefficients
  
  var_beta_hat <- summary(glm1)$coefficients[,2]**2
  
  gamma_hat<-vector()
  var_gamma_hat<-vector()
  gamma_bar<-vector()
  var_gamma_bar<-vector()
  
  
  for (ind_s in 1: ncol(s)){
    #glm for validation set validation set only with surrogate
    glm2 <- glm(val_s[,ind_s]~val_x[,-1],family="binomial")
    gamma_hat_temp <- glm2$coefficients
    var_gamma_hat_temp <- summary(glm2)$coefficients[,2]**2
    gamma_hat = rbind(gamma_hat, gamma_hat_temp)
    var_gamma_hat = rbind(var_gamma_hat, var_gamma_hat_temp)
    
    #glm for all surrogate
    glm3 <- glm(s[,ind_s]~x[,-1],family="binomial")
    gamma_bar_temp <- glm3$coefficients
    var_gamma_bar_temp <- summary(glm3)$coefficients[,2]**2
    gamma_bar = rbind(gamma_bar, gamma_bar_temp)
    var_gamma_bar = rbind(var_gamma_bar, var_gamma_bar_temp)
  }
  
  ####################################################
  # to get the augmented estimator 
  Z1 <- beta_hat %*% t(x)
  Z2 <- gamma_hat %*% t(x)
  Z3 <- gamma_bar %*% t(x)
  
  val_Z1 <- Z1[id]
  val_Z2 <- Z2[,id]
  val_Z3 <- Z3[,id]
  
  info_matrix1 <- matrix(0 , nrow = nr, ncol = nr)     # information matrix
  for (i in 1:nr){
    for (j in 1:nr){
      if (i == j){
        info_matrix1[i,j] <- sum( val_x[,i]^2 * exp(val_Z1) / (1+exp(val_Z1))^2 ) /m
      }else
        info_matrix1[i,j] <- sum( val_x[,i] * val_x[,j] * exp(val_Z1) / (1+exp(val_Z1))^2) /m
      info_matrix1[j,i] <- info_matrix1[i,j]
    }
  }
  
  #dim of info_matrix1: nr nr
  
  info_matrix2 <- matrix(0 , nrow = nr*ncol(s), ncol = nr*ncol(s))     # information matrix
  
  for(k in 1:ncol(s)){
    for (i in 1:nr){
      for (j in 1:nr){
        if (i == j){
          info_matrix2[i+(k-1)*nr,j+(k-1)*nr] <- rowSums(x[,i]^2 * exp(Z3) / (1+exp(Z3))^2 )[k] /n
        }else {
          info_matrix2[i+(k-1)*nr,j+(k-1)*nr] <- rowSums(x[,i] * x[,j] * exp(Z3) / (1+exp(Z3))^2)[k] /n
          info_matrix2[j+(k-1)*nr,i+(k-1)*nr] <- info_matrix2[i+(k-1)*nr,j+(k-1)*nr]
        }
      }
    }
  }
  
  p1 <- matrix(0, nrow = nr, ncol =nr)
  p2 <- matrix(0, nrow = nr*ncol(s), ncol =nr*ncol(s))
  p3 <- matrix(0, nrow = nr*ncol(s), ncol =nr)
  
  for (i in 1:m){
    
    phi1 <- (val_y[i] - (exp(val_Z1[i])/(1+exp(val_Z1[i])))) %*% as.numeric(val_x[i,])
    phi1 <- inv(info_matrix1) %*% as.vector(phi1)
    p1 <- p1 + phi1 %*% t(phi1)
    
    phi2 <- (val_s[i,] - (exp(val_Z3[,i])/(1+exp(val_Z3[,i])))) %*% t(as.numeric(val_x[i,]))
    phi2 <- inv(info_matrix2) %*% as.vector(phi2)
    p2 <- p2 + as.vector(phi2) %*% t(as.vector(phi2))
    
    p3 <- p3 + phi2 %*% t(phi1)
  }
  
  sigma <- p1 / m 
  sigma_star <- (1-rho) * p2 / m 
  omega <- (1-rho) * p3 / m 
  
  beta_aug_mul = t(beta_hat - t(omega) %*% inv(sigma_star) %*% as.vector(gamma_hat - gamma_bar))
  beta_aug_mul = as.vector(beta_aug_mul)
  var_aug_mul = (sigma - t(omega) %*% inv(sigma_star) %*% omega) /m
  var_aug_mul = diag(var_aug_mul)
  ci_aug_mul <- c(t(c(beta_aug_mul) + qnorm(0.975)*sqrt(var_aug_mul)%*%t(c(-1,1))))
  
  return(list(beta_aug_mul, var_aug_mul, ci_aug_mul))
}


## ===================================================================================== ##
## ===================================== Example ======================================= ##
## ===================================================================================== ##

sample_data <- read.csv("sample_data.csv")
y <- as.matrix(sample_data[,1])
s <- as.matrix(sample_data[,c(2:4)])
x <- as.matrix(sample_data[,c(5:6)])
  
augmented_result <- augmented_est(y, s, x)
