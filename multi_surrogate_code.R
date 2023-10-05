library(matlib)
library(MASS)

expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(p/(1-p))}

####################################################################################
############################# Proposed method ######################################
####################################################################################
# n: total size of data set
# m: size of valiadation set
# x: Covatiates (characteristics of patients)
# val_x, val_s, val_y: validation part of covariates, surrogates, and gold outcome
# beta_hat: estimation using val_y and val_x
# gamma_hat: estimation using val_s and val_x
# gamma_bar: estimation using s and x
####################################################################################

augmented_est <- function(n, m, x, val_x, val_s, val_y, beta_hat, gamma_hat, gamma_bar){
  rho <- m/n #validation set ratio
  nr <- dim(x)[2]
  
  Z1 <- beta_hat %*% t(x)
  Z2 <- gamma_hat %*% t(x)
  Z3 <- gamma_bar %*% t(x)
  
  val_Z1 <- Z1[id]
  val_Z2 <- Z2[,id]
  val_Z3 <- Z3[,id]
  
  # information matrix
  info_matrix1 <- matrix(0 , nrow = nr, ncol = nr)     
  for (i in 1:nr){
    for (j in 1:nr){
      if (i == j){
        info_matrix1[i,j] <- sum( val_x[,i]^2 * exp(val_Z1) / (1+exp(val_Z1))^2 ) /m
      }else
        info_matrix1[i,j] <- sum( val_x[,i] * val_x[,j] * exp(val_Z1) / (1+exp(val_Z1))^2) /m
      info_matrix1[j,i] <- info_matrix1[i,j]
    }
  }
  
  info_matrix2 <- matrix(0 , nrow = nr*ncol(s), ncol = nr*ncol(s))
  
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