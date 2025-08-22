##########################################################################
####  Optimal Reinsurance Multivariate: Simulation Study 2 (Case 1)   ####
##########################################################################

## ------------------------- 1.  packages --------------------------
library(nleqslv)     
library(copula)      
library(Matrix)      
library(parallel)    
library(stats)       

## -------------------- 2.  case setting ---------------------
B        <- 500      # outer simulations
Q        <- 500      # bootstrap resamples inside each simulation
n_cores  <- max(1L, detectCores() - 1L)
set.seed(20250604)   

n         <- 5000 # sample size: 500, 2000 or 5000
m         <- 2
h_vec     <- c(0.01*n^(-1/3), 0.01*n^(-1/3))

shape1    <- 5      # Pareto-II   (X_{t,1})
scale1    <- 4
rate2     <- 1      # Exponential (X_{t,2})

phi       <- 0.5    # t-copula
nu        <- 4
mu0       <- 0
rho       <- c(0.5, 0.5)
rho_tilde <- c(1.0, 1.0)
alpha     <- 0.9
## --------------------------------------------------------------------------- ##


## -------------------- 3.  functions ---------------
## 3.1  random losses  X  ~  t-copula + Pareto-II + Exp
qParetoII <- function(u, scale = 1, shape = 2) scale * ((1 - u)^(-1/shape) - 1)
qExp      <- function(u, rate  = 1)            -log(1 - u) / rate

sim_X <- function(n) {
  tc <- tCopula(phi, dim = 2, df = nu, dispstr = "un")
  U  <- rCopula(n, tc)                 # n × 2 uniforms
  cbind(qParetoII(U[, 1], scale1, shape1),
        qExp     (U[, 2], rate2))
}

## 3.2  quota-share estimator --------------------------------------------------
est_QS <- function(X) {
  n <- nrow(X); m <- ncol(X)
  
  mu_hat    <- colMeans(X)
  X_center  <- sweep(X, 2, mu_hat, FUN = "-")
  Sigma_hat <- crossprod(X_center) / n
  
  numerator   <- sum(mu_hat * (rho_tilde - rho)) - mu0
  denominator <- as.numeric(t(mu_hat * rho_tilde) %*% solve(Sigma_hat) %*% (mu_hat * rho_tilde))
  lambda_hat  <- 2 * numerator / denominator
  
  c_hat <- 1 - (lambda_hat / 2) * (solve(Sigma_hat) %*% (mu_hat * rho_tilde))
  c_hat <- as.numeric(c_hat)
  
  H1_sq <- as.numeric(t(1 - c_hat) %*% Sigma_hat %*% (1 - c_hat))
  H1    <- sqrt(H1_sq)
  H2    <- sum(c_hat * mu_hat * rho_tilde) - sum(rho * mu_hat)
  VaR   <- sqrt(n) * (H1 * qnorm(alpha) + H2)
  
  c(c_hat, VaR_qs = VaR)
}

## 3.3  stop-loss estimator ----------------------------------------------------
est_SL <- function(X) {
  n <- nrow(X); m <- ncol(X)
  mu_hat <- colMeans(X)
  
  compute_L <- function(par) {
    d_vec  <- par[1:m]; lambda <- par[m + 1]
    
    scaled  <- -(sweep(X, 2, d_vec)) / matrix(h_vec, n, m, byrow = TRUE)
    K_mat   <- pnorm(scaled)                # n × m
    barK    <- 1 - K_mat
    
    X_wedge <- X * K_mat + matrix(d_vec, n, m, byrow = TRUE) * barK
    retained      <- rowSums(X_wedge)
    mean_retained <- mean(retained)
    mean_barK     <- colMeans(barK)
    
    A_vec <- (2 / n) * as.numeric(crossprod(retained, barK))  
    B_vec <- 2 * mean_retained * mean_barK
    C_vec <- lambda * rho_tilde * mean_barK
    L1    <- A_vec - B_vec - C_vec
    
    Delta   <- (X - matrix(d_vec, n, m, byrow = TRUE)) * barK
    term5   <- mean(rowSums(sweep(Delta, 2, rho_tilde, "*")))
    L2      <- term5 - sum(rho * mu_hat) - mu0
    
    c(L1, L2)
  }
  
  init      <- c(apply(X, 2, median), 1)
  sol       <- tryCatch(nleqslv(init, compute_L, method = "Broyden"), error = identity)
  if (inherits(sol, "error") || sol$termcd != 1) {
    return(rep(NA_real_, 1 + m + 1))  
  }
  d_hat <- sol$x[1:m]
  
  X_wdg <- pmin(X, matrix(d_hat, n, m, byrow = TRUE))
  s_wdg <- rowSums(X_wdg)
  H1_sq <- mean(s_wdg^2) - mean(s_wdg)^2
  H1    <- sqrt(H1_sq)
  Delta_plus <- pmax(X - matrix(d_hat, n, m, byrow = TRUE), 0)
  H2    <- sum(rho_tilde * colMeans(Delta_plus)) - sum(rho * mu_hat)
  VaR   <- sqrt(n) * (H1 * qnorm(alpha) + H2)
  
  c(d_hat, VaR_sl = VaR)
}

## 3.4  full point estimate + bootstrap SDs ---------------------------
one_run <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  ## --- point estimate (simulated X) -----------------------------------------
  X      <- sim_X(n)
  qs_hat <- est_QS(X)                
  sl_hat <- est_SL(X)                  
  
  ## --- bootstrap SDs --------------------------------------------------------
  boot_stat <- matrix(NA_real_, Q, 6)   
  
  for (q in 1:Q) {
    Xb <- X[sample.int(n, n, TRUE), ]   
    qs_b <- est_QS(Xb)
    sl_b <- est_SL(Xb)
    boot_stat[q, ] <- c(qs_b[1:2], sl_b[1:2], qs_b[3], sl_b[3])
  }
  sd_vec <- apply(boot_stat, 2, sd, na.rm = TRUE)
  names(sd_vec) <- c("sd_c1","sd_c2","sd_d1","sd_d2","sd_VaR_qs","sd_VaR_sl")
  
  list(c_hat   = qs_hat[1:2],
       d_hat   = sl_hat[1:2],
       VaR_qs  = qs_hat[3],
       VaR_sl  = sl_hat[3],
       sd_boot = sd_vec)
}
## --------------------------------------------------------------------------- ##


## -------------------- 4.  outer Monte-Carlo approximation ------------------
cl <- makeCluster(n_cores)
clusterExport(cl,
              varlist = c("n","m","h_vec","shape1","scale1","rate2","phi","nu",
                          "rho","rho_tilde","mu0","alpha","Q",
                          "qParetoII","qExp","sim_X",
                          "est_QS","est_SL","one_run"),
              envir = environment())
clusterEvalQ(cl, {
  library(copula)
  library(nleqslv)
  library(Matrix)})
mc_results <- parLapply(cl, 1:B, function(b) one_run(seed = 1e6 + b))
stopCluster(cl)
## --------------------------------------------------------------------------- ##


## -------------------- 5.  consolidate outputs --------------------------------
c_mat  <- t(sapply(mc_results, function(x) x$c_hat))           # B × 2
d_mat  <- t(sapply(mc_results, function(x) x$d_hat))           # B × 2
VaR_qs_vec <- sapply(mc_results, function(x) x$VaR_qs)         # length B
VaR_sl_vec <- sapply(mc_results, function(x) x$VaR_sl)         # length B

boot_sd_mat <- t(sapply(mc_results, function(x) x$sd_boot))    # B × 6
colnames(c_mat)       <- c("c1","c2")
colnames(d_mat)       <- c("d1","d2")
colnames(boot_sd_mat) <- c("sd_c1","sd_c2","sd_d1","sd_d2","sd_VaR_qs","sd_VaR_sl")

## -------------------- 6.  summary -----------------------
summary_df <- data.frame(
  mean_c1      = mean(c_mat[,1], na.rm = TRUE),
  mean_c2      = mean(c_mat[,2], na.rm = TRUE),
  mean_d1      = mean(d_mat[,1], na.rm = TRUE),
  mean_d2      = mean(d_mat[,2], na.rm = TRUE),
  mean_VaR_qs  = mean(VaR_qs_vec, na.rm = TRUE),
  mean_VaR_sl  = mean(VaR_sl_vec, na.rm = TRUE)
)
print(summary_df)
c(mean(VaR_qs_vec-VaR_sl_vec)-1.96*sd(VaR_qs_vec-VaR_sl_vec),
  mean(VaR_qs_vec-VaR_sl_vec)+1.96*sd(VaR_qs_vec-VaR_sl_vec))
## actual vs estimated 
print(apply(cbind(c_mat,d_mat,VaR_qs_vec,VaR_sl_vec),2,sd))
print(apply(boot_sd_mat,2,mean))

## Output objects:
##   c_mat           : 500 × 2  matrix of quota-share  c-vectors
##   d_mat           : 500 × 2  matrix of stop-loss  d-vectors
##   VaR_qs_vec      : length-500 vector of approximate VaR (quota-share)
##   VaR_sl_vec      : length-500 vector of approximate VaR (stop-loss)
##   boot_sd_mat     : 500 × 6  matrix of bootstrap SDs (in column order)
###############################################################################
save(c_mat, d_mat, VaR_qs_vec, VaR_sl_vec, boot_sd_mat, file="sim2a_3.Rda")
