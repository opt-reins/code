##########################################################################
####  Optimal Reinsurance Multivariate: Simulation Study 1            ####
##########################################################################

## ------------------ setting ------------------------
library(copula)
library(Matrix)
library(nleqslv)
library(ggplot2)
library(patchwork)
B      <- 50000        # number of replicates
n      <- 100          # policies per replicate
## Marginal parameters
shape1 <- 5            # Pareto-II shape  (X_{t,1})
scale1 <- 4            # Pareto-II scale
rate2  <- 1            # exponential rate (X_{t,2})
## t-copula parameters
nu     <- 4            # degrees of freedom
rho    <- 0.5          # correlation
## Loading factors
rho_j     <- c(0.5, 0.5)
rho_tilde <- c(1.0, 1.0)
mu0 <- 0
## ---------------------------------------------------


##### functions #####
qParetoII <- function(u, scale = 1, shape = 2)  scale * ((1-u)^(-1/shape) - 1)
qExp      <- function(u, rate = 1)  -log(1 - u)/rate

sim_X <- function(n) {
  tc <- tCopula(rho, dim = 2, df = nu, dispstr = "un")
  U  <- rCopula(n, tc)     # n×2 matrix of uniforms
  X1 <- qParetoII(U[,1], scale = scale1, shape = shape1)
  X2 <- qExp     (U[,2], rate = rate2)
  cbind(X1, X2)    
}

delta  <-  rho_j      / sqrt(n)
deltaT <-  rho_tilde  / sqrt(n)

## E[I_j(X)] for each treaty
mean.eval <- function(cpar2, cpar1, treaty) { 
  mu_j.eval <- c(scale1/(shape1 - 1), 1/rate2)
  cpar <- c(cpar1,cpar2)
  if (treaty == "quota") {EI = cpar*mu_j.eval} else {
    EI = c(ifelse(cpar[1] < 0, mu_j.eval[1],
             scale1/(shape1 - 1) * (1 + cpar[1]/scale1)^(-(shape1 - 1))),
      exp(-rate2*cpar[2])/rate2)}
  sum(rho_tilde*EI)-sum(rho_j*mu_j.eval)-mu0
}

## Simulation of total cost
simulate_T <- function(treaty = c("quota", "stop"), seed, cpar1, cpar2.init, alpha.eval) {
  cpar2 <- nleqslv(x = cpar2.init, fn = mean.eval, method = "Broyden", 
                    cpar1=cpar1, treaty=treaty)$x
  cpar <- c(cpar1,cpar2)
  set.seed(seed)
  X <- sim_X(n * B)
  grp <- rep(1:B, each = n)
  if (treaty == "quota") {
    retained <- X - sweep(X, 2, cpar, `*`)
  } else { # stop-loss
    retained <- pmin(X, matrix(cpar, nrow(X), 2, byrow = TRUE))
  }
  ## compute total cost
  term1 <- rowsum(retained, grp)         
  EI <- apply(X-retained,2,mean)
  mu_j <- apply(X,2,mean)
  premR <- n * (1 + deltaT) * EI      
  premP <- n * (1 + delta ) * mu_j    
  
  Tvec <- rowSums(term1) + sum(premR)  - sum(premP)  
  
  VaR.true <- as.numeric(quantile(as.numeric(Tvec),probs=alpha.eval))
  VaR.approx <- sqrt(n)*sd(rowSums(retained))*qnorm(alpha.eval)+n*(sum(deltaT*EI)-sum(delta*mu_j))
  
  list(VaR.true=VaR.true, VaR.approx=VaR.approx, cpar=cpar)
}

##### computing VaR vs variables #####
length.out=51
c1_seq = seq(from=0, to=1, length.out=length.out)
d1_seq = seq(from=0.1, to=2.5, length.out=length.out)
alpha.eval = seq(from=0.5, to=1, length.out=length.out)
VaR.quota.true.mat=VaR.quota.approx.mat=array(NA,dim=c(length.out,length.out))
VaR.stop.true.mat=VaR.stop.approx.mat=array(NA,dim=c(length.out,length.out))
VaR.quota.cpar.mat=VaR.stop.cpar.mat=array(NA,dim=c(length.out,2))
for (s in 1:length.out) {
  T_quota <- simulate_T("quota", seed=2025, cpar1=c1_seq[s], cpar2.init=0.5, alpha.eval=alpha.eval)
  VaR.quota.true.mat[s,]=T_quota$VaR.true
  VaR.quota.approx.mat[s,]=T_quota$VaR.approx
  VaR.quota.cpar.mat[s,]=T_quota$cpar
}
for (s in 1:length.out) {
  T_stop <- simulate_T("stop", seed=2025, cpar1=d1_seq[s], cpar2.init=1, alpha.eval=alpha.eval)
  VaR.stop.true.mat[s,]=T_stop$VaR.true
  VaR.stop.approx.mat[s,]=T_stop$VaR.approx
  VaR.stop.cpar.mat[s,]=T_stop$cpar
}

## vectors for true vs approx VaR (alpha=0.9) vs loading factors/retentions/alpha
VaR.quota.true_c=VaR.quota.true.mat[,41]
VaR.quota.approx_c=VaR.quota.approx.mat[,41]
VaR.quota.true_alpha=VaR.quota.true.mat[26,-51] #remove the last one
VaR.quota.approx_alpha=VaR.quota.approx.mat[26,-51]
VaR.stop.true_c=VaR.stop.true.mat[,41]
VaR.stop.approx_c=VaR.stop.approx.mat[,41]
VaR.stop.true_alpha=VaR.stop.true.mat[26,-51] #remove the last one
VaR.stop.approx_alpha=VaR.stop.approx.mat[26,-51]


##### Plot the results computed above using ggplot2 #####
## plot function ##
build_plot <- function(df, xlab, title) {
  ggplot(df, aes(x = x, y = value, linetype = type, colour = type)) +
    geom_line(size = 1) +
    scale_colour_manual(values = c("approx" = "black",
                                   "true"   = "black"),
                        labels = c("Approx. VaR", "True VaR")) +
    scale_linetype_manual(values = c("approx" = "solid",
                                     "true"   = "dotted")) +
    guides( linetype = "none",                                 
            colour = guide_legend(
              override.aes = list(linetype = c("solid", "dotted"),  
                                  size     = 1.2),               
              keywidth = unit(2, "cm"))) +
    labs(x = xlab, y = "VaR of total cost", colour = "", linetype = "",
         title = title) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 14, hjust = .5))
}

## VaR vs c1 (quota)
df_qc <- rbind(
  data.frame(x = c1_seq, value = VaR.quota.approx_c, type = "approx"),
  data.frame(x = c1_seq, value = VaR.quota.true_c,  type = "true")
)
p_qc <- build_plot(df_qc, expression(c[1]),
                   "Value-at-Risk vs ceded share (Quota-share, n = 100, α = 0.9)")
p_qc
## VaR vs d1 (stop-loss)
df_sd <- rbind(
  data.frame(x = d1_seq, value = VaR.stop.approx_c, type = "approx"),
  data.frame(x = d1_seq, value = VaR.stop.true_c,  type = "true")
)
p_sd <- build_plot(df_sd, expression(d[1]),
                   "Value-at-Risk vs retention level (Stop-loss, n = 100, α = 0.9)")
p_sd
## VaR vs alpha (quota) 
alph <- alpha.eval[-51] 
df_qa <- rbind(
  data.frame(x = alph, value = VaR.quota.approx_alpha, type = "approx"),
  data.frame(x = alph, value = VaR.quota.true_alpha,  type = "true")
)
p_qa <- build_plot(df_qa, expression(alpha),
                   expression(paste("Value-at-Risk vs α (Quota-share, n = 100, ", c[1]," = 0.5)")))
p_qa
## VaR vs alpha (stop-loss)
df_sa <- rbind(
  data.frame(x = alph, value = VaR.stop.approx_alpha, type = "approx"),
  data.frame(x = alph, value = VaR.stop.true_alpha,  type = "true")
)
p_sa <- build_plot(df_sa, expression(alpha),
                   expression(paste("Value-at-Risk vs α (Stop-loss, n = 100, ", d[1]," = 1.3)")))
p_sa