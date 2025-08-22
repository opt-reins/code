##########################################################################
####  Optimal Reinsurance Multivariate: Data Analysis (Danish Fire Loss)##
##########################################################################

##### Load packages #####
library(Matrix)
library(nleqslv)
library(parallel)
library(stats)
library(dplyr)
library(ggplot2)

##### Retrieve Danish Fire Loss data #####
load("danishmulti.Rda") # if one has already retrieved from CASdatasets package
# Otherwise, run the following instead:
# library(CASdatasets); data(danishmulti)
X <- danishmulti[,2:3]
rownames(X)=NULL
colnames(X)=NULL
X <- as.matrix(X) # X : nxm matrix
n <- nrow(X)
m <- ncol(X)

##### settings #####
B          <- 50        # grid size for plot
Q          <- 500       # bootstrap resamples
alpha      <- 0.9
mu0_base   <- -0.25
rho_base   <- rep(0.5, m)
rho_tilde_base <- rep(1.0, m)
mu0.seq         <- seq(-0.5, 0, length.out = B)
rho_tilde.seq   <- seq(0.75, 1.25, length.out = B)   
h_vec <- rep(0.01*n^(-1/3), m) #kernel bandwidths

##### functions #####
## target function (quota)
est_QS <- function(X, mu0, rho, rho_tilde, alpha)
{
  n  <- nrow(X); m <- ncol(X)
  mu_hat    <- colMeans(X)
  Sigma_hat <- crossprod(scale(X, center = mu_hat, scale = FALSE)) / n
  
  num   <- sum(mu_hat * (rho_tilde - rho)) - mu0
  den   <- as.numeric(t(mu_hat * rho_tilde) %*% solve(Sigma_hat) %*% (mu_hat * rho_tilde))
  lambda_hat  <- 2 * num / den
  
  c_hat <- 1 - 0.5 * lambda_hat * as.numeric(solve(Sigma_hat) %*% (mu_hat * rho_tilde))
  
  H1_sq <- as.numeric(t(1 - c_hat) %*% Sigma_hat %*% (1 - c_hat))
  H1    <- sqrt(H1_sq)
  H2    <- sum(c_hat * mu_hat * rho_tilde) - sum(rho * mu_hat)
  VaR   <- sqrt(n) * (H1 * qnorm(alpha) + H2)
  
  c(c_hat, VaR_qs = VaR)
}

## target function (stop-loss)
est_SL <- function(X, mu0, rho, rho_tilde, alpha, h_vec, start=NULL)
{
  n <- nrow(X); m <- ncol(X)
  mu_hat <- colMeans(X)
  
  ##   estimating-equation function
  Lfun <- function(par) {
    d_vec  <- par[1:m]; lambda <- par[m + 1]
    scaled <- -(sweep(X, 2, d_vec)) / matrix(h_vec, n, m, byrow = TRUE)
    K_mat  <- pnorm(scaled);  barK <- 1 - K_mat
    X_wdg  <- X * K_mat + matrix(d_vec, n, m, TRUE) * barK
    retained      <- rowSums(X_wdg)
    mean_retained <- mean(retained)
    mean_barK     <- colMeans(barK)
    
    A_vec <- (2 / n) * as.numeric(crossprod(retained, barK))
    B_vec <- 2 * mean_retained * mean_barK
    C_vec <- lambda * rho_tilde * mean_barK
    L1    <- A_vec - B_vec - C_vec
    
    Delta <- (X - matrix(d_vec, n, m, TRUE)) * barK
    L2    <- mean(rowSums(sweep(Delta, 2, rho_tilde, "*"))) - sum(rho * mu_hat) - mu0
    
    c(L1, L2)
  }
  
  if(is.null(start)) {start <- c(apply(X, 2, function(x) mean(x[x>0])), 1)}
  sol <- nleqslv(start, Lfun, method = "Broyden")
  d_hat <- sol$x[1:m]
  X_wdg <- pmin(X, matrix(d_hat, n, m, TRUE))
  s_wdg <- rowSums(X_wdg)
  H1_sq <- mean(s_wdg^2) - mean(s_wdg)^2
  H1    <- sqrt(H1_sq)
  Delta_plus <- pmax(X - matrix(d_hat, n, m, TRUE), 0)
  H2 <- sum(rho_tilde * colMeans(Delta_plus)) - sum(rho * mu_hat)
  VaR <- sqrt(n) * (H1 * qnorm(alpha) + H2)
  c(d_hat, VaR_sl = VaR, lambda=sol$x[m+1])
}

## compute estimates (and bootstrap) for one grid point
compute_point <- function(mu0, rho_tilde, seed)
{
  set.seed(seed)
  ## point estimates
  qs <- est_QS(X, mu0,       rho_base, rho_tilde, alpha)
  sl <- est_SL(X, mu0, rho_base, rho_tilde, alpha, h_vec)
  
  ## bootstrap
  boot <- matrix(NA_real_, Q, (2*m+2)) 
  for (qq in 1:Q) {
    idx  <- sample.int(n, n, TRUE)
    Xb   <- X[idx, , drop = FALSE]
    qsb  <- est_QS(Xb, mu0, rho_base, rho_tilde, alpha)
    slb  <- est_SL(Xb, mu0, rho_base, rho_tilde, alpha, h_vec,
                   start=sl[c(1:m,m+2)])
    boot[qq, ] <- c(qsb[1:m], slb[1:m], qsb[m+1], slb[m+1])
  }
  sdvec <- apply(boot, 2, sd, na.rm = TRUE)
  
  list(mu0        = mu0,
       rho_tilde  = rho_tilde,
       c_hat      = qs[1:m],
       d_hat      = sl[1:m],
       VaR_qs     = qs[m+1],
       VaR_sl     = sl[m+1],
       sd_boot    = sdvec)
}

##### construct full task list #####
tasks <- list()

##  ceded share/ retention level vs mu0
for (i in seq_len(B))
  tasks[[length(tasks)+1]] <- list(mu0 = mu0.seq[i],
                                   rho_tilde = rho_tilde_base,
                                   tag = sprintf("mu0_%02d", i))

##  ceded share/ retention level vs rho_tilde
for (k in 1:m) {
  for (i in seq_len(B)) {
    rt     <- rho_tilde_base
    rt[k]  <- rho_tilde.seq[i]
    tasks[[length(tasks)+1]] <- list(mu0 = mu0_base,
                                     rho_tilde = rt,
                                     tag = sprintf("rhoT%02d_%02d", k, i))
  }
}

ntasks <- length(tasks)
cat("Total tasks:", ntasks, "\n")

##### execution with parallel computing #####
ncore <- max(1L, detectCores() - 1L)
cl    <- makeCluster(ncore)
clusterExport(cl,
              varlist = c("X","n","m","Q","alpha","rho_base","h_vec",
                          "est_QS","est_SL","compute_point","tasks"),
              envir = environment())
clusterEvalQ(cl, {library(Matrix); library(nleqslv); library(stats)})

results <- parLapply(cl, seq_along(tasks), function(j) {
  task <- tasks[[j]]
  compute_point(task$mu0, task$rho_tilde, seed = 1e6 + j)
})
stopCluster(cl)

##### output #####
extract_vec <- function(lst, name) vapply(lst, function(x) x[[name]], numeric(1))
c_mat  <- do.call(rbind, lapply(results, `[[`, "c_hat"))
d_mat  <- do.call(rbind, lapply(results, `[[`, "d_hat"))
VaR_qs <- extract_vec(results, "VaR_qs")
VaR_sl <- extract_vec(results, "VaR_sl")
sd_mat <- do.call(rbind, lapply(results, `[[`, "sd_boot"))

meta   <- do.call(rbind, lapply(results, function(x)
  cbind(mu0 = x$mu0, rho_tilde = t(x$rho_tilde))))

## Save
save(meta,c_mat,d_mat,VaR_qs,VaR_sl,sd_mat,file="real.Rda")

##### Plot #####
c_mat.l=pmax(c_mat-1.96*sd_mat[,(1:m)],0)
c_mat.u=pmin(c_mat+1.96*sd_mat[,(1:m)],1)
d_mat.l=d_mat-1.96*sd_mat[,((m+1):(2*m))]
d_mat.u=d_mat+1.96*sd_mat[,((m+1):(2*m))]
VaR_qs.l=VaR_qs-1.96*sd_mat[,(2*m+1)]
VaR_qs.u=VaR_qs+1.96*sd_mat[,(2*m+1)]
VaR_sl.l=VaR_sl-1.96*sd_mat[,(2*m+2)]
VaR_sl.u=VaR_sl+1.96*sd_mat[,(2*m+2)]

###### Plot ceded shares vs rho ######
eval.idx=1
title_eval=expression("Ceded shares vs" ~ tilde(rho)[1] ~ "(Danish fire data)")
xlab_eval=expression("Loading Factor " ~ tilde(rho)[1])
ylab_eval="Ceded shares"
sel.idx=((eval.idx*B+1):((eval.idx+1)*B))
param.eval=meta[sel.idx,(eval.idx+1)]
df1 <- data.frame(beta1   = param.eval,
                  theta   = c_mat[sel.idx,1],
                  lower   = c_mat.l[sel.idx,1],
                  upper   = c_mat.u[sel.idx,1],
                  estimator = "V1")
df2 <- data.frame(beta1   = param.eval,
                  theta   = c_mat[sel.idx,2],
                  lower   = c_mat.l[sel.idx,2],
                  upper   = c_mat.u[sel.idx,2],
                  estimator = "V2")
# Combine the three data frames
df_all <- bind_rows(df1, df2)
# Define estimator labels as expressions for the legend
estimator_labels <- c("V1" = expression(hat(c)[1]),
                      "V2" = expression(hat(c)[2]))
# Plot the estimated diversification measures against beta1 with confidence intervals
ggplot(df_all, aes(x = beta1, y = theta, group = estimator, linetype = estimator)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = estimator), alpha = 0.5) +
  geom_line(color = "black", size = 1) +
  scale_linetype_manual(
    values = c("V1" = "solid", "V2" = "dashed"),
    labels = estimator_labels
  ) +
  scale_fill_manual(
    values = c("V1" = "grey10", "V2" = "grey60"),
    labels = estimator_labels
  ) +
  labs(
    title = title_eval,
    x = xlab_eval,
    y = ylab_eval,
    linetype = "Estimator",
    fill = "Estimator"
  ) +
  ylim(0, 1)+
  theme_minimal() +
  theme(legend.position = "right") + 
  guides(linetype = guide_legend(keywidth = unit(1.5, "cm")))

###### Plot ceded shares vs mu0 ######
eval.idx=0
title_eval=expression("Ceded shares vs" ~ mu[0] ~ "(Danish fire data)")
xlab_eval=expression("Target Exposure " ~ mu[0])
ylab_eval="Ceded shares"
sel.idx=((eval.idx*B+1):((eval.idx+1)*B))
param.eval=meta[sel.idx,(eval.idx+1)]
df1 <- data.frame(beta1   = param.eval,
                  theta   = c_mat[sel.idx,1],
                  lower   = c_mat.l[sel.idx,1],
                  upper   = c_mat.u[sel.idx,1],
                  estimator = "V1")
df2 <- data.frame(beta1   = param.eval,
                  theta   = c_mat[sel.idx,2],
                  lower   = c_mat.l[sel.idx,2],
                  upper   = c_mat.u[sel.idx,2],
                  estimator = "V2")
df_all <- bind_rows(df1, df2)
estimator_labels <- c("V1" = expression(hat(c)[1]),
                      "V2" = expression(hat(c)[2]))
ggplot(df_all, aes(x = beta1, y = theta, group = estimator, linetype = estimator)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = estimator), alpha = 0.5) +
  geom_line(color = "black", size = 1) +
  scale_linetype_manual(
    values = c("V1" = "solid", "V2" = "dashed"),
    labels = estimator_labels
  ) +
  scale_fill_manual(
    values = c("V1" = "grey10", "V2" = "grey60"),
    labels = estimator_labels
  ) +
  labs(
    title = title_eval,
    x = xlab_eval,
    y = ylab_eval,
    linetype = "Estimator",
    fill = "Estimator"
  ) +
  ylim(0, 1)+
  theme_minimal() +
  theme(legend.position = "right") + 
  guides(linetype = guide_legend(keywidth = unit(1.5, "cm")))

###### Plot retention levels vs rho ######
eval.idx=1
title_eval=expression("Retention levels vs" ~ tilde(rho)[1] ~ "(Danish fire data)")
xlab_eval=expression("Loading Factor " ~ tilde(rho)[1])
ylab_eval="Retention levels"
sel.idx=((eval.idx*B+1):((eval.idx+1)*B))
param.eval=meta[sel.idx,(eval.idx+1)]
df1 <- data.frame(beta1   = param.eval,
                  theta   = d_mat[sel.idx,1],
                  lower   = d_mat.l[sel.idx,1],
                  upper   = d_mat.u[sel.idx,1],
                  estimator = "V1")
df2 <- data.frame(beta1   = param.eval,
                  theta   = d_mat[sel.idx,2],
                  lower   = d_mat.l[sel.idx,2],
                  upper   = d_mat.u[sel.idx,2],
                  estimator = "V2")
df_all <- bind_rows(df1, df2)
estimator_labels <- c("V1" = expression(hat(d)[1]),
                      "V2" = expression(hat(d)[2]))
ggplot(df_all, aes(x = beta1, y = theta, group = estimator, linetype = estimator)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = estimator), alpha = 0.5) +
  geom_line(color = "black", size = 1) +
  scale_linetype_manual(
    values = c("V1" = "solid", "V2" = "dashed"),
    labels = estimator_labels
  ) +
  scale_fill_manual(
    values = c("V1" = "grey10", "V2" = "grey60"),
    labels = estimator_labels
  ) +
  labs(
    title = title_eval,
    x = xlab_eval,
    y = ylab_eval,
    linetype = "Estimator",
    fill = "Estimator"
  ) +
  theme_minimal() +
  theme(legend.position = "right") + 
  guides(linetype = guide_legend(keywidth = unit(1.5, "cm")))

###### Plot retention levels vs mu0 ######
eval.idx=0
title_eval=expression("Retention levels vs" ~ mu[0] ~ "(Danish fire data)")
xlab_eval=expression("Target Exposure " ~ mu[0])
ylab_eval="Retention levels"
sel.idx=((eval.idx*B+1):((eval.idx+1)*B))
param.eval=meta[sel.idx,(eval.idx+1)]
df1 <- data.frame(beta1   = param.eval,
                  theta   = d_mat[sel.idx,1],
                  lower   = d_mat.l[sel.idx,1],
                  upper   = d_mat.u[sel.idx,1],
                  estimator = "V1")
df2 <- data.frame(beta1   = param.eval,
                  theta   = d_mat[sel.idx,2],
                  lower   = d_mat.l[sel.idx,2],
                  upper   = d_mat.u[sel.idx,2],
                  estimator = "V2")
df_all <- bind_rows(df1, df2)
estimator_labels <- c("V1" = expression(hat(d)[1]),
                      "V2" = expression(hat(d)[2]))
ggplot(df_all, aes(x = beta1, y = theta, group = estimator, linetype = estimator)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = estimator), alpha = 0.5) +
  geom_line(color = "black", size = 1) +
  scale_linetype_manual(
    values = c("V1" = "solid", "V2" = "dashed"),
    labels = estimator_labels
  ) +
  scale_fill_manual(
    values = c("V1" = "grey10", "V2" = "grey60"),
    labels = estimator_labels
  ) +
  labs(
    title = title_eval,
    x = xlab_eval,
    y = ylab_eval,
    linetype = "Estimator",
    fill = "Estimator"
  ) +
  theme_minimal() +
  theme(legend.position = "right") + 
  guides(linetype = guide_legend(keywidth = unit(1.5, "cm")))
###### Plot optimal VaR vs rho ######
eval.idx=1
title_eval=expression("Optimal VaR vs" ~ tilde(rho)[1] ~ "(Danish fire data)")
xlab_eval=expression("Loading Factor " ~ tilde(rho)[1])
ylab_eval="VaR of total cost"
sel.idx=((eval.idx*B+1):((eval.idx+1)*B))
param.eval=meta[sel.idx,(eval.idx+1)]
df1 <- data.frame(beta1   = param.eval,
                  theta   = VaR_qs[sel.idx],
                  lower   = VaR_qs.l[sel.idx],
                  upper   = VaR_qs.u[sel.idx],
                  estimator = "V1")
df2 <- data.frame(beta1   = param.eval,
                  theta   = VaR_sl[sel.idx],
                  lower   = VaR_sl.l[sel.idx],
                  upper   = VaR_sl.u[sel.idx],
                  estimator = "V2")
df_all <- bind_rows(df1, df2)
estimator_labels <- c("V1" = "Quota-share",
                      "V2" = "Stop-loss")
ggplot(df_all, aes(x = beta1, y = theta, group = estimator, linetype = estimator)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = estimator), alpha = 0.5) +
  geom_line(color = "black", size = 1) +
  scale_linetype_manual(
    values = c("V1" = "solid", "V2" = "dashed"),
    labels = estimator_labels
  ) +
  scale_fill_manual(
    values = c("V1" = "grey10", "V2" = "grey60"),
    labels = estimator_labels
  ) +
  labs(
    title = title_eval,
    x = xlab_eval,
    y = ylab_eval,
    linetype = "Type",
    fill = "Type"
  ) +
  theme_minimal() +
  theme(legend.position = "right") + 
  guides(linetype = guide_legend(keywidth = unit(1.5, "cm")))

###### Plot optimal VaR vs mu0 ######
eval.idx=0
title_eval=expression("Optimal VaR vs" ~ mu[0] ~ "(Danish fire data)")
xlab_eval=expression("Target Exposure " ~ mu[0])
ylab_eval="VaR of total cost"
sel.idx=((eval.idx*B+1):((eval.idx+1)*B))
param.eval=meta[sel.idx,(eval.idx+1)]
df1 <- data.frame(beta1   = param.eval,
                  theta   = VaR_qs[sel.idx],
                  lower   = VaR_qs.l[sel.idx],
                  upper   = VaR_qs.u[sel.idx],
                  estimator = "V1")
df2 <- data.frame(beta1   = param.eval,
                  theta   = VaR_sl[sel.idx],
                  lower   = VaR_sl.l[sel.idx],
                  upper   = VaR_sl.u[sel.idx],
                  estimator = "V2")
df_all <- bind_rows(df1, df2)
estimator_labels <- c("V1" = "Quota-share",
                      "V2" = "Stop-loss")
ggplot(df_all, aes(x = beta1, y = theta, group = estimator, linetype = estimator)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = estimator), alpha = 0.5) +
  geom_line(color = "black", size = 1) +
  scale_linetype_manual(
    values = c("V1" = "solid", "V2" = "dashed"),
    labels = estimator_labels
  ) +
  scale_fill_manual(
    values = c("V1" = "grey10", "V2" = "grey60"),
    labels = estimator_labels
  ) +
  labs(
    title = title_eval,
    x = xlab_eval,
    y = ylab_eval,
    linetype = "Type",
    fill = "Type"
  ) +
  theme_minimal() +
  theme(legend.position = "right") + 
  guides(linetype = guide_legend(keywidth = unit(1.5, "cm")))
