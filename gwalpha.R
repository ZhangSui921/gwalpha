# Title: Geographically Weighted Cronbach’s Alpha (GWalpha): an Exploratory Local Measure of Reliability for Scale Construction
## Sui Zhang and Ziqi Li
## Corresponce: Sui Zhang (sui.zhang@manchester.ac.uk)

## Load packages
library(sp)
library(MASS)
library(geoR)
library(psych)

## GWalpha function
gwalpha <- function(x, data, kernel = 'bisquare', adaptive = T, bw, ci = F, p = .95, nsims = 1000) {
  coords <- coordinates(data)
  data <- as.data.frame(data)
  n <- nrow(coords)
  x_v <- as.matrix(data[, x])
  n_v <- ncol(x_v)
  mv <- numeric(n)
  mcov <- numeric(n)
  gwalpha <- numeric(n)
  if (ci) {
    b_mv <- matrix(nrow = n, ncol = nsims)
    b_mcov <- matrix(nrow = n, ncol = nsims)
    b_gwalpha <- matrix(nrow = n, ncol = nsims)
  }

  d <- as.matrix(dist(coords))
  for (i in 1:n) {
    d_i <- d[, i]
    if (adaptive) {
      h <- sort(d_i)[bw]
    } else {
      h <- bw
    }
    if (kernel == 'bisquare') {
      w_i <- ifelse(d_i > h, 0, (1 - (d_i/h)^2)^2)
    } else if (kernel == 'gaussian') {
      w_i <- exp(-0.5 * (d_i/h)^2)
    } else if (kernel == 'exponential') {
      w_i <- exp(-d_i/h)
    } else {
      w_i <- ifelse(d_i > h, 0, 1)
    }

    wi <- w_i / sum(w_i)
    wcov <- cov.wt(x_v, wt = wi)$cov
    mv[i] <- sum(diag(wcov)) / n_v
    mcov[i] <- (sum(wcov) - sum(diag(wcov))) / (n_v * (n_v - 1))
    gwalpha[i] <- n_v * mcov[i] / (mv[i] + (n_v - 1) * mcov[i])

    if (ci) {
      b_wi <- wi[order(d_i)[1:bw]]
      for (j in 1:nsims) {
        b <- sample(order(d_i)[1:bw], replace = T)
        b_v <- x_v[b, , drop = F]
        b_wcov <- cov.wt(b_v, wt = b_wi)$cov
        b_mv[i, j] <- sum(diag(b_wcov)) / n_v
        b_mcov[i, j] <- (sum(b_wcov) - sum(diag(b_wcov))) / (n_v * (n_v - 1))
        b_gwalpha[i, j] <- n_v * b_mcov[i, j] / (b_mv[i, j] + (n_v - 1) * b_mcov[i, j])
      }
    }
  }

  if (ci) {
    gwalpha_u <- apply(b_gwalpha, 1, function(x) quantile(x, p))
    res <- data.frame(gwalpha, gwalpha_u, coords)
  } else {
    res <- data.frame(gwalpha, coords)
  }
  
  return(res)
}

## Bandwidth search for GWalpha

bw_gwalpha <- function(x, data, kernel = 'bisquare', adaptive = T, int, crit, p = .95, nsims = 1000) {
  bws <- data.frame(
    bw = seq(int, nrow(data), int),
    n_crit = sapply(
      seq(int, nrow(data), int),
      function(bw) {
        sum(gwalpha(x = x, data = data, kernel = kernel, adaptive = adaptive, bw = bw, ci = T, p = p, nsims = nsims)$gwalpha_u < crit)
      }
    )
  )
  
  bw_opt <- bws$bw[which.max(bws$n_crit)]
  
  bw_opt
}

## Simulation (N realisations)
N <- 1000
sim <- as.data.frame(matrix(nrow = N))
n_tot <- nrow(sim)

set.seed(000)
sim$r <- runif(N)
sim$theta <- runif(N, max = 2*pi)
sim$u <- 12.5 + 12.5 * sqrt(sim$r) * cos(sim$theta)
sim$v <- 12.5 + 12.5 * sqrt(sim$r) * sin(sim$theta)
sim <- sim[, -1]

sample <- 200
rmse <- expand.grid(u = seq(5, sample, 5), v = seq(25, n_tot, 25))

for (r in c(5, 10, 15)) {
  cov_mats <- list()
  obs <- list()
  obs_df <- data.frame(matrix(nrow = 0, ncol = 5))
  
  set.seed(111)
  sim$l1 <- grf(n_tot, sim[, c('u','v')], cov.model = 'gaussian', cov.pars = c(1, r))$data
  sim$l1 <- sim$l1 - min(sim$l1) + 2
  sim$l1 <- sim$l1 / max(sim$l1)
  
  set.seed(222)
  sim$l2 <- grf(n_tot, sim[, c('u','v')], cov.model = 'gaussian', cov.pars = c(1, r))$data
  sim$l2 <- sim$l2 - min(sim$l2) + 2
  sim$l2 <- sim$l2 / max(sim$l2)
  
  set.seed(333)
  sim$l3 <- grf(n_tot, sim[, c('u','v')], cov.model = 'gaussian', cov.pars = c(1, r))$data
  sim$l3 <- sim$l3 - min(sim$l3) + 2
  sim$l3 <- sim$l3 / max(sim$l3)
  
  for (i in 1:n_tot) {
    l <- as.matrix(sim[i, c('l1', 'l2', 'l3')])
    cov_e <- diag(0.4, 3)
    cov_mats[[i]] <- t(l) %*% l + cov_e
    sim$ta[i] <- (3 * mean(c(cov_mats[[i]][1, 2:3],
                             cov_mats[[i]][2, 3]))) / (mean(diag(cov_mats[[i]])) + 2 * 
                                                         mean(c(cov_mats[[i]][1, 2:3], cov_mats[[i]][2, 3])))
    obs[[i]] <- mvrnorm(n = sample, mu = rep(2.5, 3), Sigma = cov_mats[[i]])
    u_i <- rep(sim[i, 'u'], each = sample)
    v_i <- rep(sim[i, 'v'], each = sample)
    obs_i <- cbind(obs[[i]], u_i, v_i)
    obs_df <- rbind(obs_df, obs_i)
  }
  
  colnames(obs_df) <- c('x1', 'x2', 'x3', 'u', 'v')
  rmse_sa <- numeric(sample/5)
  rmse_gwsa <- data.frame(matrix(nrow = sample/5, ncol = n_tot/25))
  
  for (n in seq(5, sample, 5)) {
    obs_df_n <- obs_df[sapply(split(1:(sample * n_tot), rep(1:n_tot, each = sample)), 
                              function(x) sample(x, n)), ]
    obs_sp_n <- SpatialPointsDataFrame(obs_df_n[, 4:5], obs_df_n)
    
    if (n > 1) {
      sim$sa <- sapply(split(obs_df_n, rep(1:n_tot, each = n)), 
                       function(x) alpha(x[, 1:3], check.keys = F)[[1]][1])
      sim$sa <- as.numeric(sim$sa)
      rmse_sa[n/5] <- sqrt(mean((sim$ta - sim$sa)^2))
    } else {
      rmse_sa[n/5] <- Inf
    }
    
    for (b in seq(25, n_tot, 15)) {
      bw <- b * n
      sim$gwsa <- gwalpha(x = c('x1', 'x2', 'x3'), data = obs_sp_n, bw = bw)[, 1]
      rmse_gwsa[n/5, (b - 10)/15] <- sqrt(mean((sim$ta - sim$gwsa)^2))
    }
  }
  
  rmse_sa <- rep(rmse_sa, n_tot/25)
  rmse_gwsa <- unlist(rmse_gwsa)
  diff <- rmse_gwsa - rmse_sa
  rmse <- cbind(rmse, diff)
  colnames(rmse)[colnames(rmse) == 'diff'] <- paste0('diff_', r)
}

## Simulation (One realisation)
p_sim <- sim[, c('u', 'v')]
for (r in c(5, 10, 15)) {
  cov_mats <- list()
  obs <- list()
  obs_df <- data.frame(matrix(nrow = 0, ncol = 5))
  
  set.seed(111)
  p_sim$l1 <- grf(n_tot, p_sim[, c('u', 'v')], cov.model = 'gaussian', cov.pars = c(1, r))$data
  p_sim$l1 <- p_sim$l1 - min(p_sim$l1) + 2
  p_sim$l1 <- p_sim$l1 / max(p_sim$l1)
  
  set.seed(222)
  p_sim$l2 <- grf(n_tot, p_sim[, c('u', 'v')], cov.model = 'gaussian', cov.pars = c(1, r))$data
  p_sim$l2 <- p_sim$l2 - min(p_sim$l2) + 2
  p_sim$l2 <- p_sim$l2 / max(p_sim$l2)
  
  set.seed(333)
  p_sim$l3 <- grf(n_tot, p_sim[, c('u', 'v')], cov.model = 'gaussian', cov.pars = c(1, r))$data
  p_sim$l3 <- p_sim$l3 - min(p_sim$l3) + 2
  p_sim$l3 <- p_sim$l3 / max(p_sim$l3)
  
  for (i in 1:n_tot) {
    l <- as.matrix(p_sim[i, 3:5])
    cov_e <- diag(0.4, 3)
    cov_mats[[i]] <- t(l) %*% l + cov_e
    p_sim$ta[i] <- (3 * mean(c(cov_mats[[i]][1, 2:3], cov_mats[[i]][2, 3]))) / 
      (mean(diag(cov_mats[[i]])) + 2 * mean(c(cov_mats[[i]][1, 2:3], cov_mats[[i]][2, 3])))
    obs[[i]] <- mvrnorm(n = sample, mu = rep(2.5, 3), Sigma = cov_mats[[i]])
    u_i <- rep(p_sim[i, 'u'], each = sample)
    v_i <- rep(p_sim[i, 'v'], each = sample)
    obs_i <- cbind(obs[[i]], u_i, v_i)
    obs_df <- rbind(obs_df, obs_i)
  }
  
  colnames(obs_df) <- c('x1', 'x2', 'x3', 'u', 'v')
  obs_df_n <- obs_df[sapply(split(1:(sample * n_tot), rep(1:n_tot, each = sample)),
                            function(x) sample(x, 1)), ]
  obs_sp_n <- SpatialPointsDataFrame(obs_df_n[, 4:5], obs_df_n)
  
  p_sim$a <- alpha(obs_df_n[, 1:3])[[1]][[1]]
  p_sim$gwalpha <- gwalpha(x = c('x1', 'x2', 'x3'), data = obs_sp_n, bw = 30 * r)[, 1]
  
  colnames(p_sim)[colnames(p_sim) == 'ta'] <- paste0('ta', r)
  colnames(p_sim)[colnames(p_sim) == 'a'] <- paste0('a', r)
  colnames(p_sim)[colnames(p_sim) == 'gwalpha'] <- paste0('gwa', r)
}

## Real-world example
dt11 <- read.csv('data_bes2011.csv')
x <- c("willHelp", "closeKnit", "trust", "solComProb", "relGroups")
dt11_sp <- SpatialPointsDataFrame(dt11[,2:3], dt11)
bw_opt <- bw_gwalpha(x = x, data = dt11_sp, int = round(.01 * nrow(dt11_sp)), crit = .8)
gwa11_opt <- gwalpha(x = x, data = dt11_sp, bw = bw_opt, ci = T)
