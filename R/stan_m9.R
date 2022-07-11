rm(list = ls())

library(dplyr)
library(ggplot2)
library(matrixStats)
library(cmdstanr)
library(posterior)
library(bayesplot)

set_cmdstan_path('/Users/nico/GROUNDMOTION/SOFTWARE/cmdstan-2.29.2')
cmdstan_path()
cmdstan_version()

`%notin%` <- Negate(`%in%`)

func_lh <- function(x, coeff, mb = 5, delta = 0.1) {
  return(coeff[1] + coeff[2] * (x - mb) + (coeff[3] - coeff[2]) * delta *log(1 + exp((x - mb) / delta)))
}

func_lh2 <- function(x, coeff, mb = 5, delta = 0.1) {
  return(coeff[1] + coeff[2] * (x - mb) + 
           (coeff[3] - coeff[2]) * delta * sapply(x, function (r) {logSumExp(c(log(1), (r - mb) / delta))}))
  
}

q25 <- function(y) quantile(y, 0.25)
q75 <- function(y) quantile(y, 0.75)

path_base <- "/Users/nico/GROUNDMOTION/PROJECTS/SUBDUCTION/DURATION/"
path_stan <- file.path(path_base, "STAN_Jun22")

path_m9 <- file.path(path_base, "M9")
files_m9 <- list.files(path_m9, pattern = "*.dat")

n_m9 <- length(files_m9)

tmp <- read.csv(file.path(path_m9, files_m9[1]))
dim(tmp)
n_r <- nrow(tmp)

### M9 simulaions are one a grid of points 50 (lat) by 28 (lon)
grid_m9 <- tidyr::expand_grid(x = 1:50, y = 1:28)

dat_m9 <- data.frame()
for(i in 1:n_m9) {
  tmp <- readr::read_csv(file.path(path_m9, files_m9[i]))
  tmp$idx_lat <- grid_m9$x
  tmp$idx_lon <- grid_m9$y
  dat_m9 <- rbind(dat_m9, tmp)
}
#names(dat_m9) <- names(tmp)
head(dat_m9)

ggplot(dat_m9[dat_m9$idx_lat == 20 & dat_m9$idx_lon > 11,]) +
  geom_point(mapping = aes(x = RCD, y = Ds575_x, color = `# RealizationID`)) +
  scale_x_log10() + scale_y_log10()


k <- 30 # latitude index

#data_used <- dat_m9[dat_m9$idx_lat == k & dat_m9$idx_lon > 11,]
data_used <- dat_m9[dat_m9$idx_lat >= 10 & dat_m9$idx_lat <= 40 & dat_m9$idx_lon > 11,]
data_used$lnD75 <- rowMeans(cbind(log(data_used$Ds575_x), log(data_used$Ds575_y)))

ggplot() +
  geom_point(data = dat_m9, mapping = aes(x = Longitude, y = Latitude), color = 'red') +
  geom_point(data = data_used, mapping = aes(x = Longitude, y = Latitude))

#write.csv(data_used, file = file.path(path_stan, 'RESULTS', 'data_m9_Lat10-40.csv'))


n_rec <- nrow(data_used)
eq <- as.numeric(factor(data_used$`# RealizationID`))
n_eq <- max(eq)
mageq <- rep(9, n_eq)

bin_dist  <- c(0,50,100,150,200,300)
idx_bin <- as.numeric(mltools::bin_data(data_used$RCD, bins = bin_dist))
n_bin <- max(idx_bin)

data_list <- list(N = n_rec,
                  NEQ = n_eq,
                  NBIN = n_bin,
                  MEQ = mageq,
                  R = data_used$RCD,
                  Y = data_used$lnD75,
                  eq = eq,
                  idx_bin = idx_bin
)

model <- 'gmm_dsource_hier_Rlheq_logn_sigmaR3'
print(model)
file <- file.path(path_stan, sprintf('%s.stan', model))
mod <- cmdstan_model(file, include_paths = c(file.path(path_stan, 'stan_include')))

fit <- mod$sample(
  data = data_list,
  seed = 8472,
  chains = 4,
  iter_sampling = 200,
  iter_warmup = 500,
  refresh = 10,
  max_treedepth = 10,
  adapt_delta = 0.95,
  parallel_chains = 2,
  init = init_list
)
fit$cmdstan_diagnose()
#fit$save_object(file = file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat%d_%s.RDS', k, model)))
fit$save_object(file = file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat10-40_%s.RDS', model)))


draws <- fit$draws()
summarise_draws(subset(draws, variable=c('^s', 'mu', '^b', '^r1', 'r2', 'R_b','^rp', 'R_s', 'delta_s'), regex=TRUE))
summarise_draws(subset(draws, variable=c('eqterm[22]'), regex=FALSE))

mcmc_trace(draws, pars = c("lp__"))

mag_pred <- seq(4,9.5,0.25)
mu_ln_dsigma <- mean(subset(as_draws_matrix(draws), variable=c('mu_ln_dsigma'), regex=FALSE))
ln_denom <- log(3.2) + log(4.9) + log(10^6)
ln_moment <- log(10) * (1.5 * mag_pred + 16.05);
c1_pred <- exp(-1/3 * (mu_ln_dsigma - ln_moment) - ln_denom)
df_ds <- data.frame(M = mag_pred, ds = c1_pred, ln_dsigma = mu_ln_dsigma)

func_ds <- function(x) exp(-1/3 * (x - log(10) * (1.5 * 9 + 16.05)) - ln_denom)

df_eq <- data.frame(M = mageq,
                    c1_mean = colMeans(subset(as_draws_matrix(draws), variable=c('^c1'), regex=TRUE)),
                    c1_med = colMedians(subset(as_draws_matrix(draws), variable=c('^c1'), regex=TRUE)),
                    lnds_mean = colMeans(subset(as_draws_matrix(draws), variable=c('^ln_dsigma'), regex=TRUE)),
                    lnds_med = colMedians(subset(as_draws_matrix(draws), variable=c('^ln_dsigma'), regex=TRUE)),
                    dln_ds_mean = colMeans(subset(as_draws_matrix(draws), variable=c('^dln_ds'), regex=TRUE)),
                    dln_ds_med = colMedians(subset(as_draws_matrix(draws), variable=c('^dln_ds'), regex=TRUE))
)

ggplot(df_eq) +
  geom_density(mapping = aes(x = c1_mean)) +
  geom_density(mapping = aes(x = c1_med), color = 'red')

mcmc_areas_ridges(draws, 
                  regex_pars = c("^ln_dsigma*", "mu_ln_dsigma"))

mcmc_areas_ridges(draws, 
                  regex_pars = c("^ln_dsigma*", "mu_ln_dsigma"),
                  transformations = "func_ds")

mcmc_areas_ridges(draws, 
                  regex_pars = c("^shape_path"))

mcmc_areas_ridges(draws, 
                  regex_pars = c("^rp"))

mcmc_areas_ridges(draws, 
                  regex_pars = c("^R_b"))

mcmc_hist(draws, pars = c("r1", "r2"))

mcmc_areas_ridges(draws, 
                  regex_pars = c("^eqterm*"))

ggplot(df_ds) +
  geom_line(mapping = aes(x = M, y = ds)) +
  scale_y_log10()

##############
dist_pred <- seq(0,300,by = 10)
# r1_est <- mean(subset(as_draws_matrix(draws), variable=c('r1'), regex=FALSE))
# dpath_pred <- r1_est * dist_pred + exp(-1/3 * (mu_ln_dsigma - log(10) * (1.5 * 9 + 16.05)) - ln_denom)
# dpath_pred <- r1_est * dist_pred + exp(-1/3 * (mu_ln_dsigma - log(10) * (1.5 * 9 + 16.05)) - ln_denom)
coeff_r <- colMeans(subset(as_draws_matrix(draws), variable=c('r1', 'r2'), regex=FALSE))
c1_pred <- exp(-1/3 * (mu_ln_dsigma - log(10) * (1.5 * 9 + 16.05)) - ln_denom)
r_b <- mean(subset(as_draws_matrix(draws), variable=c('R_b')))
dpath_pred <- func_lh(dist_pred, c(c1_pred + coeff_r[1] * r_b, coeff_r), mb = r_b, delta = 1)

df_dp <- data.frame(R = dist_pred, lnD = log(dpath_pred))

df_rec <- data.frame(R = data_list$R, eq = eq,
                     lnD75 = data_used$lnD75,
                     resid_mean = colMeans(subset(as_draws_matrix(draws), variable=c('^resid'), regex=TRUE)),
                     resid_med = colMedians(subset(as_draws_matrix(draws), variable=c('^resid'), regex=TRUE))
)

fl <- loess(resid_mean ~ log(R), df_rec)
df_rec$smooth <- predict(fl)
ggplot() +
  geom_point(data = df_rec, mapping = aes(x = R, y = resid_mean, color = as.factor(eq))) +
  geom_line(data = df_rec, mapping = aes(x = R, y = smooth)) +
  scale_x_log10()

ggplot() +
  geom_point(data = df_rec, mapping = aes(x = R, y = exp(lnD75), color = as.factor(eq))) +
  geom_line(data = df_dp, mapping = aes(x = R, y = exp(lnD))) +
  scale_x_log10() + scale_y_log10()


###
idx <- 1:n_rec
y_rep <- log(subset(as_draws_matrix(draws), variable=c('^y_'), regex=TRUE)[,idx])
y <- data_used$lnD75[idx]

ppc_dens_overlay(y, y_rep[1:25, ])
ppc_ecdf_overlay(y, y_rep[1:25, ])
ppc_stat(y, y_rep)
ppc_stat(y, y_rep, stat = "q25")
ppc_stat(y, y_rep, stat = "q75")
ppc_stat_2d(y, y_rep, stat = c("q25", "q75"))


##### loo
k <- 20
model <- 'gmm_dsigma_hier_Rlh_logn'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat%d_%s.RDS', k, model)))
draws_logn <- fit$draws()
loo_logn <- fit$loo()

model <- 'gmm_dsigma_hier_Rlh'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat%d_%s.RDS', k, model)))
draws_gamma <- fit$draws()
loo_gamma <- fit$loo()

model <- 'gmm_dsigma_hier'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat%d_%s.RDS', k, model)))
draws_gamma1 <- fit$draws()
loo_gamma1 <- fit$loo()

model <- 'gmm_dsource_hier'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat%d_%s.RDS', k, model)))
draws_ds <- fit$draws()
loo_ds <- fit$loo()

model <- 'gmm_dsource_hier_logn'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat%d_%s.RDS', k, model)))
draws_ds_logn <- fit$draws()
loo_ds_logn <- fit$loo()

model <- 'gmm_dsource_hier_mix'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat%d_%s.RDS', k, model)))
draws_ds_mix <- fit$draws()
loo_ds_mix <- fit$loo()

loo::loo_compare(loo_gamma1, loo_gamma, loo_logn, loo_ds, loo_ds_logn, loo_ds_mix)

plot(colMeans(subset(as_draws_matrix(draws_gamma), variable=c('^c1'), regex=TRUE)),
     colMeans(subset(as_draws_matrix(draws_ds), variable=c('^c1'), regex=TRUE)))

k <- 30
model <- 'gmm_dsigma_hier_Rlh'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat%d_%s.RDS', k, model)))
draws_gamma_b <- fit$draws()


plot(colMeans(subset(as_draws_matrix(draws_gamma), variable=c('^ln_dsigma'), regex=TRUE)),
     colMeans(subset(as_draws_matrix(draws_gamma_b), variable=c('^ln_dsigma'), regex=TRUE)))


df_plot <- data.frame(post_k25 = subset(as_draws_matrix(draws_gamma), variable=c('ln_dsigma[1]'), regex=FALSE),
                      post_k30 = subset(as_draws_matrix(draws_gamma_b), variable=c('ln_dsigma[1]'), regex=FALSE)
)

ggplot(df_plot) + 
  geom_density(mapping = aes(x = ln_dsigma.1.)) +
  geom_density(mapping = aes(x = ln_dsigma.1..1), color = 'red')

