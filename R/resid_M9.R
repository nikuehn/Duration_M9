########## Select subsets
tect <- "if"
if(tect == "if") {
  data_used_sub <- data_kbcg2[data_kbcg2$Vs30_Selected_for_Analysis_m_s > 0
                          & data_kbcg2$Intra_Inter_Flag == 0
                          ,]
} else {
  data_used_sub <- data_kbcg2[data_kbcg2$Vs30_Selected_for_Analysis_m_s > 0
                          & data_kbcg2$Intra_Inter_Flag != 0
                          ,]
}

eq_sub <- as.numeric(factor(data_used_sub$NGAsubEQID))
mageq <- unique(cbind(data_used_sub$NGAsubEQID, data_used_sub$Earthquake_Magnitude))[,2]

data_list_sub <- list(N = length(eq_sub),
                  NEQ = max(eq_sub),
                  MEQ = mageq,
                  R = data_used_sub$ClstD_km,
                  Y = log(data_used_sub$D_75_old),
                  eq = eq_sub
)

model <- 'gmm_dsource_hier_Rlh_shapeR'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat10-40_%s.RDS', model)))
draws_tmp <- fit$draws()
fit$diagnostic_summary()

summarise_draws(subset(draws_tmp, variable=c('sigma_ln_dsigma'), regex=FALSE))

model <- 'gmm_dsource_hier_Rlh_shapeR_gqresid'
print(model)
file <- file.path(path_stan, sprintf('%s.stan', model))
mod <- cmdstan_model(file, include_paths = c(file.path(path_stan, 'stan_include')))

fit <- mod$generate_quantities(
  fitted_params = draws_tmp,
  data = data_list_sub,
  seed = 8472
)
draws_resid <- fit$draws()

df_rec <- data.frame(R = data_used_sub$ClstD_km, eq = eq_sub,
                     M = mageq[eq_sub],
                     VS = data_used_sub$Vs30_Selected_for_Analysis_m_s,
                     lnD75 = data_used_sub$D_75,
                     #resid_1 = colMeans(subset(as_draws_matrix(draws_resid),
                     #                             variable=c('^resid_1'), regex=TRUE)),
                     #resid_2 = colMeans(subset(as_draws_matrix(draws_resid),
                     #                          variable=c('^resid_2'), regex=TRUE)),
                     resid = colMeans(subset(as_draws_matrix(draws_resid),
                                               variable=c('^resid'), regex=TRUE))
)
resid_mean <- mean(df_rec$resid)

# ggplot(df_rec2) +
#   geom_point(aes(x = R, y = resid_1)) +
#   geom_point(aes(x = R, y = resid_2), color = 'red') +
#   geom_point(aes(x = R, y = resid_3), color = 'blue')
# 
# 
# ggplot(df_rec2) +
#   geom_point(aes(x = M, y = resid_1)) +
#   geom_point(aes(x = M, y = resid_2), color = 'red') +
#   geom_point(aes(x = M, y = resid_3), color = 'blue')
# 
# plot(df_rec2$resid_2, df_rec2$resid_3)
# df_rec$c1 - df_rec$c1c
# 
# mean(df_rec2$resid_1)
# mean(df_rec2$resid_2)
# mean(df_rec2$resid_3)
# 
# mean(df_rec$resid_1)
# mean(df_rec$resid_2)
# mean(df_rec$resid_3)

model <- "Gamma"
xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- "residual"
fl <- loess(resid ~ R, df_rec)
df_rec$smooth <- predict(fl)
p <- ggplot(df_rec) +
  geom_point(aes(x = R, y = resid), color = 'gray') +
  geom_line(aes(x = R, y = smooth)) +
  geom_hline(yintercept = resid_mean, linetype = "dashed") +
  labs(x = xlab, y = ylab, title = model) +
  theme(
    legend.position = "none"
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_residNGA_Rrup_%s.pdf', model)), p,
       width = wid, height = asp * wid)

xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- "residual"
fl <- loess(resid ~ R, df_rec)
df_rec$smooth <- predict(fl)
p <- ggplot(df_rec) +
  geom_point(aes(x = R, y = resid), color = 'gray') +
  geom_line(aes(x = R, y = smooth)) +
  geom_hline(yintercept = resid_mean, linetype = "dashed") +
  labs(x = xlab, y = ylab, title = model) +
  theme(
    legend.position = "none"
  ) +
  scale_x_log10()
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_residNGA_lnRrup_%s.pdf', model)), p,
       width = wid, height = asp * wid)

xlab <- "M"
ylab <- "residual"
fl <- loess(resid ~ M, df_rec)
df_rec$smooth <- predict(fl)
p <- ggplot(df_rec) +
  geom_point(aes(x = M, y = resid), color = 'gray') +
  geom_line(aes(x = M, y = smooth)) +
  geom_hline(yintercept = resid_mean, linetype = "dashed") +
  labs(x = xlab, y = ylab, title = model) +
  theme(
    legend.position = "none"
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_residNGA_M_%s.pdf', model)), p,
       width = wid, height = asp * wid)

xlab <- expression(atop(paste(V[S30] (m/s))))
ylab <- "residual"
fl <- loess(resid ~ VS, df_rec)
df_rec$smooth <- predict(fl)
p <- ggplot(df_rec) +
  geom_point(aes(x = VS, y = resid), color = 'gray') +
  geom_line(aes(x = VS, y = smooth)) +
  geom_hline(yintercept = resid_mean, linetype = "dashed") +
  labs(x = xlab, y = ylab, title = model) +
  theme(
    legend.position = "none"
  ) +
  scale_x_log10()
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_residNGA_Vs30_%s.pdf', model)), p,
       width = wid, height = asp * wid)


####################################
model <- 'gmm_dsource_hier_Rlh_logn_sigmaR3'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat10-40_%s.RDS', model)))
draws_tmp <- fit$draws()
fit$diagnostic_summary()

summarise_draws(subset(draws_tmp, variable=c('sigma_ln_dsigma'), regex=FALSE))

model <- 'gmm_dsource_hier_Rlh_logn_sigmaR3_gqresid'
print(model)
file <- file.path(path_stan, sprintf('%s.stan', model))
mod <- cmdstan_model(file, include_paths = c(file.path(path_stan, 'stan_include')))

fit <- mod$generate_quantities(
  fitted_params = draws_tmp,
  data = data_list_sub,
  seed = 8472
)
draws_resid <- fit$draws()

df_rec <- data.frame(R = data_used_sub$ClstD_km, eq = eq_sub,
                     M = mageq[eq_sub],
                     VS = data_used_sub$Vs30_Selected_for_Analysis_m_s,
                     lnD75 = data_used_sub$D_75,
                     #resid_1 = colMeans(subset(as_draws_matrix(draws_resid),
                     #                             variable=c('^resid_1'), regex=TRUE)),
                     #resid_2 = colMeans(subset(as_draws_matrix(draws_resid),
                     #                          variable=c('^resid_2'), regex=TRUE)),
                     resid = colMeans(subset(as_draws_matrix(draws_resid),
                                             variable=c('^resid_2\\['), regex=TRUE))
)
resid_mean <- mean(df_rec$resid)

model <- "Lognormal"
xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- "residual"
fl <- loess(resid ~ R, df_rec)
df_rec$smooth <- predict(fl)
p <- ggplot(df_rec) +
  geom_point(aes(x = R, y = resid), color = 'gray') +
  geom_line(aes(x = R, y = smooth)) +
  geom_hline(yintercept = resid_mean, linetype = "dashed") +
  labs(x = xlab, y = ylab, title = model) +
  theme(
    legend.position = "none"
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_residNGA_Rrup_%s.pdf', model)), p,
       width = wid, height = asp * wid)

xlab <- expression(atop(paste(R[RUP] (km))))
fl <- loess(resid ~ R, df_rec)
df_rec$smooth <- predict(fl)
p <- ggplot(df_rec) +
  geom_point(aes(x = R, y = resid), color = 'gray') +
  geom_line(aes(x = R, y = smooth)) +
  geom_hline(yintercept = resid_mean, linetype = "dashed") +
  labs(x = xlab, y = ylab, title = model) +
  theme(
    legend.position = "none"
  ) +
  scale_x_log10()
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_residNGA_lnRrup_%s.pdf', model)), p,
       width = wid, height = asp * wid)

xlab <- "M"
fl <- loess(resid ~ M, df_rec)
df_rec$smooth <- predict(fl)
p <- ggplot(df_rec) +
  geom_point(aes(x = M, y = resid), color = 'gray') +
  geom_line(aes(x = M, y = smooth)) +
  geom_hline(yintercept = resid_mean, linetype = "dashed") +
  labs(x = xlab, y = ylab, title = model) +
  theme(
    legend.position = "none"
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_residNGA_M_%s.pdf', model)), p,
       width = wid, height = asp * wid)

xlab <- expression(atop(paste(V[S30] (m/s))))
fl <- loess(resid ~ VS, df_rec)
df_rec$smooth <- predict(fl)
p <- ggplot(df_rec) +
  geom_point(aes(x = VS, y = resid), color = 'gray') +
  geom_line(aes(x = VS, y = smooth)) +
  geom_hline(yintercept = resid_mean, linetype = "dashed") +
  labs(x = xlab, y = ylab, title = model) +
  theme(
    legend.position = "none"
  ) +
  scale_x_log10()
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_residNGA_Vs30_%s.pdf', model)), p,
       width = wid, height = asp * wid)

####################################
model <- 'gmm_dsource_hier_Rlh_mix_sigmaR3'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat10-40_%s.RDS', model)))
draws_tmp <- fit$draws()
fit$diagnostic_summary()

summarise_draws(subset(draws_tmp, variable=c('sigma_ln_dsigma'), regex=FALSE))

model <- 'gmm_dsource_hier_Rlh_mix_sigmaR3_gqresid'
print(model)
file <- file.path(path_stan, sprintf('%s.stan', model))
mod <- cmdstan_model(file, include_paths = c(file.path(path_stan, 'stan_include')))

fit <- mod$generate_quantities(
  fitted_params = draws_tmp,
  data = data_list_sub,
  seed = 8472
)
draws_resid <- fit$draws()

df_rec <- data.frame(R = data_used_sub$ClstD_km, eq = eq_sub,
                     M = mageq[eq_sub],
                     VS = data_used_sub$Vs30_Selected_for_Analysis_m_s,
                     lnD75 = data_used_sub$D_75,
                     #resid_1 = colMeans(subset(as_draws_matrix(draws_resid),
                     #                             variable=c('^resid_1'), regex=TRUE)),
                     #resid_2 = colMeans(subset(as_draws_matrix(draws_resid),
                     #                          variable=c('^resid_2'), regex=TRUE)),
                     resid = colMeans(subset(as_draws_matrix(draws_resid),
                                             variable=c('^resid\\['), regex=TRUE))
)
resid_mean <- mean(df_rec$resid)

model <- "Mixed"
xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- "residual"
fl <- loess(resid ~ R, df_rec)
df_rec$smooth <- predict(fl)
p <- ggplot(df_rec) +
  geom_point(aes(x = R, y = resid), color = 'gray') +
  geom_line(aes(x = R, y = smooth)) +
  geom_hline(yintercept = resid_mean, linetype = "dashed") +
  labs(x = xlab, y = ylab, title = model) +
  theme(
    legend.position = "none"
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_residNGA_Rrup_%s.pdf', model)), p,
       width = wid, height = asp * wid)

xlab <- expression(atop(paste(R[RUP] (km))))
fl <- loess(resid ~ R, df_rec)
df_rec$smooth <- predict(fl)
p <- ggplot(df_rec) +
  geom_point(aes(x = R, y = resid), color = 'gray') +
  geom_line(aes(x = R, y = smooth)) +
  geom_hline(yintercept = resid_mean, linetype = "dashed") +
  labs(x = xlab, y = ylab, title = model) +
  theme(
    legend.position = "none"
  ) +
  scale_x_log10()
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_residNGA_lnRrup_%s.pdf', model)), p,
       width = wid, height = asp * wid)

xlab <- "M"
fl <- loess(resid ~ M, df_rec)
df_rec$smooth <- predict(fl)
p <- ggplot(df_rec) +
  geom_point(aes(x = M, y = resid), color = 'gray') +
  geom_line(aes(x = M, y = smooth)) +
  geom_hline(yintercept = resid_mean, linetype = "dashed") +
  labs(x = xlab, y = ylab, title = model) +
  theme(
    legend.position = "none"
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_residNGA_M_%s.pdf', model)), p,
       width = wid, height = asp * wid)

xlab <- expression(atop(paste(V[S30] (m/s))))
fl <- loess(resid ~ VS, df_rec)
df_rec$smooth <- predict(fl)
p <- ggplot(df_rec) +
  geom_point(aes(x = VS, y = resid), color = 'gray') +
  geom_line(aes(x = VS, y = smooth)) +
  geom_hline(yintercept = resid_mean, linetype = "dashed") +
  labs(x = xlab, y = ylab, title = model) +
  theme(
    legend.position = "none"
  ) +
  scale_x_log10()
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_residNGA_Vs30_%s.pdf', model)), p,
       width = wid, height = asp * wid)

####################################
####################################
####################################

psis_gam_source <- loo_gam_source$psis_object
lw_g <- weights(psis_gam_source) # normalized log weights

psis_logn_source <- loo_logn_source$psis_object
lw_ln <- weights(psis_logn_source) # normalized log weights

psis_mix_source <- loo_mix_source$psis_object
lw_m <- weights(psis_mix_source) # normalized log weights

# marginal predictive check using LOO probability integral transform
ppc_loo_pit_overlay(y, y_rep_gam, lw = lw_g)
ppc_loo_pit_overlay(y, y_rep_logn, lw = lw_ln)
ppc_loo_pit_overlay(y, y_rep_mix, lw = lw_m)


ppc_loo_pit_qq(y, y_rep_gam, lw = lw_g)
ppc_loo_pit_qq(y, y_rep_logn, lw = lw_ln)
ppc_loo_pit_qq(y, y_rep_mix, lw = lw_m)





