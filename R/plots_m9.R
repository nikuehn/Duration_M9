library(tidyverse)
library(patchwork)

wid <- 8
asp <- 0.8

theme_set(theme_bw() + theme(panel.grid.minor = element_blank(),
                             axis.title = element_text(size = 30),
                             axis.text = element_text(size = 20),
                             plot.title = element_text(size = 30),
                             legend.text = element_text(size = 20),
                             legend.title = element_text(size = 20),
                             legend.key.width = unit(1, "cm"),
                             legend.box.background = element_rect(colour = "black")
                             ))


# Plot of duration vs. distance
# fl <- loess(exp(lnD75) ~ RCD, data_used)
# pfl <- predict(fl, se = TRUE)
# data_used$smooth<- pfl$fit
# data_used$smooth_q1 <- pfl$fit - qt(0.975,pfl$df)*pfl$se
# data_used$smooth_q2 <- pfl$fit + qt(0.975,pfl$df)*pfl$se

l <- msir::loess.sd(data_used$RCD, exp(data_used$lnD75), nsigma = 1.96)
l_df <- data.frame(x = l$x, smooth = l$y, lower = l$lower, upper = l$upper)

xlab <- expression(atop(paste(R[RUP]," (km)")))
ylab <- expression(atop(paste(D[5-75])))

p <- ggplot(data_used) +
  geom_point(mapping = aes(x = RCD, y = exp(lnD75), color = `# RealizationID`)) +
  geom_line(data = l_df, aes(x = x, y = smooth), color = 'black') +
  #geom_line(data = l_df, aes(x = x, y = lower), color = 'black', linetype = 'dashed') +
  #geom_line(data = l_df, aes(x = x, y = upper), color = 'black', linetype = 'dashed') +
  #scale_x_continuous(trans='log10', limits = c(1,500)) +
  #scale_y_continuous(trans='log10', limits = c(1,300)) +
  labs(x = xlab, y = ylab) +
  theme(
    legend.position = "none"
  ) #+
  annotation_logticks()

ggsave(file.path(path_base,'PLOTS_M9','plot_D_Rrup.pdf'), p,
       width = wid, height = asp * wid)

p <- ggplot(data_used) +
  geom_point(mapping = aes(x = RCD, y = exp(lnD75), color = `# RealizationID`)) +
  geom_line(data = l_df, aes(x = x, y = smooth), color = 'black') +
  #geom_line(data = l_df, aes(x = x, y = lower), color = 'black', linetype = 'dashed') +
  #geom_line(data = l_df, aes(x = x, y = upper), color = 'black', linetype = 'dashed') +
  scale_x_continuous(trans='log10', limits = c(1,500)) +
  scale_y_continuous(trans='log10', limits = c(1,300)) +
  labs(x = xlab, y = ylab) +
  theme(
    legend.position = "none"
  ) +
  annotation_logticks()

ggsave(file.path(path_base,'PLOTS_M9','plot_lnD_lnRrup.pdf'), p,
       width = wid, height = asp * wid)

#########################
# Results
model <- 'gmm_dsource_hier_Rlh'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat10-40_%s.RDS', model)))
draws_gam <- fit$draws()
fit$diagnostic_summary()
loo_gam_source <- fit$loo()

model <- 'gmm_dsource_hier_Rlh_logn'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat10-40_%s.RDS', model)))
draws_logn <- fit$draws()
fit$diagnostic_summary()
loo_logn_source <- fit$loo()

model <- 'gmm_dsource_hier_Rlh_mix'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat10-40_%s.RDS', model)))
draws_mix <- fit$draws()
fit$diagnostic_summary()
loo_mix_source <- fit$loo()


#########################################
param <- "r1"
xlab <- expression(atop(paste(r[1])))
df_plot <- data.frame(gamma = as_draws_matrix(subset_draws(draws_gam, variable = param)),
                      logn = as_draws_matrix(subset_draws(draws_logn, variable = param)),
                      mix = as_draws_matrix(subset_draws(draws_mix, variable = param))
)
names(df_plot) <- c("gamma", "logn", "mix")
p <- df_plot %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  labs(x = xlab, y = "density") +
  scale_x_continuous(limits = c(0.4,0.6)) +
  theme(
    legend.position = c(0.8, 0.8),
    legend.box.background = element_rect(colour = "black")
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE,title=NULL)) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plot_dens_%s.pdf', param)), p,
       width = wid, height = asp * wid)

#########################################
param <- "r2"
xlab <- expression(atop(paste(r[2])))
df_plot <- data.frame(gamma = as_draws_matrix(subset_draws(draws_gam, variable = param)),
                      logn = as_draws_matrix(subset_draws(draws_logn, variable = param)),
                      mix = as_draws_matrix(subset_draws(draws_mix, variable = param))
)
names(df_plot) <- c("gamma", "logn", "mix")
p <- df_plot %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  labs(x = xlab, y = "density") +
  #scale_x_continuous(limits = c(0.4,0.6)) +
  theme(
    legend.position = "none"
  ) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plot_dens_%s.pdf', param)), p,
       width = wid, height = asp * wid)

#########################################
param <- "R_b"
xlab <- expression(atop(paste(R[b])))
df_plot <- data.frame(gamma = as_draws_matrix(subset_draws(draws_gam, variable = param)),
                      logn = as_draws_matrix(subset_draws(draws_logn, variable = param)),
                      mix = as_draws_matrix(subset_draws(draws_mix, variable = param))
)
names(df_plot) <- c("gamma", "logn", "mix")
p <- df_plot %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  labs(x = xlab, y = "density") +
  #scale_x_continuous(limits = c(0.4,0.6)) +
  theme(
    legend.position = "none"
  ) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plot_dens_%s.pdf', param)), p,
       width = wid, height = asp * wid)

#########################################
param <- "mu_ln_c1"
xlab <- expression(atop(paste(exp(mu[lnD[S]]))))
df_plot <- data.frame(gamma = exp(as_draws_matrix(subset_draws(draws_gam, variable = param))),
                      logn = exp(as_draws_matrix(subset_draws(draws_logn, variable = param))),
                      mix = as_draws_matrix(subset_draws(draws_mix, variable = "c1"))
)
names(df_plot) <- c("gamma", "logn", "mix")
p <- df_plot %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  labs(x = xlab, y = "density") +
  #scale_x_continuous(limits = c(0.4,0.6)) +
  theme(
    legend.position = "none"
  ) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plot_dens_%s.pdf', param)), p,
       width = wid, height = asp * wid)


#########################################
param <- "sigma_ln_c1"
xlab <- expression(atop(paste(sigma[S])))
df_plot <- data.frame(gamma = as_draws_matrix(subset_draws(draws_gam, variable = param)),
                      logn = as_draws_matrix(subset_draws(draws_logn, variable = param)),
                      mix = as_draws_matrix(subset_draws(draws_mix, variable = "sigma_source"))
)
names(df_plot) <- c("gamma", "logn", "mix")
p <- df_plot %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  labs(x = xlab, y = "density") +
  scale_x_continuous(limits = c(0.,0.3)) +
  theme(
    legend.position = "none"
  ) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plot_dens_%s.pdf', param)), p,
       width = wid, height = asp * wid)

#########################################
param <- "sigma_path"
xlab <- expression(atop(paste(sigma[P])))
df_plot <- data.frame(
  logn = as_draws_matrix(subset_draws(draws_logn, variable = param)),
  mix = as_draws_matrix(subset_draws(draws_mix, variable = param))
)
names(df_plot) <- c("logn", "mix")
p <- df_plot %>%
  pivot_longer(cols = c("logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  labs(x = xlab, y = "density") +
  scale_x_continuous(limits = c(0.45,0.55)) +
  theme(
    legend.position = "none"
  ) +
  scale_color_manual(values=c("red", "blue"),
                     labels=c(
                       "Lognormal",
                       "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plot_dens_%s.pdf', param)), p,
       width = wid, height = asp * wid)

#########################################
param <- "shape_path"
xlab <- expression(atop(paste(alpha[P])))
df_plot <- data.frame(
  gamma = as_draws_matrix(subset_draws(draws_gam, variable = param))
)
names(df_plot) <- c("gamma")
p <- df_plot %>%
  pivot_longer(cols = c("gamma"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  labs(x = xlab, y = "density") +
  #scale_x_continuous(limits = c(0.45,0.55)) +
  theme(
    legend.position = "none"
  ) +
  scale_color_manual(values=c("black"),
                     labels=c(
                       "Gamma"
                       )
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plot_dens_%s.pdf', param)), p,
       width = wid, height = asp * wid)



#########################################
dist_pred <- seq(0,300,by = 10)
n_eq

coeffs_gam <- colMeans(subset(as_draws_matrix(draws_gam), variable=c("r1", "r2", "R_b", "mu_ln_c1", "c1"), regex=FALSE))
coeffs_logn <- colMeans(subset(as_draws_matrix(draws_logn), variable=c("r1", "r2", "R_b", "mu_ln_c1", "c1"), regex=FALSE))
coeffs_mix <- colMeans(subset(as_draws_matrix(draws_mix), variable=c("r1", "r2", "R_b", "c1", "eqterm"), regex=FALSE))

pred_gam <- matrix(nrow = length(dist_pred), ncol =n_eq)
pred_logn <- matrix(nrow = length(dist_pred), ncol =n_eq)
pred_mix <- matrix(nrow = length(dist_pred), ncol =n_eq)

for(k in 1:n_eq) {
  pred_gam[,k] <- func_lh(dist_pred, c(coeffs_gam[4 + k] + coeffs_gam[1] * coeffs_gam[3], coeffs_gam), 
                          mb = coeffs_gam[3], delta = 1)
  pred_logn[,k] <- func_lh(dist_pred, c(coeffs_logn[4 + k] + coeffs_logn[1] * coeffs_logn[3], coeffs_logn), 
                           mb = coeffs_logn[3], delta = 1)
  pred_mix[,k] <- exp(log(func_lh(dist_pred, c(coeffs_mix[4] + coeffs_mix[1] * coeffs_mix[3], coeffs_mix), 
                                  mb = coeffs_mix[3], delta = 1)) + coeffs_mix[4 + k])
}

# Gamma
model <- "Gamma"
xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- expression(atop(paste(D[5-75] (s))))
df_1 <- as.data.frame(pred_gam)
df_1$R <- dist_pred
p <- df_1 %>%
  pivot_longer(!R, names_to = "eq", values_to = "var") %>%
  ggplot() +
  geom_line(aes(x = R, y = var, color = eq)) +
  labs(x = xlab, y = ylab, title = model) +
  ylim(0,150) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20)
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plot_pred_eq_%s.pdf', model)), p,
       width = wid, height = asp * wid)

# logn
model <- "Lognormal"
xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- expression(atop(paste(D[5-75] (s))))
df_1 <- as.data.frame(pred_logn)
df_1$R <- dist_pred
p <- df_1 %>%
  pivot_longer(!R, names_to = "eq", values_to = "var") %>%
  ggplot() +
  geom_line(aes(x = R, y = var, color = eq)) +
  labs(x = xlab, y = ylab, title = model) +
  ylim(0,150) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20)
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plot_pred_eq_%s.pdf', model)), p,
       width = wid, height = asp * wid)

# mixed
model <- "Mixed"
xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- expression(atop(paste(D[5-75] (s))))
df_1 <- as.data.frame(pred_mix)
df_1$R <- dist_pred
p <- df_1 %>%
  pivot_longer(!R, names_to = "eq", values_to = "var") %>%
  ggplot() +
  geom_line(aes(x = R, y = var, color = eq)) +
  labs(x = xlab, y = ylab, title = model) +
  ylim(0,150) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20)
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plot_pred_eq_%s.pdf', model)), p,
       width = wid, height = asp * wid)


#########################################
# Residuals
xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- "residual"

# Gamma
model = "Gamma"
df_rec <- data.frame(R = data_used$RCD, eq = eq,
                     lnD75 = data_used$lnD75,
                     resid_mean = colMeans(subset(as_draws_matrix(draws_gam), variable=c('^resid'), regex=TRUE))
)
fl <- loess(resid_mean ~ log(R), df_rec)
df_rec$smooth <- predict(fl)
p <- df_rec %>%
  ggplot() +
  geom_point(aes(x = R, y = resid_mean, color = as.factor(eq))) +
  geom_line(data = df_rec, mapping = aes(x = R, y = smooth)) +
  labs(x = xlab, y = ylab, title = model) +
  ylim(-2.5,2.5) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20)
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plot_resid_%s.pdf', model)), p,
       width = wid, height = asp * wid)

# Lognormal
model = "Lognormal"
df_rec <- data.frame(R = data_used$RCD, eq = eq,
                     lnD75 = data_used$lnD75,
                     resid_mean = colMeans(subset(as_draws_matrix(draws_logn), variable=c('^resid'), regex=TRUE))
)
fl <- loess(resid_mean ~ log(R), df_rec)
df_rec$smooth <- predict(fl)
p <- df_rec %>%
  ggplot() +
  geom_point(aes(x = R, y = resid_mean, color = as.factor(eq))) +
  geom_line(data = df_rec, mapping = aes(x = R, y = smooth)) +
  labs(x = xlab, y = ylab, title = model) +
  ylim(-2.5,2.5) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20)
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plot_resid_%s.pdf', model)), p,
       width = wid, height = asp * wid)


# Mixed
model = "Mixed"
df_rec <- data.frame(R = data_used$RCD, eq = eq,
                     lnD75 = data_used$lnD75,
                     resid_mean = colMeans(subset(as_draws_matrix(draws_mix), variable=c('^resid'), regex=TRUE))
)
fl <- loess(resid_mean ~ log(R), df_rec)
df_rec$smooth <- predict(fl)
p <- df_rec %>%
  ggplot() +
  geom_point(aes(x = R, y = resid_mean, color = as.factor(eq))) +
  geom_line(data = df_rec, mapping = aes(x = R, y = smooth)) +
  labs(x = xlab, y = ylab, title = model) +
  ylim(-2.5,2.5) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20)
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plot_resid_%s.pdf', model)), p,
       width = wid, height = asp * wid)



#########################################
# distribuion prediction
post_gam <- subset(as_draws_matrix(draws_gam), variable=c("r1", "r2", "R_b",
                                                          "mu_ln_c1", "sigma_ln_c1", "shape_path"), regex=FALSE)
post_logn <- subset(as_draws_matrix(draws_logn), variable=c("r1", "r2", "R_b",
                                                            "mu_ln_c1", "sigma_ln_c1", "sigma_path"), regex=FALSE)
post_mix <- subset(as_draws_matrix(draws_mix), variable=c("r1", "r2", "R_b", 
                                                          'c1', "sigma_source", 'sigma_path'), regex=FALSE)

write.csv(post_gam, file = file.path(path_stan, 'RESULTS', 'posterior_gamma.csv'))
write.csv(post_logn, file = file.path(path_stan, 'RESULTS', 'posterior_logn.csv'))
write.csv(post_mix, file = file.path(path_stan, 'RESULTS', 'posterior_mix.csv'))


df_pred <- read.csv(file.path(path_stan, 'RESULTS', 'predict_distri.csv'))

xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- expression(atop(paste(D[5-75] (s))))
p <- ggplot(df_pred) +
  geom_point(data_used, mapping = aes(x = RCD, y = exp(lnD75)), color = 'gray') +
  geom_line(mapping = aes(x = R, y = Mean_gam)) +
  geom_line(mapping = aes(x = R, y = Med_gam), linetype = "dotted") +
  geom_line(mapping = aes(x = R, y = Q05_gam), linetype = "dashed") +
  geom_line(mapping = aes(x = R, y = Q95_gam), linetype = "dashed") +
  geom_line(mapping = aes(x = R, y = Mean_logn), color = 'red') +
  geom_line(mapping = aes(x = R, y = Med_logn), color = 'red', linetype = "dotted") +
  geom_line(mapping = aes(x = R, y = Q05_logn), color = 'red', linetype = "dashed") +
  geom_line(mapping = aes(x = R, y = Q95_logn), color = 'red', linetype = "dashed") +
  geom_line(mapping = aes(x = R, y = Mean_mix), color = 'blue') +
  geom_line(mapping = aes(x = R, y = Med_mix), color = 'blue', linetype = "dotted") +
  geom_line(mapping = aes(x = R, y = Q05_mix), color = 'blue', linetype = "dashed") +
  geom_line(mapping = aes(x = R, y = Q95_mix), color = 'blue', linetype = "dashed") +
  labs(x = xlab, y = ylab)
ggsave(file.path(path_base,'PLOTS_M9', 'plot_predR_distri.pdf'), p,
       width = wid, height = asp * wid)


#######################################
### Posterior predictive checks
idx <- 1:n_rec

model <- "Gamma"
y_rep <- subset(as_draws_matrix(draws_gam), variable=c('^y_'), regex=TRUE)[,idx]
y <- exp(data_used$lnD75[idx])

ppc_stat(y, y_rep) +
  labs(title = model) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20)
  )
ppc_stat_2d(y, y_rep, stat = c("q25", "q75")) +
  labs(title = model) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20)
  )


###########################################################################
###########################################################################
###########################################################################
###########################################################################
########################## Models wih shapeR/sigmaR3
########################## Results
model <- 'gmm_dsource_hier_Rlh_shapeR'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat10-40_%s.RDS', model)))
draws_gam <- fit$draws()
fit$diagnostic_summary()
loo_gam_source <- fit$loo(save_psis = TRUE)

model <- 'gmm_dsource_hier_Rlh_logn_sigmaR3'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat10-40_%s.RDS', model)))
draws_logn <- fit$draws()
fit$diagnostic_summary()
loo_logn_source <- fit$loo(save_psis = TRUE)

model <- 'gmm_dsource_hier_Rlh_mix_sigmaR3'
fit <- readRDS(file.path(path_stan, 'RESULTS', sprintf('fit_m9_Lat10-40_%s.RDS', model)))
draws_mix <- fit$draws()
fit$diagnostic_summary()
loo_mix_source <- fit$loo(save_psis = TRUE)

loo::loo_compare(loo_gam_source, loo_logn_source, loo_mix_source)


#########################################
param <- "r1"
xlab <- expression(atop(paste(r[1])))
df_plot <- data.frame(gamma = as_draws_matrix(subset_draws(draws_gam, variable = param)),
                      logn = as_draws_matrix(subset_draws(draws_logn, variable = param)),
                      mix = as_draws_matrix(subset_draws(draws_mix, variable = param))
)
names(df_plot) <- c("gamma", "logn", "mix")
p <- df_plot %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model), key_glyph = draw_key_path) +
  labs(x = xlab, y = "density") +
  scale_x_continuous(limits = c(0.4,0.65)) +
  theme(
    legend.position = c(0.77, 0.8)
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE,title=NULL)) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_dens_%s.pdf', param)), p,
       width = wid, height = asp * wid)

#########################################
param <- "r2"
xlab <- expression(atop(paste(r[2])))
df_plot <- data.frame(gamma = as_draws_matrix(subset_draws(draws_gam, variable = param)),
                      logn = as_draws_matrix(subset_draws(draws_logn, variable = param)),
                      mix = as_draws_matrix(subset_draws(draws_mix, variable = param))
)
names(df_plot) <- c("gamma", "logn", "mix")
p <- df_plot %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  labs(x = xlab, y = "density") +
  scale_x_continuous(limits = c(0.24,0.3)) +
  theme(
    legend.position = "none"
  ) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_dens_%s.pdf', param)), p,
       width = wid, height = asp * wid)

#########################################
param <- "R_b"
xlab <- expression(atop(paste(R[b])))
df_plot <- data.frame(gamma = as_draws_matrix(subset_draws(draws_gam, variable = param)),
                      logn = as_draws_matrix(subset_draws(draws_logn, variable = param)),
                      mix = as_draws_matrix(subset_draws(draws_mix, variable = param))
)
names(df_plot) <- c("gamma", "logn", "mix")
p <- df_plot %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  labs(x = xlab, y = "density") +
  scale_x_continuous(limits = c(50,90)) +
  theme(
    legend.position = "none"
  ) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_dens_%s.pdf', param)), p,
       width = wid, height = asp * wid)

#########################################
param <- "mu_ln_c1"
xlab <- expression(atop(paste(exp(mu[lnD[S]]))))
df_plot <- data.frame(gamma = exp(as_draws_matrix(subset_draws(draws_gam, variable = param))),
                      logn = exp(as_draws_matrix(subset_draws(draws_logn, variable = param))),
                      mix = as_draws_matrix(subset_draws(draws_mix, variable = "c1"))
)
names(df_plot) <- c("gamma", "logn", "mix")
p <- df_plot %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  labs(x = xlab, y = "density") +
  #scale_x_continuous(limits = c(0.4,0.6)) +
  theme(
    legend.position = "none"
  ) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_dens_%s.pdf', param)), p,
       width = wid, height = asp * wid)


xlab <- "average source duration"
df_plot <- rbind(setNames(data.frame(value = exp(as_draws_matrix(subset_draws(draws_gam, variable = param))),
                                     model = "gamma", par = "mean"), c("value", "model", "par")),
                 setNames(data.frame(value = exp(as_draws_matrix(subset_draws(draws_logn, variable = param))),
                                     model = "logn", par = "median"), c("value", "model", "par")),
                 setNames(data.frame(value = as_draws_matrix(subset_draws(draws_mix, variable = "c1")),
                                     model = "mix", par = "median"), c("value", "model", "par")),
                 setNames(data.frame(value = exp(as_draws_matrix(subset_draws(
                   mutate_variables(coeffs_logn, mean_ln_c1 = mu_ln_c1 + 0.5 * (rp_0 + rp_1 / (1 + exp(-rp_s * (ref_dist - R_s))))^2),
                   variable = "mean_ln_c1"))),
                   model = "logn", par = "mean"), c("value", "model", "par")),
                 setNames(data.frame(value = exp(as_draws_matrix(subset_draws(
                   mutate_variables(coeffs_mix, mean_ln_c1 = log(c1) + 0.5 * (rp_0 + rp_1 / (1 + exp(-rp_s * (ref_dist - R_s))))^2),
                   variable = "mean_ln_c1"))),
                   model = "mix", par = "mean"), c("value", "model", "par"))
)

p <- df_plot %>%
  ggplot() +
  geom_density(aes(x = value, color = model, linetype = par), key_glyph = draw_key_path) +
  labs(x = xlab, y = "density") +
  #scale_x_continuous(limits = c(0.4,0.6)) +
  theme(
    #legend.position = "none"
    legend.position = c(0.82,0.76)
  ) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed"),
                     name = "Model"
  ) +
  scale_linetype_manual(values=c("solid",  "dashed"),
                        labels=c("Mean", "Median")
  ) +
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE,title=NULL),
         color = guide_legend(nrow = 3, byrow = TRUE,title=NULL))
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_dens_%s.pdf', "c1")), p,
       width = wid, height = asp * wid)


#########################################
param <- "sigma_ln_c1"
xlab <- expression(atop(paste(sigma[S])))
df_plot <- data.frame(gamma = as_draws_matrix(subset_draws(draws_gam, variable = param)),
                      logn = as_draws_matrix(subset_draws(draws_logn, variable = param)),
                      mix = as_draws_matrix(subset_draws(draws_mix, variable = "sigma_source"))
)
names(df_plot) <- c("gamma", "logn", "mix")
p <- df_plot %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  labs(x = xlab, y = "density") +
  scale_x_continuous(limits = c(0.,0.45)) +
  theme(
    legend.position = "none"
  ) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_dens_%s.pdf', param)), p,
       width = wid, height = asp * wid)



#########################################
dist_pred <- seq(0,300,by = 10)
n_eq

coeffs_gam <- subset(as_draws_matrix(draws_gam), variable=c("rp_0", "rp_1"), regex=FALSE)
coeffs_logn <- subset(as_draws_matrix(draws_logn), variable=c("rp_0", "rp_1", "rp_s", "R_s"), regex=FALSE)
coeffs_mix <- subset(as_draws_matrix(draws_mix), variable=c("rp_0", "rp_1", "rp_s", "R_s", "sigma_source"), regex=FALSE)

shape_gam <- matrix(nrow = length(dist_pred), ncol =3)
shape_logn <- matrix(nrow = length(dist_pred), ncol =3)
shape_mix <- matrix(nrow = length(dist_pred), ncol =3)
shape_mix_total <- matrix(nrow = length(dist_pred), ncol =3)
k <- 1
for(k in 1:length(dist_pred)) {
  tmp <- mutate_variables(coeffs_gam, shape = rp_0 + rp_1 * dist_pred[k])
  shape_gam[k,] <- c(mean(subset_draws(tmp, variable = "shape")),
                     quantile(subset_draws(tmp, variable = "shape"), 0.05),
                     quantile(subset_draws(tmp, variable = "shape"), 0.95))
  
  tmp <- mutate_variables(coeffs_logn, shape = rp_0 + rp_1 / (1 + exp(-rp_s * (dist_pred[k] - R_s))))
  shape_logn[k,] <- c(mean(subset_draws(tmp, variable = "shape")),
                     quantile(subset_draws(tmp, variable = "shape"), 0.05),
                     quantile(subset_draws(tmp, variable = "shape"), 0.95))
  
  tmp <- mutate_variables(coeffs_mix, shape = rp_0 + rp_1 / (1 + exp(-rp_s * (dist_pred[k] - R_s))))
  shape_mix[k,] <- c(mean(subset_draws(tmp, variable = "shape")),
                      quantile(subset_draws(tmp, variable = "shape"), 0.05),
                      quantile(subset_draws(tmp, variable = "shape"), 0.95))
  
  tmp <- mutate_variables(coeffs_mix, 
                          shape = sqrt((rp_0 + rp_1 / (1 + exp(-rp_s * (dist_pred[k] - R_s))))^2 + sigma_source^2))
  shape_mix_total[k,] <- c(mean(subset_draws(tmp, variable = "shape")),
                     quantile(subset_draws(tmp, variable = "shape"), 0.05),
                     quantile(subset_draws(tmp, variable = "shape"), 0.95))
}

model <- "Gamma"
xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- expression(atop(paste(alpha[P])))
shape_gam <- as.data.frame(shape_gam)
shape_gam$R <- dist_pred
p <- shape_gam %>%
  pivot_longer(!R, names_to = "frac", values_to = "var") %>%
  ggplot() +
  geom_line(aes(x = R, y = var, linetype = frac)) +
  labs(x = xlab, y = ylab) +
  ylim(0,25) +
  theme(
    legend.position = c(0.3,0.77),
    plot.title = element_text(size = 20)
  ) +
  scale_linetype_manual(values=c("solid", "dotted", "dashed"),
                     labels=c("Mean", "5% fractile", "95% fractile")
  ) +
  guides(linetype = guide_legend(nrow = 3, byrow = TRUE,title=NULL))
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_shapeP_%s.pdf', model)), p,
       width = wid, height = asp * wid)


shape_mix <- as.data.frame(shape_mix)
shape_mix$R <- dist_pred
shape_mix <- shape_mix %>%
  pivot_longer(!R, names_to = "frac", values_to = "var")
shape_mix$model <- "Mixed"

shape_logn <- as.data.frame(shape_logn)
shape_logn$R <- dist_pred
shape_logn <- shape_logn %>%
  pivot_longer(!R, names_to = "frac", values_to = "var")
shape_logn$model <- "Lognormal"

ylab <- expression(atop(paste(sigma[P])))
df_plot <- rbind(shape_logn, shape_mix)
p <- ggplot(df_plot) +
  geom_line(aes(x = R, y = var, linetype = frac, color = model), key_glyph = draw_key_path) +
  labs(x = xlab, y = ylab) +
  ylim(0,1) +
  theme(
    legend.position = c(0.7,0.8)
  ) +
  scale_linetype_manual(values=c("solid", "dotted", "dashed"),
                        labels=c("Mean", "5% fractile", "95% fractile")
  ) +
  scale_color_manual(values=c("red", "blue"),
                        labels=c("Lognormal", "Mixed")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE,title=NULL),
         linetype = "none")
ggsave(file.path(path_base,'PLOTS_M9','plotR_sigmaP.pdf'), p,
       width = wid, height = asp * wid)

shape_mix_total <- as.data.frame(shape_mix_total)
shape_mix_total$R <- dist_pred
shape_mix_total <- shape_mix_total %>%
  pivot_longer(!R, names_to = "frac", values_to = "var")
shape_mix_total$model <- "Mix_total"

ylab <- expression(atop(paste(sigma)))
df_plot <- rbind(shape_mix, shape_mix_total)
p <- ggplot(df_plot) +
  geom_line(aes(x = R, y = var, color = model, linetype = frac), key_glyph = draw_key_path) +
  labs(x = xlab, y = ylab) +
  ylim(0,1) +
  theme(
    legend.position = c(0.8,0.8)
  ) +
  scale_linetype_manual(values=c("solid", "dotted", "dashed"),
                        labels=c("Mean", "5% fractile", "95% fractile")
  ) +
  scale_color_manual(values=c("cyan", "blue"),
                     labels=c(expression(atop(paste(sigma[T]))), expression(atop(paste(sigma[P]))))
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE,title=NULL),
         linetype = "none")
ggsave(file.path(path_base,'PLOTS_M9','plotR_sigmaPT_mix.pdf'), p,
       width = wid, height = asp * wid)
  

#########################################
dist_pred <- seq(0,300,by = 10)
n_eq

coeffs_gam <- colMeans(subset(as_draws_matrix(draws_gam), variable=c("r1", "r2", "R_b", "mu_ln_c1", "c1"), regex=FALSE))
coeffs_logn <- colMeans(subset(as_draws_matrix(draws_logn), variable=c("r1", "r2", "R_b", "mu_ln_c1", "c1"), regex=FALSE))
coeffs_mix <- colMeans(subset(as_draws_matrix(draws_mix), variable=c("r1", "r2", "R_b", "c1", "eqterm"), regex=FALSE))

pred_gam <- matrix(nrow = length(dist_pred), ncol =n_eq)
pred_logn <- matrix(nrow = length(dist_pred), ncol =n_eq)
pred_mix <- matrix(nrow = length(dist_pred), ncol =n_eq)

for(k in 1:n_eq) {
  pred_gam[,k] <- func_lh(dist_pred, c(coeffs_gam[4 + k] + coeffs_gam[1] * coeffs_gam[3], coeffs_gam), 
                          mb = coeffs_gam[3], delta = 1)
  pred_logn[,k] <- func_lh(dist_pred, c(coeffs_logn[4 + k] + coeffs_logn[1] * coeffs_logn[3], coeffs_logn), 
                           mb = coeffs_logn[3], delta = 1)
  pred_mix[,k] <- exp(log(func_lh(dist_pred, c(coeffs_mix[4] + coeffs_mix[1] * coeffs_mix[3], coeffs_mix), 
                                  mb = coeffs_mix[3], delta = 1)) + coeffs_mix[4 + k])
}

# Gamma
model <- "Gamma"
xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- expression(atop(paste(D[5-75] (s))))
df_1 <- as.data.frame(pred_gam)
df_1$R <- dist_pred
p <- df_1 %>%
  pivot_longer(!R, names_to = "eq", values_to = "var") %>%
  ggplot() +
  geom_line(aes(x = R, y = var, color = eq)) +
  labs(x = xlab, y = ylab, title = model) +
  ylim(0,150) +
  theme(
    legend.position = "none"
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_pred_eq_%s.pdf', model)), p,
       width = wid, height = asp * wid)

# logn
model <- "Lognormal"
xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- expression(atop(paste(D[5-75] (s))))
df_1 <- as.data.frame(pred_logn)
df_1$R <- dist_pred
p <- df_1 %>%
  pivot_longer(!R, names_to = "eq", values_to = "var") %>%
  ggplot() +
  geom_line(aes(x = R, y = var, color = eq)) +
  labs(x = xlab, y = ylab, title = model) +
  ylim(0,150) +
  theme(
    legend.position = "none"
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_pred_eq_%s.pdf', model)), p,
       width = wid, height = asp * wid)

# mixed
model <- "Mixed"
xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- expression(atop(paste(D[5-75] (s))))
df_1 <- as.data.frame(pred_mix)
df_1$R <- dist_pred
p <- df_1 %>%
  pivot_longer(!R, names_to = "eq", values_to = "var") %>%
  ggplot() +
  geom_line(aes(x = R, y = var, color = eq)) +
  labs(x = xlab, y = ylab, title = model) +
  ylim(0,150) +
  theme(
    legend.position = "none"
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_pred_eq_%s.pdf', model)), p,
       width = wid, height = asp * wid)


#########################################
# Residuals
xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- "residual"

# Gamma
model = "Gamma"
df_rec <- data.frame(R = data_used$RCD, eq = eq,
                     lnD75 = data_used$lnD75,
                     resid_mean = colMeans(subset(as_draws_matrix(draws_gam), variable=c('^resid'), regex=TRUE))
)
fl <- loess(resid_mean ~ log(R), df_rec)
df_rec$smooth <- predict(fl)
p <- df_rec %>%
  ggplot() +
  geom_point(aes(x = R, y = resid_mean, color = as.factor(eq))) +
  geom_line(data = df_rec, mapping = aes(x = R, y = smooth)) +
  labs(x = xlab, y = ylab, title = model) +
  ylim(-2.5,2.5) +
  theme(
    legend.position = "none"
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_resid_%s.pdf', model)), p,
       width = wid, height = asp * wid)

# Lognormal
model = "Lognormal"
df_rec <- data.frame(R = data_used$RCD, eq = eq,
                     lnD75 = data_used$lnD75,
                     resid_mean = colMeans(subset(as_draws_matrix(draws_logn), variable=c('^resid'), regex=TRUE))
)
fl <- loess(resid_mean ~ log(R), df_rec)
df_rec$smooth <- predict(fl)
p <- df_rec %>%
  ggplot() +
  geom_point(aes(x = R, y = resid_mean, color = as.factor(eq))) +
  geom_line(data = df_rec, mapping = aes(x = R, y = smooth)) +
  labs(x = xlab, y = ylab, title = model) +
  ylim(-2.5,2.5) +
  theme(
    legend.position = "none"
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_resid_%s.pdf', model)), p,
       width = wid, height = asp * wid)


# Mixed
model = "Mixed"
df_rec <- data.frame(R = data_used$RCD, eq = eq,
                     lnD75 = data_used$lnD75,
                     resid_mean = colMeans(subset(as_draws_matrix(draws_mix), variable=c('^resid'), regex=TRUE))
)
fl <- loess(resid_mean ~ log(R), df_rec)
df_rec$smooth <- predict(fl)
p <- df_rec %>%
  ggplot() +
  geom_point(aes(x = R, y = resid_mean, color = as.factor(eq))) +
  geom_line(data = df_rec, mapping = aes(x = R, y = smooth)) +
  labs(x = xlab, y = ylab, title = model) +
  ylim(-2.5,2.5) +
  theme(
    legend.position = "none"
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_resid_%s.pdf', model)), p,
       width = wid, height = asp * wid)



#########################################
# distribuion prediction
post_gam <- subset(as_draws_matrix(draws_gam), variable=c("r1", "r2", "R_b",
                                                          "mu_ln_c1", "sigma_ln_c1", "rp_0", "rp_1"), regex=FALSE)
post_logn <- subset(as_draws_matrix(draws_logn), variable=c("r1", "r2", "R_b",
                                                            "mu_ln_c1", "sigma_ln_c1",
                                                            "rp_0", "rp_1", "rp_s", "R_s"), regex=FALSE)
post_mix <- subset(as_draws_matrix(draws_mix), variable=c("r1", "r2", "R_b", 
                                                          'c1', "sigma_source",
                                                          "rp_0", "rp_1", "rp_s", "R_s"), regex=FALSE)

write.csv(post_gam, file = file.path(path_stan, 'RESULTS', 'posterior_gamma_R.csv'))
write.csv(post_logn, file = file.path(path_stan, 'RESULTS', 'posterior_logn_R.csv'))
write.csv(post_mix, file = file.path(path_stan, 'RESULTS', 'posterior_mix_R.csv'))


df_pred <- read.csv(file.path(path_stan, 'RESULTS', 'predict_distri_R.csv'))
dim(df_pred)

xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- expression(atop(paste(D[5-75] (s))))
p <- ggplot(df_pred) +
  geom_point(data_used, mapping = aes(x = RCD, y = exp(lnD75)), color = 'gray') +
  geom_line(mapping = aes(x = R, y = Mean_gam)) +
  geom_line(mapping = aes(x = R, y = Med_gam), linetype = "dotted") +
  geom_line(mapping = aes(x = R, y = Q05_gam), linetype = "dashed") +
  geom_line(mapping = aes(x = R, y = Q95_gam), linetype = "dashed") +
  geom_line(mapping = aes(x = R, y = Mean_logn), color = 'red') +
  geom_line(mapping = aes(x = R, y = Med_logn), color = 'red', linetype = "dotted") +
  geom_line(mapping = aes(x = R, y = Q05_logn), color = 'red', linetype = "dashed") +
  geom_line(mapping = aes(x = R, y = Q95_logn), color = 'red', linetype = "dashed") +
  geom_line(mapping = aes(x = R, y = Mean_mix), color = 'blue') +
  geom_line(mapping = aes(x = R, y = Med_mix), color = 'blue', linetype = "dotted") +
  geom_line(mapping = aes(x = R, y = Q05_mix), color = 'blue', linetype = "dashed") +
  geom_line(mapping = aes(x = R, y = Q95_mix), color = 'blue', linetype = "dashed") +
  labs(x = xlab, y = ylab) #+
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA))
ggsave(file.path(path_base,'PLOTS_M9', 'plotR_pred_distri.pdf'), p,
       width = wid, height = asp * wid)

ggsave(file.path(path_base,'PLOTS_M9', 'plotR_pred_distri_log.pdf'), p + scale_x_log10(limits = c(5,300)) + scale_y_log10(),
       width = wid, height = asp * wid)


#-------------------------------------------
df_pred <- read.csv(file.path(path_stan, 'RESULTS', 'predict_distri_R_long.csv'))
dim(df_pred)
names(df_pred)

xlab <- expression(atop(paste(R[RUP] (km))))
ylab <- expression(atop(paste(D[5-75] (s))))
p <- ggplot(df_pred) +
  geom_point(data_used, mapping = aes(x = RCD, y = exp(lnD75)), color = 'gray') +
  geom_line(mapping = aes(x = R, y = value, color = model, linetype = par)) +
  labs(x = xlab, y = ylab) +
  theme(
    legend.position = c(0.35, 0.88)
  ) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma", "Lognormal", "Mixed")
  ) +
  scale_linetype_manual(values=c("solid",  "dotted", "dashed", "dotdash"),
                        labels=c("Mean", "Median", "5% fractile", "95% fractile")
  ) +
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE,title=NULL),
         #color = guide_legend(nrow = 3, byrow = TRUE,title=NULL)
         color = "none"
         )
ggsave(file.path(path_base,'PLOTS_M9', 'plotR_pred_distri2.pdf'), p,
       width = wid, height = asp * wid)

p <- ggplot(df_pred) +
  geom_point(data_used, mapping = aes(x = RCD, y = exp(lnD75)), color = 'gray') +
  geom_line(mapping = aes(x = R, y = value, color = model, linetype = par)) +
  labs(x = xlab, y = ylab) +
  theme(
    legend.position = c(0.8, 0.2)
  ) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma", "Lognormal", "Mixed")
  ) +
  scale_linetype_manual(values=c("solid",  "dotted", "dashed", "dotdash"),
                        labels=c("Mean", "Median", "5% fractile", "95% fractile")
  ) +
  guides(linetype = "none",
         color = guide_legend(nrow = 3, byrow = TRUE,title=NULL)
  ) +
  scale_x_log10(limits = c(5,300)) + scale_y_log10()
ggsave(file.path(path_base,'PLOTS_M9', 'plotR_pred_distri2_log.pdf'), p,
       width = wid, height = asp * wid)

#######################################
### Posterior predictive checks
idx <- 1:n_rec

q25 <- function(y) quantile(y, 0.05)
q75 <- function(y) quantile(y, 0.95)

y_rep_gam <- subset(as_draws_matrix(draws_gam), variable=c('^y_'), regex=TRUE)[,idx]
y_rep_logn <- subset(as_draws_matrix(draws_logn), variable=c('^y_'), regex=TRUE)[,idx]
y_rep_mix <- subset(as_draws_matrix(draws_mix), variable=c('^y_'), regex=TRUE)[,idx]
y <- exp(data_used$lnD75[idx])

### Median
par <- "Median"
xlab <- par
p <- data.frame(gamma = rowMedians(y_rep_gam),
           logn = rowMedians(y_rep_logn),
           mix = rowMedians(y_rep_mix)
) %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
    geom_density(aes(x = var, color = model), key_glyph = draw_key_path) +
  vline_at(median(y), size = 2) +
  labs(x = xlab, y = "density") +
  scale_x_continuous(limits = c(57,62), breaks=seq(57, 62, 1)) +
  theme(
    legend.position = c(0.77, 0.8)
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE,title=NULL)) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_PPC_%s.pdf', par)), p,
       width = wid, height = asp * wid)

### Mean
par <- "Mean"
xlab <- par
p <- data.frame(gamma = rowMeans(y_rep_gam),
                logn = rowMeans(y_rep_logn),
                mix = rowMeans(y_rep_mix)
) %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model), key_glyph = draw_key_path) +
  vline_at(mean(y), size = 2) +
  labs(x = xlab, y = "density") +
  scale_x_continuous(limits = c(61,64), breaks=seq(61, 64, 1)) +
  theme(
    legend.position = "none"
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE,title=NULL)) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_PPC_%s.pdf', par)), p,
       width = wid, height = asp * wid)

qu <- 0.95
par <- sprintf("Qu95")
xlab <- "95% fractile"
p <- data.frame(gamma = rowQuantiles(y_rep_gam, probs = qu),
           logn = rowQuantiles(y_rep_logn, probs = qu),
           mix = rowQuantiles(y_rep_mix, probs = qu)
) %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  vline_at(quantile(y, qu), size = 2) +
  labs(x = xlab, y = "density") + 
  scale_x_continuous(limits = c(119,131), breaks=seq(120, 130, 5)) +
  theme(
    legend.position = "none"
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE,title=NULL)) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_PPC_%s.pdf', par)), p,
       width = wid, height = asp * wid)

qu <- 0.05
par <- sprintf("Qu05")
xlab <- "5% fractile"
p <- data.frame(gamma = rowQuantiles(y_rep_gam, probs = qu),
                logn = rowQuantiles(y_rep_logn, probs = qu),
                mix = rowQuantiles(y_rep_mix, probs = qu)
) %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  vline_at(quantile(y, qu), size = 2) +
  labs(x = xlab, y = "density") + 
  scale_x_continuous(limits = c(12,15), breaks=seq(12, 15, 1)) +
  theme(
    legend.position = "none"
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE,title=NULL)) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_PPC_%s.pdf', par)), p,
       width = wid, height = asp * wid)

qu <- 0.25
par <- sprintf("Qu25")
xlab <- "25% fractile"
p <- data.frame(gamma = rowQuantiles(y_rep_gam, probs = qu),
                logn = rowQuantiles(y_rep_logn, probs = qu),
                mix = rowQuantiles(y_rep_mix, probs = qu)
) %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  vline_at(quantile(y, qu), size = 2) +
  labs(x = xlab, y = "density") + 
  scale_x_continuous(limits = c(32,37), breaks=seq(32, 37, 1)) +
  theme(
    legend.position = "none"
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE,title=NULL)) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_PPC_%s.pdf', par)), p,
       width = wid, height = asp * wid)


qu <- 0.75
par <- sprintf("Qu75")
xlab <- "75% fractile"
p <- data.frame(gamma = rowQuantiles(y_rep_gam, probs = qu),
                logn = rowQuantiles(y_rep_logn, probs = qu),
                mix = rowQuantiles(y_rep_mix, probs = qu)
) %>%
  pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
  ggplot() +
  geom_density(aes(x = var, color = model)) +
  vline_at(quantile(y, qu), size = 2) +
  labs(x = xlab, y = "density") + 
  scale_x_continuous(limits = c(83,88), breaks=seq(83, 88, 1)) +
  theme(
    legend.position = "none"
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE,title=NULL)) +
  scale_color_manual(values=c("black", "red", "blue"),
                     labels=c("Gamma",
                              "Lognormal",
                              "Mixed")
  )
ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_PPC_%s.pdf', par)), p,
       width = wid, height = asp * wid)

##################
# PPC with bins
bin_dist  <- seq(0,300, by = 25)
idx_bin <- as.numeric(mltools::bin_data(data_used$RCD, bins = bin_dist))
n_bin <- max(idx_bin)

k <- 1
for(k in 1:n_bin) {
  idx <- which(idx_bin == k)
  y_rep_gam <- subset(as_draws_matrix(draws_gam), variable=c('^y_'), regex=TRUE)[,idx]
  y_rep_logn <- subset(as_draws_matrix(draws_logn), variable=c('^y_'), regex=TRUE)[,idx]
  y_rep_mix <- subset(as_draws_matrix(draws_mix), variable=c('^y_'), regex=TRUE)[,idx]
  y <- exp(data_used$lnD75[idx])
  
  ### Median
  par <- "Median"
  xlab <- par
  p <- data.frame(gamma = rowMedians(y_rep_gam),
                  logn = rowMedians(y_rep_logn),
                  mix = rowMedians(y_rep_mix)
  ) %>%
    pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
    ggplot() +
    geom_density(aes(x = var, color = model), key_glyph = draw_key_path) +
    vline_at(median(y), size = 2) +
    labs(x = xlab, y = "density", title = sprintf("%d < R <= %d", bin_dist[k], bin_dist[k+1])) +
    theme(
      legend.position = "none"
    ) +
    guides(color = guide_legend(nrow = 3, byrow = TRUE,title=NULL)) +
    scale_color_manual(values=c("black", "red", "blue"),
                       labels=c("Gamma",
                                "Lognormal",
                                "Mixed")
    )
  ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_PPC_%s_bin%d.pdf', par,k)), p,
         width = wid, height = asp * wid)
  
  ### Mean
  par <- "Mean"
  xlab <- par
  p <- data.frame(gamma = rowMeans(y_rep_gam),
                  logn = rowMeans(y_rep_logn),
                  mix = rowMeans(y_rep_mix)
  ) %>%
    pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
    ggplot() +
    geom_density(aes(x = var, color = model)) +
    vline_at(mean(y), size = 2) +
    labs(x = xlab, y = "density", title = sprintf("%d < R <= %d", bin_dist[k], bin_dist[k+1])) +
    theme(
      legend.position = c(0.77, 0.8)
    ) +
    guides(color = guide_legend(nrow = 3, byrow = TRUE,title=NULL)) +
    scale_color_manual(values=c("black", "red", "blue"),
                       labels=c("Gamma",
                                "Lognormal",
                                "Mixed")
    )
  ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_PPC_%s_bin%d.pdf', par,k)), p,
         width = wid, height = asp * wid)
  
  qu <- 0.95
  par <- sprintf("Qu95")
  xlab <- "95% fractile"
  p <- data.frame(gamma = rowQuantiles(y_rep_gam, probs = qu),
                  logn = rowQuantiles(y_rep_logn, probs = qu),
                  mix = rowQuantiles(y_rep_mix, probs = qu)
  ) %>%
    pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
    ggplot() +
    geom_density(aes(x = var, color = model)) +
    vline_at(quantile(y, qu), size = 2) +
    labs(x = xlab, y = "density", title = sprintf("%d < R <= %d", bin_dist[k], bin_dist[k+1])) +
    theme(
      legend.position = "none"
    ) +
    guides(color = guide_legend(nrow = 3, byrow = TRUE,title=NULL)) +
    scale_color_manual(values=c("black", "red", "blue"),
                       labels=c("Gamma",
                                "Lognormal",
                                "Mixed")
    )
  ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_PPC_%s_bin%d.pdf', par,k)), p,
         width = wid, height = asp * wid)
  
  qu <- 0.05
  par <- sprintf("Qu05")
  xlab <- "5% fractile"
  p <- data.frame(gamma = rowQuantiles(y_rep_gam, probs = qu),
                  logn = rowQuantiles(y_rep_logn, probs = qu),
                  mix = rowQuantiles(y_rep_mix, probs = qu)
  ) %>%
    pivot_longer(cols = c("gamma", "logn", "mix"), names_to = "model", values_to = "var") %>%
    ggplot() +
    geom_density(aes(x = var, color = model)) +
    vline_at(quantile(y, qu), size = 2) +
    labs(x = xlab, y = "density", title = sprintf("%d < R <= %d", bin_dist[k], bin_dist[k+1])) +
    theme(
      legend.position = "none"
    ) +
    guides(color = guide_legend(nrow = 3, byrow = TRUE,title=NULL)) +
    scale_color_manual(values=c("black", "red", "blue"),
                       labels=c("Gamma",
                                "Lognormal",
                                "Mixed")
    )
  ggsave(file.path(path_base,'PLOTS_M9',sprintf('plotR_PPC_%s_bin%d.pdf', par,k)), p,
         width = wid, height = asp * wid)
}


