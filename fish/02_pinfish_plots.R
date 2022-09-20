#------------------------------------------------------------------------------------
#
# Code to reproduce the pinfish growth plots for the manuscript:
# "Cross validation for model selection: a review with examples from ecology" 
# - this script is contain function definitions for use with main script scat.R
#
# Author: Luke Yates
# Date: 29/08/2022
#
#------------------------------------------------------------------------------------

library(ggplot2)
library(purrr)
library(dplyr)
library(future)
library(fishmethods)
library(stringr)
library(tidyr)
library(RColorBrewer)
library(brms)

rm(list=ls())

run_date <- "2021_10_07"
sv_dir <- paste0("/media/layates/CVprimer/fits_",run_date)
if(!dir.exists(sv_dir)) dir.create(sv_dir) else message(paste0("Saving results to ",sv_dir,"/"))

# load PSIS-LOO and LOGO CV estimates
m.loo <- readRDS("files_fish/m.loo.rds")
m.logo <- readRDS("files_fish/m.logo.rds")

make_plot_data <- function(metric_data, levels = names(metric_data)){
  best_model <- metric_data %>% map_dbl(mean) %>% which.max() %>% names()
  tibble(model = factor(names(metric_data), levels = levels), 
         metric = metric_data %>% map_dbl(mean),
         metric_diff = metric - max(metric),
         se = metric_data %>% map_dbl(~ .x %>% {sd(.)/sqrt(length(.))}),
         se_diff = metric_data %>% map(~ .x - metric_data[[best_model]]) %>% map_dbl(~ .x %>% {sd(.)/sqrt(length(.))}),
         se_mod = sqrt(1 -cor(metric_data)[best_model,])*se[best_model])
}


## PSIS-LOO OSE plot

models_grid <- expand.grid(fun = c("vB","G","log"), K = c(0,1), L = c(0,1), t  = c(0,1), stringsAsFactors = F) %>% 
  mutate(fe = paste0(ifelse(K,"K",""),ifelse(L,"L",""),ifelse(t,"t",""), ifelse(K + L + t, "","0")),
        name = paste0(fun,".",fe),
        dim = K + L + t) %>% 
  arrange(fun,dim)
models_grid

loo_pointwise <- m.loo %>% map("pointwise") %>% map_dfc(~ .[,"elpd_loo"]) %>% 
  relocate(all_of(models_grid$name))

se_type <- "se_mod"

loo_df <- make_plot_data(loo_pointwise) %>% 
            mutate(fe = factor(models_grid$fe, levels = models_grid$fe[1:8]),
                   fun = models_grid$fun,
                   dim = models_grid$dim,
                   se_ose = se,
                   se = .[[se_type]])

loo_ose_models <- loo_df %>% filter(metric + se >= max(metric)) %>% filter(dim == min(dim))
fun_labels = c(G = "Gompertz", log = "logistic", vB = "von Bertalanffy")

# faceted plot
loo_plot <- loo_df %>% 
  ggplot(aes(x = fe)) +
  #geom_line(aes(y = metric_diff, group = fun), col = "grey40", lty = "dashed", show.legend = F) +
  geom_point(aes(y = metric_diff, col = fun), size = 1.5, show.legend = F) +
  geom_linerange(aes(ymin = metric_diff - se, ymax = metric_diff + se, col = fun), show.legend = F) +
  theme_classic() +
  theme(strip.placement = "outside", panel.border = element_blank(), 
        strip.background = element_rect(), axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8), strip.background.x = element_rect(linetype = 0, fill = "grey90"),
        axis.title.y = element_text(size = 8)) +
  facet_wrap(~fun,nrow = 1, strip.position = "bottom", labeller = labeller(fun = fun_labels)) +
  geom_point(aes(y = metric_diff), shape = 1, size = 4, col = "black", data = loo_ose_models) +
  scale_color_brewer(type = "div", palette = "Dark2")+
  scale_y_continuous(labels=function(x)x*1000, limits = c(-0.022,0.002)) +
  labs(x = NULL, subtitle = "Pinfish PSIS-LOO estimates", y =  expression(Delta*hat(S)*" "~group("(", 10^- 3,")")))

loo_plot


# LOGO plots
logo_pointwise <- m.logo %>% map("pointwise") %>% map_dfc(~ .[,"elpd_kfold"]) %>% 
  relocate(all_of(models_grid$name))

se_type <- "se_mod"

logo_df <- make_plot_data(logo_pointwise) %>% 
  mutate(fe = factor(models_grid$fe, levels = models_grid$fe[1:8]),
         fun = models_grid$fun,
         dim = models_grid$dim,
         se_ose = se,
         se = .[[se_type]])

logo_ose_models <- logo_df %>% filter(metric + se >= max(metric)) %>% filter(dim == min(dim))
fun_labels = c(G = "Gompertz", log = "logistic", vB = "von Bertalanffy")

# faceted plot
logo_plot <- logo_df %>%
  filter(metric_diff > -0.012) %>% 
  ggplot(aes(x = fe)) +
  #geom_line(aes(y = metric_diff, group = fun), col = "grey40", lty = "dashed", show.legend = F) +
  geom_point(aes(y = metric_diff, col = fun), size = 1.5, show.legend = F) +
  geom_linerange(aes(ymin = metric_diff - se, ymax = metric_diff + se, col = fun), show.legend = F) +
  theme_classic() +
  theme(strip.placement = "outside", panel.border = element_blank(), 
        strip.background = element_rect(), axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8), strip.background.x = element_rect(linetype = 0, fill = "grey90"),
        axis.title.y = element_text(size = 8)) +
  facet_wrap(~fun,nrow = 1, strip.position = "bottom", labeller = labeller(fun = fun_labels)) +
  geom_point(aes(y = metric_diff), shape = 1, size = 4, col = "black", data = logo_ose_models) +
  scale_color_brewer(type = "div", palette = "Dark2")+
  scale_y_continuous(labels=function(x)x*1000, limits = c(-0.012,0.002)) +
  labs(x = NULL, subtitle = "Pinfish LOGO estimates", y =  expression(Delta*"Score "~group("(", 10^{- 3},")")))

logo_plot

tt <- theme(axis.text.x = element_text(size = 5))

ggpubr::ggarrange(loo_plot + tt,logo_plot + tt, labels = "AUTO", nrow = 1)

##-----------
## DATA PLOTS
##-----------

# Plot 2 - conditional model

m.vB.0 <- readRDS("files_fish/m.vB.0.rds")

m.vB.0
hauls <- m.vB.0$data$haul %>% unique
ce_data <- hauls %>% map_dfr(~ tibble(age = seq(0.5,6.3,0.04), haul = .x))

g.vB <- function(age, Linf, K, t0) (152 + Linf) * (1 - exp(-1*((1.3 + K) * (age - t0))))

ps.0 <- m.vB.0 %>% as_tibble %>%
  rename_with(~ str_replace(.x, "r_haul__Linf\\[","")) %>% 
  rename_with(~ str_replace(.x, ",Intercept\\]",""), ends_with("Intercept]")) %>% 
  rename(K = b_K_Intercept, t0 = b_t0_Intercept, Linf0 = b_Linf_Intercept, tau = sd_haul__Linf_Intercept) %>% 
  select(-lp__) %>% sample_n(500) %>% mutate(s = 1:500)

ce_plot_data <- 
  ps.0 %>% 
  pivot_longer(-any_of(c("sigma","tau","s","K","Linf0","t0")),names_to = "haul", values_to = "Linf_haul") %>% 
  mutate(Linf = Linf0 + Linf_haul) %>% 
  full_join(ce_data, by = "haul") %>% 
  group_by(s) %>% 
  mutate(sig_err = rnorm(1,0,sigma)) %>%  # generate residual error
  mutate(length_e = g.vB(age,Linf,K,t0), # expected length
         length = length_e + sig_err) %>% 
  group_by(haul) %>% 
  mutate(haul_col = mean(Linf) + 152) # haul colour index for plots

#ribbon_data <- 
##  ce_plot_data %>% 
#  group_by(age) %>% 
#  summarise(l_low = quantile(length, 0.025), l_hi = quantile(length, 0.975)) 

obs_data <- m.vB.0$data %>% left_join(ce_plot_data %>% select(haul, haul_col) %>% distinct, by = "haul")

ce_plot_cond <- 
  ce_plot_data %>% 
  group_by(haul_col,age) %>% 
  summarise(length = mean(length_e)) %>% 
  ggplot(aes(x = age)) +
  #geom_ribbon(aes(ymin = l_low, ymax = l_hi), fill = "black", alpha = 0.1, data = ribbon_data) +
  geom_line(aes(y = length, group = haul_col, col = haul_col), size = 0.25, alpha = 0.4) +
  geom_point(aes(y = sl, group = haul_col, col = haul_col), size = 0.5, alpha = 1, data = obs_data) +
  theme_classic() +
  theme(legend.position = "none", legend.key.height =  unit(5,"mm")) 

ce_plot_cond



# Plot 2 - marginal model

m.log.0 <- readRDS("files_fish/m.log.0.rds")
#hauls <- m.log.t$data$haul %>% unique
ce_data_2 <-tibble(age = seq(0.5,6.3,0.04), J = 1)

#m.log.t$formula
g.log.0 <- function(age, Linf, K, t0) (152 + Linf)/(1 + exp(-(1.3 + K) * (age - t0)))
m.log.0$formula
ps.0 <- m.log.0 %>% as_tibble %>% 
  select(K = b_K_Intercept, t0 = b_t0_Intercept, Linf0 = b_Linf_Intercept, 
         tau = sd_haul__Linf_Intercept, sigma) %>% 
  sample_n(500) %>% mutate(s = 1:500)


ce_plot_data2 <- full_join(ps.0 %>% mutate(J = 1), ce_data_2, by = "J") %>% 
  group_by(s) %>% 
  mutate(sig_err = rnorm(1,0,mean(sigma))) %>%  ## draw from residual error model
  mutate(tau_err = rnorm(1,0,mean(tau))) %>% ## draw from hierarchical model for Linf
  mutate(length = g.log.0(age,Linf0 + tau_err,K,t0) + sig_err,
         length_mean = g.log.0(age,Linf0,K,t0))


ce_plot_marg <- ce_plot_data2 %>% 
  group_by(age) %>% 
  summarise(l_low = quantile(length, 0.025),
            l_hi = quantile(length, 0.975),
            length = mean(length_mean)) %>% 
  ggplot(aes(x = age)) +
  geom_ribbon(aes(ymin = l_low, ymax = l_hi), fill = blues9[3], alpha = 0.4) +
  geom_line(aes(y = length), size = 0.7, col = blues9[9], lty = "solid") +
  geom_point(aes(y = sl), size = 0.4, alpha = 1, data = m.log.0$data) +
  theme_classic() +
  theme(legend.position = "none", legend.key.height =  unit(10,"mm"))

ce_plot_marg

yl2 <- ylim(c(45,206))

ggpubr::ggarrange(ce_plot_cond + labs(subtitle = "Pinfish data: conditional focus") + yl2, 
                  ce_plot_marg + yl2 + labs(y = "", subtitle = "Pinfish data: marginal focus"), nrow = 1, labels = "auto")

ggpubr::ggarrange(loo_plot + labs(subtitle = "Leave-one-out cross validation") + tt,
                  logo_plot + labs(subtitle = "Leave-one-group-out cross validation", y = " ") + tt,
                  ce_plot_cond + labs(subtitle = "Conditional prediction: vB|0", x = "age (years)", y = "length (mm)") + yl2, 
                  ce_plot_marg + yl2 + labs(y = "", subtitle = "Marginal prediction: log|0", x = "age (years)"), 
                  nrow = 2, ncol = 2, labels = c("(a)","(b)","(c)","(d)"), heights = c(3,3))

#ggsave(paste0("plots/pinfish_panel_2022_08_29",".pdf"), width = 180, height = 150, units = "mm", device = cairo_pdf()); dev.off()  

