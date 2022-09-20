#------------------------------------------------------------------------------------
#
# Code to reproduce the pinfish growth analysis for the manuscript:
# "Cross validation for model selection: a review with examples from ecology" 
# - this script is contain function definitions for use with main script scat.R
#
# Author: Luke Yates
# Date: 29/08/2022
#
#------------------------------------------------------------------------------------

#library(tidyverse)
library(dplyr)
library(purrr)
library(ggplot2)
library(future)
library(fishmethods)
#library(FlexParamCurve)
library(RColorBrewer)
library(brms)

rm(list=ls())

#source("fish_functions_nlme.R")
options(brms.file_refit = "on_change", auto_write = TRUE)

run_date <- "2021_10_07"

sv_dir <- paste0("/media/layates/CVprimer/fits_",run_date)
if(!dir.exists(sv_dir)) dir.create(sv_dir) else message(paste0("Saving results to ",sv_dir,"/"))

#-------------------------
# Load and prepare data
#-------------------------

data(pinfish); pinfish # 670 obs
names(pinfish)[1] <- "haul"

# remove samples with undetermined sex
fishData <- pinfish %>% filter(sex != 0) %>% mutate(sex = factor(sex)) %>% as_tibble

# restrict minimum haul size
fishData %>% group_by(haul) %>% summarise(n = n()) %>% pull(n) %>% table() 
min.size <- 5 # set minimum haul size
fishData <- fishData %>% group_by(haul) %>% filter(n()>= min.size) %>% ungroup

# remove outlier
fishData <- fishData %>% filter(haul != "TBD930066") %>% mutate(haul = factor(haul)) 

fishData 


#----------------------
# Define and fit models
#----------------------

# growth functions
nlf.vB <- sl ~ (152 + Linf) * (1 - exp(-1*((1.3 + K) * (age - t0))))
nlf.G <- sl~ (152 + Linf) * exp(-exp(-(1.3 + K) * (age - t0)))
nlf.log <- sl ~ (152 + Linf)/(1 + exp(-(1.3 + K) * (age - t0)))

# specify brms formula
make_form <- function(fun, K, t, L) {
  if(K) {f.K <- K ~ sex} else {f.K <- K ~ 1}
  if(t) {f.t <- t0 ~ sex} else {f.t <- t0 ~ 1}
  if(L) {f.L <- Linf ~ sex + (1|haul)} else {f.L <- Linf ~ 1 + (1|haul)}
  bf(get(paste0("nlf.",fun)), f.K, f.t, f.L, nl = T)
}

# specify priors
make_prior <- function(K, t, L, sig_Linf = 10, sig_sex2 = 0.3){
  p <-  
    prior(normal(0, 1), nlpar = "K") +  
    prior(student_t(3, 0, 10), nlpar = "Linf", class = "sd") +
    prior(normal(0, 1), nlpar = "t0") +
    set_prior(paste0("normal(0, ",sig_Linf,")"), nlpar = "Linf")
  if(K) p <-  p + set_prior(paste0("normal(0, ",sig_sex2,")"), nlpar = "K", coef = "sex2") 
  if(t) p <-  p + set_prior(paste0("normal(0, ",sig_sex2,")"), nlpar = "t0", coef = "sex2") 
  if(L) p <-  p + set_prior(paste0("normal(0, ",sig_sex2,")"), nlpar = "Linf", coef = "sex2") 
  return(p)
}

# characterise model set
models_grid <- expand.grid(fun = c("vB","G","log"), K = c(0,1), L = c(0,1), t  = c(0,1), stringsAsFactors = F) %>% 
  mutate(name = paste0(fun,".",ifelse(K,"K",""),ifelse(L,"L",""),ifelse(t,"t",""), ifelse(K + L + t, "","0")),
         dim = K + L + t)

# fit and save models
if(F){
  seed <- 60869 #sample(1e5,1)
  m <- list()
  plan(multisession(workers = 24))
  for(i in 1:nrow(models_grid)){
    print(paste("Fitting model", i, models_grid[i,"name"])) 
    m[[models_grid[i,"name"]]] <- futureCall(brm, args = list(formula = with(models_grid[i,], make_form(fun,K,t,L)),
                                                              prior = with(models_grid[i,], make_prior(K,t,L)),
                                                              data = fishData,
                                                              future = F,
                                                              chains = 4,
                                                              seed = seed,
                                                              iter = 4000, 
                                                              file = paste0(sv_dir,"/",models_grid[i,"name"]),
                                                              file_refit = "always"),
                                             seed = T,
                                             earlySignal = T
    )
  }
}

# load all models and check convergence
m.fits <- models_grid %>% pull(name, name = name) %>% map(~ paste0(sv_dir,"/",.x,".rds") %>% readRDS)

rhat_initial <- m.fits %>% map_dbl(~ .x %>% rhat() %>% map_dbl(mean) %>% mean)
rhat_initial %>% saveRDS("files_fish/rhat_initial.rds")


# 2 models have convergence issues: refit them with narrower priors
if(F){
  m.fits$vB.KLt %>% update(cores = 4, inits = "0", seed = seed, prior = make_prior(1,1,1,sig_Linf = 5, sig_sex2 = 0.2), file = paste0(sv_dir,"/vB.KLt_update"))
  m.fits$vB.Kt %>% update(cores = 4, inits = "0", seed = seed, prior = make_prior(1,1,0,sig_Linf = 5, sig_sex2 = 0.2), file = paste0(sv_dir,"/vB.Kt_update"))
}

# refits successful: replace them in the list of fits
m.fits$vB.KLt <- readRDS(paste0(sv_dir,"/vB.KLt_update.rds"))
m.fits$vB.Kt <- readRDS(paste0(sv_dir,"/vB.Kt_update.rds"))


# compute PSIS-LOO
if(F){
  m.loo <- m.fits %>% furrr::future_map(loo)
  #saveRDS(m.loo,paste0(sv_dir,"/m.loo.rds"))
}
m.loo <- readRDS(paste0(sv_dir,"/m.loo.rds"))
m.loo$vB.0$diagnostics$pareto_k
m.loo %>% map("diagnostics") %>% map("pareto_k") %>% map_dbl(~ sum(.x>0.7))

# compute LOGO CV estimates
if(F){
  plan(multisession(workers = 45))
  m.logo.cv <- lapply(m.fits, kfold, chains = 2, future = T, group = "haul")
  #m.logo.cv %>% saveRDS(paste0(sv_dir,"/m.logo.rds"))
}

m.logo.cv<- readRDS(paste0(sv_dir,"/m.logo.rds"))


#-----------------------------------------
# fit a fixed-effect model for comparison
#-----------------------------------------

make_form_fe <- function(fun, K, t, L) {
  if(K) {f.K <- K ~ sex} else {f.K <- K ~ 1}
  if(t) {f.t <- t0 ~ sex} else {f.t <- t0 ~ 1}
  if(L) {f.L <- Linf ~ sex + haul} else {f.L <- Linf ~ 0 + haul}
  bf(get(paste0("nlf.",fun)), f.K, f.t, f.L, nl = T)
}

make_prior_fe <- function(K, t, L, sig_Linf = 10, sig_sex2 = 0.3){
    prior(normal(0, 1), nlpar = "K") +  
    prior(normal(0, 1), nlpar = "t0") +
    set_prior(paste0("student_t(3, 0, ", sig_Linf,")"), nlpar = "Linf")
}

m.vB.0.fe <- brm(formula = make_form_fe("vB",0,0,0),
                 prior = make_prior_fe(0,0,0,sig_Linf = 200),
                 data = fishData,
                 cores =4,
                 chains = 4,
                 seed = seed,
                 iter = 4000, 
                 file = paste0(sv_dir,"/m.vB.0.fe.wide"))

m.vB.0.fe
m.vB.0.fe %>% prior_summary()

make_form_pop <- function(fun, K, t, L) {
  if(K) {f.K <- K ~ sex} else {f.K <- K ~ 1}
  if(t) {f.t <- t0 ~ sex} else {f.t <- t0 ~ 1}
  if(L) {f.L <- Linf ~ sex } else {f.L <- Linf ~ 1}
  bf(get(paste0("nlf.",fun)), f.K, f.t, f.L, nl = T)
}

m.vB.0.pop <- brm(formula = make_form_pop("vB",0,0,0),
                 prior = make_prior_fe(0,0,0),
                 data = fishData,
                 cores = 4,
                 chains = 4,
                 seed = seed,
                 iter = 4000, 
                 file = paste0(sv_dir,"/m.vB.0.pop"),
                 file_refit = "always")


m.vB.0 <- readRDS(paste0(sv_dir,"/vB.0.rds"))

m.vB.0.pop <- readRDS(paste0(sv_dir,"/m.vB.0.pop.rds"))
m.vB.0.fe <- readRDS(paste0(sv_dir,"/m.vB.0.fe.wide.rds"))

m.list <- list(FE = m.vB.0.fe, POP = m.vB.0.pop, RE = m.vB.0)
m.list.loo <- m.list %>% map(loo, cores = 10)
#saveRDS(m.list.loo, paste0(sv_dir,"/m.list.loo"))
m.list.loo <- readRDS(paste0(sv_dir,"/m.list.loo"))
m.list.loo %>% loo_compare() %>% as_tibble()


# Extract fixed (FE) and random (RE) effects for Linf
L_haul_FE <- m.vB.0.fe %>% as_tibble %>% select(starts_with("b_Linf_")) %>% 
  rename_with(~ str_replace(.x, "b_Linf_haul","")) %>% 
  map_dbl(mean) %>% 
  {tibble(haul = names(.), FE = unname(.))}

L_haul_FE <- m.vB.0.fe %>% as_tibble %>% select(starts_with("b_Linf_")) %>% 
  rename_with(~ str_replace(.x, "b_Linf_haul","")) %>% 
  #mutate(across(everything(), ~ .x + (m.vB.0.fe %>% as_tibble %>% pull("b_Linf_Intercept")))) %>% 
  mutate(across(-b_Linf_Intercept, ~ .x + b_Linf_Intercept)) %>%   
  rename(BFC970009 = b_Linf_Intercept) %>% 
  map_dbl(mean) %>% 
  {tibble(haul = names(.), FE = unname(.))}


L_haul_RE <- m.vB.0 %>% as_tibble %>% select(starts_with("r_haul")) %>% 
  rename_with(~ str_replace(.x, "r_haul__Linf\\[","")) %>% 
  rename_with(~ str_replace(.x, ",Intercept\\]","")) %>% 
  mutate(across(everything(), ~ .x + (m.vB.0 %>% as_tibble %>% pull("b_Linf_Intercept")))) %>% 
  map_dbl(mean) %>% 
  {tibble(haul = names(.), RE = unname(.))}

lm <- m.vB.0 %>% as_tibble %>% pull(b_Linf_Intercept) %>% mean
m.vB.0

dark2 <- RColorBrewer::brewer.pal(8,"Dark2")
pie(rep(1, length(dark2)), col = dark2 , main="") 


# combine data and make shrinkage plot
full_join(L_haul_FE,L_haul_RE, by = "haul") %>% 
  arrange(RE) %>% 
  mutate(haul = factor(haul, levels = haul),
         FE = FE - lm) %>% 
  mutate(across(-haul, ~ .x + 152)) %>% 
  ggplot(aes(x = haul)) +
  #geom_hline(aes(yintercept = 0), lty = "dashed") +
  geom_linerange(aes(ymin = RE, ymax = FE), col = "grey60")+
  geom_point(aes(y = FE, col = "Fixed"), size = 0.5)+
  geom_point(aes(y = RE, col = "Random"), size = 0.5) +
  theme_classic() +
  scale_color_manual(values = c(Fixed = "black", Random = dark2[2])) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  labs(y = expression(L[0*", haul"]), subtitle = expression("Shrinkage of haul-level"~L[0]~"estimates"),
       col = "", x = "haul")

ggsave(paste0("plots/re_shrinkage_plot_2_",run_date,".pdf"), width = 90, height = 80, unit = "mm")


# manually compute p_loo to check summary: YES
RE.lppd <- log_lik(m.vB.0) %>% exp %>% apply(2,mean) %>% log 
RE.lppd_loo <- m.list.loo$RE$pointwise[,"elpd_loo"]
(RE.lppd - RE.lppd_loo) %>% sum
m.list.loo$RE$pointwise[,"p_loo"] %>% sum

FE.lppd <- log_lik(m.vB.0.fe) %>% exp %>% apply(2,mean) %>% log 
FE.lppd_loo <- m.list.loo$FE$pointwise[,"elpd_loo"]
(FE.lppd - FE.lppd_loo) %>% sum
m.list.loo$FE$pointwise[,"p_loo"] %>% sum

# actual parameter count for FE
m.vB.0.fe %>% as_tibble %>% names %>% length -1

(m.list.loo$RE$pointwise[,"p_loo"] - m.list.loo$FE$pointwise[,"p_loo"]) %>% 
  {c(p_loo_diff = sum(.), se_p_loo_diff = sd(.)*sqrt(length(.)))}


