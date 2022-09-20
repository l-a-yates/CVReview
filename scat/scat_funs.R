#------------------------------------------------------------------------------------
#
# Code to reproduce the scat classification analysis for the manuscript:
# "Cross validation for model selection: a review with examples from ecology" 
# - this script is contain function definitions for use with main script scat.R
#
# Author: Luke Yates
# Date: 29/08/2022
#
#------------------------------------------------------------------------------------


##---------------
## prepare data
##--------------- 

# select and transform initial variables; remove NAs
# merge response into two categories rather than three
prep_data <- function(data){
  data %>% as_tibble %>% 
    select(y = Species, Number, Length, Diameter, Taper, TI, Mass, CN, Location, ropey, segmented) %>% 
    rename_with(str_to_lower) %>% 
    mutate(across(c(mass,cn),log)) %>% # following original data publication
    mutate(number = as.double(number)) %>% 
    drop_na() %>% 
    mutate(across(where(is_double), ~ .x %>% scale %>% as.numeric)) %>% # centre and scale predictors
    mutate(y = fct_collapse(y, canid = c("coyote","gray_fox"), felid = c("bobcat")) %>% 
             fct_relevel("canid","felid"),
           location = fct_collapse(location, mid = c("middle"), edge = c("edge","off_edge")) %>% 
             factor(labels = c(0,1))) %>% 
    mutate(across(c(ropey,segmented), factor))
}



##---------------------
## Plotting functions
##---------------------


# computes summary statistics for model score estimates
make_plot_data <- function(metric_data, levels){
  best_model <- metric_data %>% map_dbl(mean) %>% which.max() %>% names()
  tibble(model = factor(names(metric_data), levels = levels), 
         metric = metric_data %>% map_dbl(mean),
         se = metric_data %>% map_dbl(sd),
         se_diff = metric_data %>% map(~ .x - metric_data[[best_model]]) %>% map_dbl(sd),
         se_mod = sqrt(1 -cor(metric_data)[best_model,])*se[best_model])
}

# plots estimates of the model scores and the selected uncertainty measure
plot_model_comparisons <- function(plot_data, se_type = c("se_mod", "se_diff", "se")){
  plot_data <- plot_data %>% mutate(se = plot_data[[se_type[1]]])
  plot_data %>% 
    ggplot(aes(model)) +
    geom_point(aes(y = metric)) +
    geom_linerange(aes(ymin = metric - se, ymax = metric + se)) +
    theme_bw() 
}


##-------------------
## Part 1 functions
##-------------------


# perform variable step selection using MLE and a supplied metric
# @vars: a vector variable names
# @cv_data: a vfold_cv object from the rsample package
# @metric: the name of a function which computes a metric from a 2x2 confusion-matrix 
step_glm <- function(vars, cv_data, metric){
  fit_fun <- fit_confusion_glm
  if(metric == "log_density") fit_fun <- fit_ll_glm else metric <- get(metric)
  sel <- c()
  metric_step <- list()
  make_formula <- function(vars) paste("y ~ 1 + ",paste(vars,collapse=" + "))
  
  for(step in 1:length(vars)){
    message(paste("Step:",step))
    newfits = pbmclapply(setdiff(vars,sel), function(var){
      fit_fun(cv_data, make_formula(c(sel,var)), metric)
    }, mc.cores = length(setdiff(vars,sel)))
    best = newfits %>% map("metric") %>% map_dbl(mean) %>% which.max()
    sel[step] <- setdiff(vars,sel)[best]
    message(paste("Selected:",sel[step]))
    metric_step[[step]] <- newfits[[best]]
  }
  list(sel = sel, metric_step = metric_step)
}


# Computes (L) estimates of the selected confusion-matrix metric
# @cv_data: a repeated vfold_cv object from the rsample package with L repeats
# @formula: if supplied, the same formula is applied across all CV samples, otherwise the formulas 
#    must be defined row-by-row in a column called `form_glm` in the vfold_cv object
# @metric: the name of a function which computes a metric from a 2x2 confusion-matrix 
fit_confusion_glm <- function(cv_data, formula = NULL, metric){
  get_confusion_matrix <- function(form, split, thresh){
    glm(form, data = analysis(split), family = binomial()) %>% 
      {tibble(pred = predict(.,newdata = assessment(split), type = "response") >= thresh, 
              y = assessment(split)$y != levels(assessment(split)$y)[1])} %>% 
      mutate(across(everything(), ~ .x %>% as.numeric() %>% factor(levels = c(1,0)))) %>% 
      table %>% {c(tp = .[1,1], fp = .[1,2], fn = .[2,1], tn = .[2,2])}
  }
  if(!is.null(formula)) form <- list(formula) else form <- cv_data$form_glm
  with(cv_data, mapply(get_confusion_matrix, form, splits, thresh, SIMPLIFY = T)) %>% t %>% 
    as_tibble %>% 
    bind_cols(rep = cv_data$id, .) %>% 
    group_by(rep) %>% 
    summarise(metric = matrix(c(sum(tp),sum(fn),sum(fp),sum(tn)),2,2) %>% metric) %>% 
    mutate(rep = 1:n())
}

# Computes the L estimates of the log probability of the logistic model
# @cv_data: a repeated vfold_cv object from the rsample package with L repeats
# @formula: formula to be applied across all CV samples
fit_ll_glm <- function(cv_data, formula, metric = NULL){
  calc_pred_ll <- function(form, split){
    fit <- glm(form, data = analysis(split), family = binomial())
    tibble(y_pred = predict(fit, newdata = assessment(split),type = "response"),
           y_data = as.numeric(assessment(split)$y != levels(assessment(split)$y)[1]),
           log_p = log(abs(y_data - 1 + y_pred))) %>% pull(log_p) %>% sum
  }
  with(cv_data, mapply(calc_pred_ll, list(formula), splits, SIMPLIFY = T)) %>% 
    bind_cols(rep = cv_data$id, ll = .) %>% 
    group_by(rep) %>% 
    summarise(metric = sum(ll)) %>% 
    mutate(rep = 1:n())
}

# computes TSS from a confusion matrix `cm`
tss <- function (cm) {
  sens <- cm[1,1]/(cm[2,1] + cm[1,1])
  spec <- cm[2,2]/(cm[1,2] + cm[2,2])
  sens + spec - 1
}

# computes MCC from a confusion matrix `cm`
mcc <- function(cm){
  (cm[1,1]*cm[2,2] - cm[1,2]*cm[2,1])/
    (sqrt((cm[1,1]+cm[1,2])*(cm[1,1]+cm[2,1])*(cm[2,2]+cm[1,2])*(cm[2,2]+cm[2,1])))
}


##-------------------
## Part 2 functions
##-------------------



# compute out-sample-sample (i.e., CV) confusion-matrix entries for all lambda values for a given split
# @split an rsample object
# @rep name of outer repetition
# @fold index of inner fold
# @alpha as per elastic net parameter
# @lambda as per elastic net parameter
get_cm <- function(split,rep,fold, alpha, lambda, thresh){
  glmnet(x = analysis(split) %>% select(-y) %>% as.data.frame() %>% makeX(),
         y = analysis(split) %>% pull(y),
         family = "binomial",
         alpha = alpha,
         lambda = lambda) %>% 
    predict(newx = assessment(split) %>% select(-y) %>% as.data.frame() %>% makeX(),
            type = "response") %>%
    as_tibble %>%
    rename_with(~str_remove(.x,stringr::fixed("s"))) %>% 
    imap_dfr(~ tibble(y_pred = as.numeric(.x > thresh), 
                 y = assessment(split) %>% pull(y) %>% as.numeric() %>% {.-1},
                 tp = as.numeric(y==1 & y_pred==1),
                 fp = as.numeric(y==0 & y_pred==1),
                 fn = as.numeric(y==1 & y_pred==0),
                 tn = as.numeric(y==0 & y_pred==0),
                 rep = rep, fold = fold, lambda = lambda[as.numeric(.y)+1]))
}
  

##-------------------
## Part 3 functions
##-------------------


# determine probability threshold value that maximises the metric
metric_thresh <- function(tbl, metric){
  probs<- tbl$probs
  y <- tbl$y
  probs_to_metric <- function(thresh){
    pred  = as.numeric(probs >= thresh) %>% as.numeric %>% factor(levels = c(1,0))
    y = y %>% factor(levels = c(1,0))
    metric(table(pred,y))
  }
  threshold = (1:100)/100
  metric =  threshold %>% map_dbl(probs_to_metric)
  tibble(threshold,metric) %>% arrange(-metric) %>% slice(1) %>% as_vector
}

# Tunes tree depth and probability threshold for a random forest model
# @cv_data: nested cv object from rsamples
# @mtry_vec: vector of candidate mtry values
# @ntree: number of trees (fixed)
tune_rf <- function(cv_data, mtry_vec, ntree, metric){
  pbmclapply(cv_data$inner_resamples, function(sample){
    lapply(mtry_vec, function(m){
      sample[["splits"]] %>% map(
        ~ randomForest(y ~ ., mtry = m, data = analysis(.x), ntree = ntree) %>% 
          predict(newdata = assessment(.x), type = "prob") %>% 
          {tibble(probs = .[,"felid"], y = as.numeric(assessment(.x)$y =="felid"))}
      ) %>% bind_rows %>% metric_thresh(metric)
    }) %>% bind_rows %>% mutate(mtry = mtry_vec) %>% arrange(-metric) %>% slice(1)
  }, mc.cores = MAX_CORES) %>% bind_rows
}


# tune probability threshold and step-select variables for each (outer) training set
tune_glm_step <- function(cv_data, vars, metric){
  make_formula <- function(vars) paste("y ~ 1 + ",paste(vars,collapse=" + "))
    pbmclapply(cv_data_2$inner_resamples, function(sample){
    sel <- thresh <- c()
    metric_step = list()
    
    for(step in 1:length(vars)){
      message(paste("Step:",step))
      newfits = lapply(setdiff(vars,sel), function(var){
        sample[["splits"]] %>% map(
          ~ glm(make_formula(c(sel,var)), data = analysis(.x), family = binomial()) %>% 
            {tibble(probs = predict(.,newdata = assessment(.x), type = "response"), 
                    y = as.numeric(assessment(.x)$y != levels(assessment(.x)$y)[1]))}
        ) %>% bind_rows %>% metric_thresh(metric)
      })
      best = newfits %>% map("metric") %>% map_dbl(mean) %>% which.max()
      thresh[step] = newfits[[best]][["threshold"]]
      sel[step] = setdiff(vars,sel)[best]
      message(paste("Selected:",sel[step]))
      metric_step[[step]] = newfits[[best]]
    }
    best_model = metric_step %>% map_dbl("metric") %>% which.max
    list(form = make_formula(sel[1:best_model]), 
         threshold = thresh[best_model])
  }, mc.cores = MAX_CORES) %>% bind_rows
}

# tune probability threshold and select variables for each (outer) training set
tune_glm_all <- function(cv_data, models, metric){
  pbmclapply(cv_data$inner_resamples, function(sample){
    sapply(models, function(form){
      sample$splits %>% 
        map(~ glm(form, data = analysis(.x), family = binomial()) %>% 
              {tibble(probs = predict(.,newdata = assessment(.x), type = "response"), 
                      y = as.numeric(assessment(.x)$y != levels(assessment(.x)$y)[1]))}
        ) %>% bind_rows %>% metric_thresh(metric)
    }) %>% t %>% as_tibble %>% mutate(form = unname(models), model = names(models)) %>% arrange(-metric) %>% slice(1)
  }, mc.cores = MAX_CORES) %>% bind_rows()
}


# Fit the random forest model and return metric estimates using tuned hyperparameters for each sample
fit_rf <- function(cv_data, ntree = 500, type = c("all","best"), metric){
  make_formula = function(vars) paste("y ~ 1 + ",paste(vars,collapse=" + ")) %>% as.formula()
  type = type[1]
  get_confusion_matrix <- function(split, thresh, mtry, ntree){
    success = levels(assessment(split)$y)[2]
    if(type == "best"){ # select the most important vars up to mtry terms
      fit_init = randomForest(y ~ ., mtry = mtry, data = analysis(split), ntree = ntree)
      form = fit_init %>% importance() %>% as.data.frame %>% 
        arrange(-MeanDecreaseGini) %>% slice(1:mtry) %>% row.names() %>% make_formula
    } else form = formula(y~.)
    randomForest(form, mtry = mtry, data = analysis(split), ntree = ntree) %>% 
      {tibble(pred = predict(., newdata = assessment(split), type = "prob")[,success] >= thresh,
              y = assessment(split)$y == success)} %>% 
      mutate(across(everything(), ~ .x %>% as.numeric() %>% factor(levels = c(1,0)))) %>% 
      table %>% {c(tp = .[1,1], fp = .[1,2], fn = .[2,1], tn = .[2,2])}
  }
  with(cv_data, pbmcmapply(get_confusion_matrix, splits, thresh_rf, mtry_rf, list(ntree), mc.cores = MAX_CORES)) %>% 
    simplify %>% t %>% 
    as_tibble %>% 
    mutate(rep = cv_data$id) %>% 
    group_by(rep) %>% 
    summarise(metric = matrix(c(sum(tp),sum(fn),sum(fp),sum(tn)),2,2) %>% metric) %>% 
    mutate(rep = 1:n())
}

