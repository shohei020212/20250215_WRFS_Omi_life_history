
# Modified by Noda on 2024.11.28

# 1. setup ---------------------------------------------------------------------

# package
{
  library("dplyr")
  require(ggplot2)
  library(lme4)
  library(RColorBrewer)
  library(stringr)
  library(flextable)
  library(partitions)
  library(MuMIn)
  library(rcompanion)
  options(na.action = "na.fail")
  library(ggeffects)
  library(multcomp)
  }

# reset environment
rm(list=ls(all.names=T))

# 2. prepare data --------------------------------------------------------------

# read raw data
d0 <- read.csv("./output/rawdata_filtered.csv")
idlist <- read.csv("./output/id_list.csv")
gc <- read.csv("./output/k4_data.csv")


## id list
# filtering
d1 <- idlist %>% 
  left_join(gc, by="id") %>% 
  mutate(early_mature = case_when(age_at_mature == "0" ~ "1",
                                  age_at_mature == "1" ~ "1",
                                  age_at_mature == "2" ~ "0",
                                  age_at_mature == "3" ~ "0",
                                  age_at_mature == ">1" ~ "0",
                                  age_at_mature == ">2" ~ "0",
                                  age_at_mature == ">3" ~ "0",
                                  TRUE ~ "nd"),
         growth_cluster = case_when(ocluster_rename == 1 ~ '1',
                                    ocluster_rename == 2 ~ '2',
                                    ocluster_rename == 3 ~ '3',
                                    ocluster_rename == 4 ~ '4',
                                    TRUE ~ "nd")) %>% 
  filter(treatment == 'early' | treatment == 'late' | treatment == 'control') %>% 
  mutate(treatment = factor(treatment, levels = c("early", "late", "control")),
         growth_cluster = as.factor(growth_cluster))

## mark-recapture data
mr_data <- d0 %>% 
  mutate(early_mature = case_when(age_at_mature == "0" ~ "1",
                                  age_at_mature == "1" ~ "1",
                                  age_at_mature == "2" ~ "0",
                                  age_at_mature == "3" ~ "0",
                                  age_at_mature == ">1" ~ "0",
                                  age_at_mature == ">2" ~ "0",
                                  age_at_mature == ">3" ~ "0",
                                  TRUE ~ "nd")) %>% 
  filter(between(cohort,2016,2021),
         !str_detect(fl, "\\?"),
         !str_detect(fl, "nd"),
         treatment == 'early' | treatment == 'late' | treatment == 'control',
         early_mature != "nd",
         age_at_mature != "0") %>% 
  mutate(fl = as.numeric(fl),
         treatment = factor(treatment, levels = c("control", "early", "late")),
         early_mature = as.numeric(early_mature))


# 3. Statistical_analysis ------------------------------------------------------
## 3-1. GLM: early_mature ~ growth_cluster * cohort-----------------------------

# input data
data <- d1 %>% 
  filter(early_mature != "nd",
         growth_cluster != "nd",
         cohort != "2017",
         !(cohort == 2016 & treatment == "early"),
         !(cohort == 2016 & treatment == "late"),) %>% 
  mutate(early_mature = as.numeric(early_mature),
         cohort = as.factor(cohort))

### 3-1-1. model selection -----------------------------------------------------

# AIC
models <- glm(early_mature ~ growth_cluster * cohort, 
              data=data, 
              family=binomial(link="logit"))

model.list <- dredge(models,rank="AIC")
model.list

# best model
best.model <- get.models(model.list, subset = 1)[[1]]
summary(best.model)

summary(models)

# ANOVA
# analysis of deviance

# null model
model.null <- glm(early_mature ~ 1, 
               data=data, family=binomial(link="logit"))

# full model
model.a <- glm(early_mature ~ growth_cluster * cohort, 
               data=data, family=binomial(link="logit"))

anova(model.a, test = "Chisq")

# remove co-efficient
model.b <- glm(early_mature ~ growth_cluster + cohort, 
               data=data, family=binomial(link="logit"))

anova(model.b, test = "Chisq")

# compare models
anova(model.a, model.b, test = "Chisq")
anova(model.null, model.b, test = "Chisq")

# best model
result <- glm(early_mature ~ growth_cluster + cohort, 
               data=data, family=binomial(link="logit"))

# odds ratio
odds <- as.data.frame(exp(result$coefficients))
odds$odds <- formatC(odds$`exp(result$coefficients)`, format = "fg")
odds

# summary to data frame
summary_table <- as.data.frame(summary(result)[[12]])
summary_table <- summary_table %>%
  mutate(Variable = rownames(summary_table)) %>% 
  dplyr::select(., "Variable", "Estimate", SE = "Std. Error", z = "z value",  p = "Pr(>|z|)") %>% 
  mutate(Estimate = round(Estimate, 3),
         SE = round(SE, 3),
         z = round(z, 3),
         p = round(p, 3),
         Signif = case_when(p < 0.001 ~ "**",
                            p <= 0.05 ~ "*",
                            TRUE ~ ""),
         p = as.character(p),
         p = case_when(Signif == "**" ~ "<0.001",
                       Signif != "**" ~ p)) %>% 
  rename("z value" = z,
         "p value" = p,
         "Signif. code" = Signif)

# make table
table <-
  summary_table %>%
  flextable::regulartable() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "sans", part = "all") %>%
  flextable::autofit() %>% 
  flextable::bold(part = "header") %>% 
  flextable::align(align = "right", part = "header") %>% 
  flextable::align(align = "right", part = "body") %>% 
  flextable::bg(bg = 'white', part = "all")
table

# output
save_as_image(table, path = "output/summary_glm_early_mature-growth_cluster+cohort.png")
write.csv(summary_table, "output/summary_glm_early_mature-growth_cluster+cohort.csv")

# Multiple comparison test 
model <- glm(early_mature ~ growth_cluster + cohort, 
             data=data, family=binomial(link="logit"))

multicomp <- glht(model, mcp(growth_cluster="Tukey"))
summary(multicomp)

model <- glm(early_mature ~ growth_cluster + cohort, 
             data=data, family=binomial(link="logit"))

multicomp <- glht(model, mcp(cohort="Tukey"))
summary(multicomp)

### 3-1-2. plot ----------------------------------------------------------------
## model
model <- glm(early_mature ~ cohort + growth_cluster, 
             data=data, 
             family=binomial(link="logit"))

# setting
cols = brewer.pal(6, "BrBG")

# prediction
preds <- ggpredict(model, terms = c("growth_cluster","cohort"), ci_level = 0.95)

# plot prediction
preds <- tibble(
  growth_cluster = as.numeric(preds$x),
  cohort = as.character(preds$group),
  predicted = preds$predicted,
  conf.low = preds$conf.low,
  conf.high = preds$conf.high
)

### Warning: SD of gorwth cluster4 were removed
preds <- preds %>% 
  mutate(conf.low = case_when(growth_cluster == 4 ~ 0,
                              TRUE ~ conf.low),
         conf.high = case_when(growth_cluster == 4 ~ 0,
                               TRUE ~ conf.high))

# convert as numeric
data$growth_cluster2 <- as.numeric(data$growth_cluster)

g <- 
  ggplot(NULL) +
  geom_jitter(data = data, aes(x = growth_cluster2, y = early_mature, colour = cohort), 
              size = 2, height=0, width = 0.3, alpha=0.3) +
  # prediction
  geom_errorbar(data = preds, 
                aes(x = growth_cluster, y = predicted, 
                    ymin = conf.low, ymax = conf.high,
                    colour = cohort),
                width = 0, linewidth = 1.5, position = position_dodge(0.5)) +
  geom_point(data = preds, aes(x = growth_cluster, y = predicted, colour = cohort),
             size = 5, position = position_dodge(0.5)) +
  
  scale_color_manual('cohort: ', values = cols) +
  scale_x_continuous(limits = c(0.5, 4.5), breaks=seq(0,5,1)) +
  scale_y_continuous(breaks=seq(0,1.5,0.5)) +
  theme_bw() +
  theme(axis.title=element_text(size=30, family="sans", face="bold"),
        axis.text=element_text(size=20, family="sans"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = 'top') +
  labs(x = "Growth cluster", y = "Prob. of maturation")
g

# save
ggsave(file = "./output/figure_logistic_early_mature-growth_cluster+cohort.png", 
       plot = g, dpi = 200, width = 9, height = 6)


### 3-1-3. grouping for model selection ----------------------------------------
#### cohort --------------------------------------------------------------------
# make AIC data frame
models_data <- data.frame(
  "lank" = NULL,
  "group" = NULL,
  "model" = NULL,
  "AIC" = NULL)

# input data
indata <- data %>% 
  mutate(subgr = case_when(cohort == "2016" ~ 1,
                           cohort == "2018" ~ 2,
                           cohort == "2019" ~ 3,
                           cohort == "2020" ~ 4,
                           cohort == "2021" ~ 5,
  ),
  subgr = as.factor(subgr),
  y = as.numeric(early_mature),
  x = as.factor(growth_cluster))

# check combination
listParts(length(unique(indata$subgr)))
level <- levels(indata$subgr)

# read model for apartition file
source("code/model_for_apartition_binominal_glm_subgr+x.R")

subsets <- listParts(length(level))
partseq <- 2:length(subsets)
result <- purrr::map(partseq, glm4apartition)

# Best model
bestresult <- result[[which.min(purrr::map(result, function(ap) ap$model$aic))]]
best_subset <- subsets[bestresult$partition]
bestmodel <- bestresult$model
summary(bestmodel)

sorted_indices <- order(purrr::map_dbl(result, ~.x$model$aic))

# Second model
second_result <- result[[sorted_indices[2]]]
second_subset <- subsets[second_result$partition]
second_model <- second_result$model

# AIC
null <- glm(y ~ subgr + x,
            family = binomial(link="logit"),
            data = indata)

est = c(AIC(null, bestmodel, second_model)$AIC)
models_data_sub <- data.frame(
  "lank" = c("NULL",
             "1",
             "2"),
  "group" = c("NULL", as.character(best_subset), as.character(second_subset)),
  "model" = "cohort+growth_cluster",
  "AIC" = c(est)
)
# merge aic data
models_data <- rbind(models_data, models_data_sub)
models_data

# summary to data frame
summary_table <- as.data.frame(summary(bestmodel)[[12]])
summary_table <- summary_table %>%
  mutate(Variable = rownames(summary_table)) %>% 
  dplyr::select(., "Variable", "Estimate", SE = "Std. Error", z = "z value",  p = "Pr(>|z|)") %>% 
  mutate(Estimate = round(Estimate, 3),
         SE = round(SE, 3),
         z = round(z, 3),
         p = round(p, 3),
         Signif = case_when(p < 0.001 ~ "**",
                            p <= 0.05 ~ "*",
                            TRUE ~ ""),
         p = as.character(p),
         p = case_when(Signif == "**" ~ "<0.001",
                       Signif != "**" ~ p)) %>% 
  rename("z value" = z,
         "p value" = p,
         "Signif. code" = Signif) %>% 
  mutate(Variable = gsub(Variable, pattern = "subgr2", replacement = "cohort_subgr2"),
         Variable = gsub(Variable, pattern = "x", replacement = "growth_culster"),)

# make table
table <-
  summary_table %>%
  flextable::regulartable() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "sans", part = "all") %>%
  flextable::autofit() %>% 
  flextable::bold(part = "header") %>% 
  flextable::align(align = "right", part = "header") %>% 
  flextable::align(align = "right", part = "body") %>% 
  flextable::bg(bg = 'white', part = "all")
table

# output
save_as_image(table, path = "output/summary_glm_early_mature-growth_cluster+cohort_best_grouping_cohort.png")

#### growth cluster ------------------------------------------------------------

# input data
indata <- data %>% 
  mutate(cohort2 = case_when(cohort == "2016" ~ 1,
                             cohort == "2018" ~ 1,
                             cohort == "2019" ~ 1,
                             cohort == "2020" ~ 2,
                             cohort == "2021" ~ 2,
  ),
  subgr = as.factor(growth_cluster),
  y = as.numeric(early_mature),
  x = as.factor(cohort2))

# check combination
listParts(length(unique(indata$subgr)))
level <- levels(indata$subgr)

# read model for apartition file
source("code/model_for_apartition_binominal_glm_subgr+x.R")

subsets <- listParts(length(level))
partseq <- 2:length(subsets)
result <- purrr::map(partseq, glm4apartition)

# Best model
bestresult <- result[[which.min(purrr::map(result, function(ap) ap$model$aic))]]
best_subset <- subsets[bestresult$partition]
best_subset
bestmodel <- bestresult$model
summary(bestmodel)

sorted_indices <- order(purrr::map_dbl(result, ~.x$model$aic))

# Second model
second_result <- result[[sorted_indices[2]]]
second_subset <- subsets[second_result$partition]
second_model <- second_result$model

# AIC
null <- glm(y ~ subgr + x,
            family = binomial(link="logit"),
            data = indata)
est = c(AIC(null, bestmodel, second_model)$AIC)
models_data_sub <- data.frame(
  "lank" = c("NULL",
             "1",
             "2"),
  "group" = c("NULL", as.character(best_subset), as.character(second_subset)),
  "model" = "cohort2+growth_cluster",
  "AIC" = c(est)
)
# merge aic data
models_data <- rbind(models_data, models_data_sub)
models_data

# summary to data frame
summary_table <- as.data.frame(summary(bestmodel)[[12]])
summary_table <- summary_table %>%
  mutate(Variable = rownames(summary_table)) %>% 
  dplyr::select(., "Variable", "Estimate", SE = "Std. Error", z = "z value",  p = "Pr(>|z|)") %>% 
  mutate(Estimate = round(Estimate, 3),
         SE = round(SE, 3),
         z = round(z, 3),
         p = round(p, 3),
         Signif = case_when(p < 0.001 ~ "**",
                            p <= 0.05 ~ "*",
                            TRUE ~ ""),
         p = as.character(p),
         p = case_when(Signif == "**" ~ "<0.001",
                       Signif != "**" ~ p)) %>% 
  rename("z value" = z,
         "p value" = p,
         "Signif. code" = Signif) %>% 
  mutate(Variable = gsub(Variable, pattern = "subgr", replacement = "growth_culster_subgr"),
         Variable = gsub(Variable, pattern = "x", replacement = "cohort_subgr"),)

# make table
table <-
  summary_table %>%
  flextable::regulartable() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "sans", part = "all") %>%
  flextable::autofit() %>% 
  flextable::bold(part = "header") %>% 
  flextable::align(align = "right", part = "header") %>% 
  flextable::align(align = "right", part = "body") %>% 
  flextable::bg(bg = 'white', part = "all")
table

# output
save_as_image(table, path = "output/summary_glm_early_mature-growth_cluster+cohort_best_grouping_growth_cluster.png")

# output model selecion data
write.csv(models_data, "output/model_select_glm_early_mature-growth_cluster+cohort.csv")

## 3-2. GLM: early_mature ~ growth_cluster * treatment--------------------------

# input data
data <- d1 %>% 
  filter(early_mature != "nd",
         growth_cluster != "nd",
         cohort == "2016") %>% 
  mutate(early_mature = as.numeric(early_mature),
         growth_cluster = as.factor(growth_cluster),
         cohort = as.factor(cohort))

### 3-2-1. model selection -----------------------------------------------------
# AIC
models <- glm(early_mature ~ growth_cluster * treatment, 
              data=data, 
              family=binomial(link="logit"))

model.list <- dredge(models,rank="AIC")
model.list

# best model
best.model <- get.models(model.list, subset = 1)[[1]]
summary(best.model)

## model
result <- glm(early_mature ~ growth_cluster + treatment, 
              data=data, 
              family=binomial(link="logit"))
summary(result)

# odds ratio
odds <- as.data.frame(exp(result$coefficients))
odds$odds <- formatC(odds$`exp(result$coefficients)`, format = "fg")
odds

# summary to data frame
summary_table <- as.data.frame(summary(result)[[12]])
summary_table <- summary_table %>%
  mutate(Variable = rownames(summary_table)) %>% 
  dplyr::select(., "Variable", "Estimate", SE = "Std. Error", z = "z value",  p = "Pr(>|z|)") %>% 
  mutate(Estimate = round(Estimate, 3),
         SE = round(SE, 3),
         z = round(z, 3),
         p = round(p, 3),
         Signif = case_when(p < 0.001 ~ "**",
                            p <= 0.05 ~ "*",
                            TRUE ~ ""),
         p = as.character(p),
         p = case_when(Signif == "**" ~ "<0.001",
                       Signif != "**" ~ p)) %>% 
  rename("z value" = z,
         "p value" = p,
         "Signif. code" = Signif)

# make table
table <-
  summary_table %>%
  flextable::regulartable() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "sans", part = "all") %>%
  flextable::autofit() %>% 
  flextable::bold(part = "header") %>% 
  flextable::align(align = "right", part = "header") %>% 
  flextable::align(align = "right", part = "body") %>% 
  flextable::bg(bg = 'white', part = "all")
table

# output
save_as_image(table, path = "output/summary_glm_early_mature-growth_cluster+treatment.png")
write.csv(summary_table, "output/summary_glm_early_mature-growth_cluster+treatment.csv")


### 3-2-2. plot ----------------------------------------------------------------
## model
model <- glm(early_mature ~ growth_cluster + treatment, 
             data=data, 
             family=binomial(link="logit"))

# setting
cols = c("green4","peru","dimgray")

# prediction
preds <- ggpredict(model, terms = c("growth_cluster","treatment"), ci_level = 0.95)

# plot prediction
preds <- tibble(
  growth_cluster = as.numeric(preds$x),
  treatment = as.character(preds$group),
  predicted = preds$predicted,
  conf.low = preds$conf.low,
  conf.high = preds$conf.high
)

# convert as numeric
data$growth_cluster2 <- as.numeric(data$growth_cluster)

g <- 
  ggplot(NULL) +
  geom_jitter(data = data, aes(x = growth_cluster2, y = early_mature, colour = treatment), 
              size = 2, height=0, width = 0.3, alpha=0.3) +
  # prediction
  geom_errorbar(data = preds, 
                aes(x = growth_cluster, y = predicted, 
                    ymin = conf.low, ymax = conf.high,
                    colour = treatment),
                width = 0, linewidth = 1.5, position = position_dodge(0.5)) +
  geom_point(data = preds, aes(x = growth_cluster, y = predicted, colour = treatment),
             size = 5, position = position_dodge(0.5)) +
  
  scale_color_manual('treatment: ', values = cols) +
  scale_x_continuous(limits = c(0.5, 4.5), breaks=seq(0,5,1)) +
  scale_y_continuous(breaks=seq(0,1.5,0.5)) +
  theme_bw() +
  theme(axis.title=element_text(size=30, family="sans", face="bold"),
        axis.text=element_text(size=20, family="sans"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = 'top') +
  labs(x = "Growth cluster", y = "Prob. of maturation")
g

# save
ggsave(file = "./output/figure_logistic_early_mature-growth_cluster+treatment.png", 
       plot = g, dpi = 200, width = 9, height = 6)



## 3-3. GLM: early_mature ~ fl * treatment -------------------------------------

# filtering
d2 <- mr_data %>%
  filter(cohort == "2016",
         between(month, 9,11),
         year == cohort + 1) %>% 
  mutate(cohort = as.factor(cohort),
         early_mature = as.numeric(early_mature))

# check data
count(d2, treatment, early_mature)

### 3-3-1. model selection -----------------------------------------------------
## logistic regression
models <- glm(early_mature ~ fl * treatment, 
              data=d2, 
              family=binomial(link="logit"))

model.list <- dredge(models,rank="AIC")
model.list

# best model
best.model <- get.models(model.list, subset = 1)[[1]]
summary(best.model)

# model
result <- glm(early_mature ~ fl + treatment, 
              data=d2, 
              family=binomial(link="logit"))
summary(result)

# odds ratio
odds <- as.data.frame(exp(result$coefficients))
odds$odds <- formatC(odds$`exp(result$coefficients)`, format = "fg")
odds

# summary to data frame
summary_table <- as.data.frame(summary(result)[[12]])
summary_table <- summary_table %>%
  mutate(Variable = rownames(summary_table),
         Variable = gsub(Variable, pattern = "treatment", replacement = "treatment_")) %>% 
  dplyr::select(., "Variable", "Estimate", SE = "Std. Error", z = "z value",  p = "Pr(>|z|)") %>% 
  mutate(Estimate = round(Estimate, 3),
         SE = round(SE, 3),
         z = round(z, 3),
         p = round(p, 3),
         Signif = case_when(p < 0.001 ~ "**",
                            p <= 0.05 ~ "*",
                            TRUE ~ ""),
         p = as.character(p),
         p = case_when(Signif == "**" ~ "<0.001",
                       Signif != "**" ~ p)) %>% 
  rename("z value" = z,
         "p value" = p,
         "Signif. code" = Signif)

# make table
table <-
  summary_table %>%
  flextable::regulartable() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "sans", part = "all") %>%
  flextable::autofit() %>% 
  flextable::bold(part = "header") %>% 
  flextable::align(align = "right", part = "header") %>% 
  flextable::align(align = "right", part = "body") %>% 
  flextable::bg(bg = 'white', part = "all")
table

# output
save_as_image(table, path = "output/summary_glm_early_mature-fl+treatment.png")
write.csv(summary_table, "output/summary_glm_early_mature-fl+treatment.csv")


### 3-3-2. plot ----------------------------------------------------------------

## model
model <- glm(early_mature ~ fl + treatment, data=d2, family=binomial(link="logit"))

# setting
cols = c("dimgray","green4","peru")

# prediction
preds <- ggpredict(model, terms = c("fl[all]","treatment"), ci_level = 0.95)

# plot prediction

preds <- tibble(
  fl = as.numeric(preds$x),
  treatment = as.character(preds$group),
  predicted = preds$predicted,
  conf.low = preds$conf.low,
  conf.high = preds$conf.high
)

# 50%  probability of maturing
len50.c <- -(model$coefficients["(Intercept)"])/model$coefficients["fl"]
len50.e <- -(model$coefficients["(Intercept)"] + model$coefficients["treatmentearly"])/model$coefficients["fl"]
len50.l <- -(model$coefficients["(Intercept)"] + model$coefficients["treatmentlate"])/model$coefficients["fl"]

# plot
g <- 
  ggplot(NULL) +
  geom_jitter(data = d2, aes(x = fl, y = early_mature, colour = treatment), 
              size = 2, height=0, width = 0.2, alpha=0.3) +
  # prediction
  #geom_ribbon(data = preds, aes(x = fl, xend = fl, y = conf.low, yend = conf.high, colour = treatment), alpha=0.1) +
  geom_line(data = preds, aes(x = fl, y = predicted, colour = treatment), linewidth = 2) +
  
  # lines of 50% probability of maturing
  geom_vline(xintercept=len50.c, linewidth=0.8, linetype="solid", color = cols[1]) +
  geom_vline(xintercept=len50.e, linewidth=0.8, linetype="solid", color = cols[2]) +
  geom_vline(xintercept=len50.l, linewidth=0.8, linetype="solid", color = cols[3]) +
  
  scale_color_manual('treatment: ', values = cols) +
  scale_x_continuous(limits = c(min(d2$fl)-10, max(d2$fl)+10)) +
  scale_y_continuous(breaks=seq(0,1.5,0.5)) +
  theme_bw() +
  theme(axis.title=element_text(size=30, family="sans", face="bold"),
        axis.text=element_text(size=20, family="sans"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = 'right') +
  labs(x = "Fork length (mm)", y = "Prob. of maturation")
g

# save
ggsave(file = "./output/figure_logistic_early_mature-fl+treatment.png", 
       plot = g, dpi = 200, width = 9, height = 5)


## 3-4. GLM: early_mature ~ fl * cohort ----------------------------------------

# filtering
d2 <- mr_data %>%
  filter(!(cohort == 2016 & treatment == "early"),
         !(cohort == 2016 & treatment == "late"),
         !(cohort == 2018 & treatment == "early"),
         !(cohort == 2018 & treatment == "late"),
         between(month, 9,11),
         year == cohort + 1) %>% 
  mutate(cohort = as.factor(cohort),
         early_mature = as.numeric(early_mature))

# check data
count(d2, cohort, early_mature)

### 3-4-1. model selection -----------------------------------------------------
# AIC
models <- glm(early_mature ~ fl * cohort, 
              data=d2, 
              family=binomial(link="logit"))

model.list <- dredge(models,rank="AIC")
model.list

# best model
best.model <- get.models(model.list, subset = 1)[[1]]
summary(best.model)

# ANOVA
# analysis of deviance

# null model
model.null <- glm(early_mature ~ 1, 
                  data=d2, family=binomial(link="logit"))

# full model
model.a <- glm(early_mature ~ fl * cohort, 
               data=d2, family=binomial(link="logit"))

anova(model.a, test = "Chisq")

# remove co-efficient
model.b <- glm(early_mature ~ fl + cohort, 
               data=d2, family=binomial(link="logit"))

anova(model.b, test = "Chisq")

# compare models
anova(model.a, model.b, test = "Chisq")
anova(model.null, model.b, test = "Chisq")

# best model
result <- glm(early_mature ~ fl + cohort, 
              data=d2, family=binomial(link="logit"))

# odds ratio
odds <- as.data.frame(exp(result$coefficients))
odds$odds <- formatC(odds$`exp(result$coefficients)`, format = "fg")
odds

# summary to data frame
summary_table <- as.data.frame(summary(result)[[12]])
summary_table <- summary_table %>%
  mutate(Variable = rownames(summary_table)) %>% 
  dplyr::select(., "Variable", "Estimate", SE = "Std. Error", z = "z value",  p = "Pr(>|z|)") %>% 
  mutate(Estimate = round(Estimate, 3),
         SE = round(SE, 3),
         z = round(z, 3),
         p = round(p, 3),
         Signif = case_when(p < 0.001 ~ "**",
                            p <= 0.05 ~ "*",
                            TRUE ~ ""),
         p = as.character(p),
         p = case_when(Signif == "**" ~ "<0.001",
                       Signif != "**" ~ p)) %>% 
  rename("z value" = z,
         "p value" = p,
         "Signif. code" = Signif)

# make table
table <-
  summary_table %>%
  flextable::regulartable() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "sans", part = "all") %>%
  flextable::autofit() %>% 
  flextable::bold(part = "header") %>% 
  flextable::align(align = "right", part = "header") %>% 
  flextable::align(align = "right", part = "body") %>% 
  flextable::bg(bg = 'white', part = "all")
table

# output
save_as_image(table, path = "output/summary_glm_early_mature-fl+cohort.png")
write.csv(summary_table, "output/summary_glm_early_mature-fl+cohort.csv")

### 3-4-2. plot ----------------------------------------------------------------

## model
model <- glm(early_mature ~ fl + cohort, 
             data=d2,
             family=binomial(link="logit"))

# setting
cols = brewer.pal(5, "Set1")

# prediction
preds <- ggpredict(model, terms = c("fl[all]","cohort"), ci_level = 0.95)

# plot prediction

preds <- tibble(
  fl = as.numeric(preds$x),
  cohort = as.character(preds$group),
  predicted = preds$predicted,
  conf.low = preds$conf.low,
  conf.high = preds$conf.high
)

# 50%  probability of maturing
len50.2016 <- -(model$coefficients["(Intercept)"])/model$coefficients["fl"]
len50.2018 <- -(model$coefficients["(Intercept)"] + model$coefficients["cohort2018"])/model$coefficients["fl"]
len50.2019 <- -(model$coefficients["(Intercept)"] + model$coefficients["cohort2019"])/model$coefficients["fl"]
len50.2020 <- -(model$coefficients["(Intercept)"] + model$coefficients["cohort2020"])/model$coefficients["fl"]
len50.2021 <- -(model$coefficients["(Intercept)"] + model$coefficients["cohort2021"])/model$coefficients["fl"]

# plot
g <- 
  ggplot(NULL) +
  geom_jitter(data = d2, aes(x = fl, y = early_mature, colour = cohort), 
              size = 2, height=0, width = 0.2, alpha=0.3) +
  # prediction
  #geom_ribbon(data = preds, aes(x = fl, xend = fl, y = conf.low, yend = conf.high, colour = cohort), alpha=0.1) +
  geom_line(data = preds, aes(x = fl, y = predicted, colour = cohort), linewidth = 2) +
  
  # lines of 50% probability of maturing
  geom_vline(xintercept=len50.2016, linewidth=0.8, linetype="solid", color = cols[1]) +
  geom_vline(xintercept=len50.2018, linewidth=0.8, linetype="solid", color = cols[2]) +
  geom_vline(xintercept=len50.2019, linewidth=0.8, linetype="solid", color = cols[3]) +
  geom_vline(xintercept=len50.2020, linewidth=0.8, linetype="solid", color = cols[4]) +
  geom_vline(xintercept=len50.2021, linewidth=0.8, linetype="solid", color = cols[5]) +
  
  scale_color_manual('cohort: ', values = cols) +
  scale_x_continuous(limits = c(min(d2$fl)-10, max(d2$fl)+10)) +
  scale_y_continuous(breaks=seq(0,1.5,0.5)) +
  theme_bw() +
  theme(axis.title=element_text(size=30, family="sans", face="bold"),
        axis.text=element_text(size=20, family="sans"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = 'right') +
  labs(x = "Fork length (mm)", y = "Prob. of maturation")
g

# save
ggsave(file = "./output/figure_logistic_early_mature-fl+cohort.png", 
       plot = g, dpi = 200, width = 9, height = 5)


### 3-4-3. grouping for model selection ----------------------------------------

# input data
indata <- d2 %>% 
  mutate(subgr = case_when(cohort == "2016" ~ 1,
                           cohort == "2018" ~ 2,
                           cohort == "2019" ~ 3,
                           cohort == "2020" ~ 4,
                           cohort == "2021" ~ 5,
  ),
  subgr = as.factor(subgr),
  y = as.numeric(early_mature),
  x = as.numeric(fl))

# check combination
listParts(length(unique(indata$subgr)))
level <- levels(indata$subgr)

# read model for apartition file
source("code/model_for_apartition_binominal_glm_subgr+x.R")

subsets <- listParts(length(level))
partseq <- 2:length(subsets)
result <- purrr::map(partseq, glm4apartition)

# Best model
bestresult <- result[[which.min(purrr::map(result, function(ap) ap$model$aic))]]
best_subset <- subsets[bestresult$partition]
best_subset
bestmodel <- bestresult$model
summary(bestmodel)

sorted_indices <- order(purrr::map_dbl(result, ~.x$model$aic))

# Second model
second_result <- result[[sorted_indices[2]]]
second_subset <- subsets[second_result$partition]
second_model <- second_result$model

# AIC
null <- glm(y ~ subgr + x,
            family = binomial(link="logit"),
            data = indata)
est = c(AIC(null, bestmodel, second_model)$AIC)
models_data_sub <- data.frame(
  "lank" = c("NULL",
             "1",
             "2"),
  "group" = c("NULL", as.character(best_subset), as.character(second_subset)),
  "model" = "cohort*fl",
  "AIC" = c(est)
)
models_data_sub

# summary to data frame
summary_table <- as.data.frame(summary(bestmodel)[[12]])
summary_table <- summary_table %>%
  mutate(" " = rownames(summary_table)) %>% 
  dplyr::select(., " ", "Estimate", SE = "Std. Error", z = "z value",  p = "Pr(>|z|)") %>% 
  mutate(Estimate = round(Estimate, 3),
         SE = round(SE, 3),
         z = round(z, 3),
         p = round(p, 3),
         Signif = case_when(p < 0.001 ~ "**",
                            p <= 0.05 ~ "*",
                            TRUE ~ ""),
         p = as.character(p),
         p = case_when(Signif == "**" ~ "<0.001",
                       Signif != "**" ~ p)) %>% 
  rename("z value" = z,
         "p value" = p,
         "Signif. code" = Signif)

# make table
table <-
  summary_table %>%
  flextable::regulartable() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "sans", part = "all") %>%
  flextable::autofit() %>% 
  flextable::bold(part = "header") %>% 
  flextable::align(align = "right", part = "body") %>% 
  flextable::bg(bg = 'white', part = "all")
table

# output
save_as_image(table, path = "output/summary_glm_early_mature-fl+cohort_best_grouping_cohort.png")


## 3-5. GLM: early_mature ~ fl + cohort + treatment ----------------------------

# filtering
d2 <- mr_data %>%
  filter(!(cohort == 2016 & treatment == "early"),
         !(cohort == 2016 & treatment == "late"),
         cohort != 2018,
         between(month, 9,11),
         year == cohort + 1) %>% 
  mutate(cohort = as.factor(cohort))

### 3-5-1 logistic regression --------------------------------------------------
# model
result <- glm(early_mature ~ fl + cohort + treatment, data=d2, family=binomial(link="logit"))
summary(result)

# odds ratio
odds <- as.data.frame(exp(result$coefficients))
odds$odds <- formatC(odds$`exp(result$coefficients)`, format = "fg")
odds

# summary to data frame
summary_table <- as.data.frame(summary(result)[[12]])
summary_table <- summary_table %>%
  mutate(Variable = rownames(summary_table)) %>% 
  dplyr::select(., "Variable", "Estimate", SE = "Std. Error", z = "z value",  p = "Pr(>|z|)") %>% 
  mutate(Estimate = round(Estimate, 3),
         SE = round(SE, 3),
         z = round(z, 3),
         p = round(p, 3),
         Signif = case_when(p < 0.001 ~ "**",
                            p <= 0.05 ~ "*",
                            TRUE ~ ""),
         p = as.character(p),
         p = case_when(Signif == "**" ~ "<0.001",
                       Signif != "**" ~ p)) %>% 
  rename("z value" = z,
         "p value" = p,
         "Signif. code" = Signif)

# make table
table <-
  summary_table %>%
  flextable::regulartable() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "sans", part = "all") %>%
  flextable::autofit() %>% 
  flextable::bold(part = "header") %>% 
  flextable::align(align = "right", part = "header") %>% 
  flextable::align(align = "right", part = "body") %>% 
  flextable::bg(bg = 'white', part = "all")
table

# output
save_as_image(table, path = "output/summary_glm_early_mature-fl+cohort+treatment.png")

### 3-5-2 plot -----------------------------------------------------------------

## model
model <- glm(early_mature ~ fl + cohort + treatment, data=d2, family=binomial(link="logit"))

# prediction
preds <- ggpredict(model, terms = c("fl[all]","treatment", "cohort"), ci_level = 0.95)

# plot prediction

preds <- tibble(
  fl = as.numeric(preds$x),
  treatment = as.character(preds$group),
  cohort = as.character(preds$facet),
  predicted = preds$predicted,
  conf.low = preds$conf.low,
  conf.high = preds$conf.high
)


# 50%  probability of maturing
len50.2016c <- -(model$coefficients["(Intercept)"])/model$coefficients["fl"]
len50.2019c <- -(model$coefficients["(Intercept)"] + model$coefficients["cohort2019"])/model$coefficients["fl"]
len50.2020c <- -(model$coefficients["(Intercept)"] + model$coefficients["cohort2020"])/model$coefficients["fl"]
len50.2021c <- -(model$coefficients["(Intercept)"] + model$coefficients["cohort2021"])/model$coefficients["fl"]
len50.2019e <- -(model$coefficients["(Intercept)"] + model$coefficients["cohort2019"] + model$coefficients["treatmentearly"])/model$coefficients["fl"]
len50.2020e <- -(model$coefficients["(Intercept)"] + model$coefficients["cohort2020"] + model$coefficients["treatmentearly"])/model$coefficients["fl"]
len50.2021e <- -(model$coefficients["(Intercept)"] + model$coefficients["cohort2021"] + model$coefficients["treatmentearly"])/model$coefficients["fl"]
len50.2019l <- -(model$coefficients["(Intercept)"] + model$coefficients["cohort2019"] + model$coefficients["treatmentlate"])/model$coefficients["fl"]
len50.2020l <- -(model$coefficients["(Intercept)"] + model$coefficients["cohort2020"] + model$coefficients["treatmentlate"])/model$coefficients["fl"]
len50.2021l <- -(model$coefficients["(Intercept)"] + model$coefficients["cohort2021"] + model$coefficients["treatmentlate"])/model$coefficients["fl"]

# setting
cols = brewer.pal(4, "Spectral")
lines = c("solid","dashed","dotdash","dotted")

# plot
g <- 
  ggplot(NULL) +
  # lines of 50% probability of maturing
  geom_vline(xintercept=len50.2016c, linewidth=0.5, linetype="solid", color = cols[1]) +
  geom_vline(xintercept=len50.2019c, linewidth=0.5, linetype="solid", color = cols[2]) +
  geom_vline(xintercept=len50.2020c, linewidth=0.5, linetype="solid", color = cols[3]) +
  geom_vline(xintercept=len50.2021c, linewidth=0.5, linetype="solid", color = cols[4]) +
  geom_vline(xintercept=len50.2019e, linewidth=0.5, linetype="solid", color = cols[2]) +
  geom_vline(xintercept=len50.2020e, linewidth=0.5, linetype="solid", color = cols[3]) +
  geom_vline(xintercept=len50.2021e, linewidth=0.5, linetype="solid", color = cols[4]) +
  geom_vline(xintercept=len50.2019l, linewidth=0.5, linetype="solid", color = cols[2]) +
  geom_vline(xintercept=len50.2020l, linewidth=0.5, linetype="solid", color = cols[3]) +
  geom_vline(xintercept=len50.2021l, linewidth=0.5, linetype="solid", color = cols[4]) +
  # prediction
  geom_line(data = preds, aes(x = fl, y = predicted, colour = cohort, linetype = treatment), linewidth = 1) +
  
  scale_color_manual('cohort: ', values = cols) +
  scale_linetype_manual('treatment: ', values = lines) +
  scale_x_continuous(limits = c(min(d2$fl)-10, max(d2$fl)+10)) +
  scale_y_continuous(limits = c(-0.05,1.05), breaks=seq(0,1.5,0.5)) +
  theme_bw() +
  theme(axis.title=element_text(size=18, family="sans", face="bold"),
        axis.text=element_text(size=16, family="sans"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = 'top') +
  labs(x = "Fork length (mm)", y = "Prob. of early maturation")
g

# save
ggsave(file = "./output/figure_logistic_early_mature-fl+cohort+treatment.png", plot = g, dpi = 100, width = 8, height = 6)


### 3-4-3 GLM: FL ~ density ----------------------------------------------------
## The density data is made by other codes.

# read density data
dens0 <- read.csv("./output/WRFS_density.csv")

# Pivot aggregation
# total area
area <- dens0 %>% 
  filter(year == 2016) %>% 
  aggregate(total.area ~ treatment, FUN = sum)
# density
dens1 <- aggregate(Nc ~ year + treatment, data = dens0, FUN = sum) %>% 
  left_join(area, by = "treatment") %>% 
  mutate(density = Nc / total.area)

# make data frame
len50.2016c <- data.frame(cohort = "2016", treatment = "control", len50 = len50.2016c)
len50.2019c <- data.frame(cohort = "2019", treatment = "control", len50 = len50.2019c)
len50.2020c <- data.frame(cohort = "2020", treatment = "control", len50 = len50.2020c)
len50.2021c <- data.frame(cohort = "2021", treatment = "control", len50 = len50.2021c)
len50.2019e <- data.frame(cohort = "2019", treatment = "early", len50 = len50.2019e)
len50.2020e <- data.frame(cohort = "2020", treatment = "early", len50 = len50.2020e)
len50.2021e <- data.frame(cohort = "2021", treatment = "early", len50 = len50.2021e)
len50.2019l <- data.frame(cohort = "2019", treatment = "late", len50 = len50.2019l)
len50.2020l <- data.frame(cohort = "2020", treatment = "late", len50 = len50.2020l)
len50.2021l <- data.frame(cohort = "2021", treatment = "late", len50 = len50.2021l)
# merge data
len50 <- rbind(len50.2016c,
               len50.2019c,
               len50.2020c,
               len50.2021c,
               len50.2019e,
               len50.2020e,
               len50.2021e,
               len50.2019l,
               len50.2020l,
               len50.2021l)
# reset rawnames
rownames(len50) <- NULL

dens2 <- len50 %>% 
  mutate(year = as.integer(cohort)) %>% 
  left_join(dens1, by = c("year","treatment"))

str(dens2)

# GLM
result <- glm(len50 ~ density, data=dens2, family = Gamma(link = "identity"))
summary(result)

# prediction
# x values
x <- seq(min(dens2$density), max(dens2$density), length = 1000)
pred <- predict(result, data.frame(density = x), type="response", se.fit = TRUE)

# make data frame
pred <- data.frame(x = x, y = pred$fit, ci = pred$se.fit)

# regression coefficients
r2 <- round(result$coefficients["density"],1)

# plot
g <- pred %>% 
  ggplot(aes(x = x, y = y)) +
  geom_jitter(data = dens2, aes(x = density, y = len50, color = cohort), size = 3) +
  geom_ribbon(aes(ymin = y-ci, ymax = y+ci), alpha = 0.2) +
  geom_line(linewidth = 3) +
  annotate("text", x=0.15,   y=150, size=6, family="sans",
           label = bquote(italic(r^2) == .(r2))) +
  theme_bw() +
  scale_color_manual('cohort: ', values = cols) +
  theme(axis.title=element_text(size=18, family="sans", face="bold"),
        axis.text=element_text(size=16, family="sans"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = 'top') +
  labs(x = "density", y = "threshold size at maturity (mm)")
g

# save
ggsave(file = "./output/figure_glm_thresholdsize-density.png", plot = g, dpi = 150, width = 8, height = 6)


## 3-6. chi-square test: treatment vs fast_growth------------------------------

# input data
data <- d1 %>% 
  mutate(cohort = as.numeric(as.character(cohort))) %>% 
  filter(growth_cluster != "nd",
         cohort <= 2018) %>% 
  mutate(fast_growth = case_when(ocluster_rename == 1 ~ '1',
                                 ocluster_rename == 2 ~ '0',
                                 ocluster_rename == 3 ~ '0',
                                 ocluster_rename == 4 ~ '0',
                                 TRUE ~ "nd"),
         fast_growth = as.factor(fast_growth))

# cross table
table <- table(data$treatment,data$fast_growth)
table

# test
chi_out <- chisq.test(table)
chi_out

# pairwise test
pairwise_result <- pairwiseNominalIndependence(table,
                                               compare = "row",
                                               fisher = FALSE,
                                               gtest = FALSE,
                                               chisq = TRUE,
                                               method = "bonferroni", # The method to adjust multiple p-values
                                               correct = "none",
                                               cramer = FALSE,
                                               digits = 3) # The number of significant digits in the output

# make table
pairwise_result <- pairwise_result %>%
  mutate(p.Chisq = round(p.Chisq, 3),
         p.adj.Chisq = round(p.adj.Chisq, 3),
         Signif = case_when(p.adj.Chisq < 0.001 ~ "**",
                            p.adj.Chisq <= 0.05 ~ "*",
                            TRUE ~ ""),
         p.adj.Chisq = as.character(p.adj.Chisq),
         p.adj.Chisq = case_when(Signif == "**" ~ "<0.001",
                                 Signif != "**" ~ p.adj.Chisq)) %>% 
  rename("Signif. code" = Signif)

table <-
  pairwise_result %>%
  flextable::regulartable() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "sans", part = "all") %>%
  flextable::autofit() %>% 
  flextable::bold(part = "header") %>% 
  flextable::align(align = "right", part = "header") %>% 
  flextable::align(align = "right", part = "body") %>% 
  flextable::bg(bg = 'white', part = "all")
table

# output
save_as_image(table, path = "output/summary_chisq_fast_growth_treatment.png")


## 3-7. chi-square test: cohort vs fast_growth --------------------------------
# input data
data <- d1 %>% 
  mutate(cohort = as.numeric(as.character(cohort))) %>% 
  filter(growth_cluster != "nd",
         !(cohort == 2016 & treatment == "early"),
         !(cohort == 2016 & treatment == "late"),
         !(cohort == 2018 & treatment == "early"),
         !(cohort == 2018 & treatment == "late")) %>% 
  mutate(fast_growth = case_when(ocluster_rename == 1 ~ '1',
                                 ocluster_rename == 2 ~ '0',
                                 ocluster_rename == 3 ~ '0',
                                 ocluster_rename == 4 ~ '0',
                                 TRUE ~ "nd"),
         fast_growth = as.factor(fast_growth))

# cross table
table <- table(data$cohort, data$fast_growth)
table

# test
chi_out <- chisq.test(table)
chi_out

# pairwise test
pairwise_result <- pairwiseNominalIndependence(table,
                                               compare = "row",
                                               fisher = TRUE,
                                               gtest = FALSE,
                                               chisq = FALSE,
                                               method = "bonferroni", # The method to adjust multiple p-values
                                               correct = "none",
                                               cramer = FALSE,
                                               digits = 3) # The number of significant digits in the output

# make table
pairwise_result <- pairwise_result %>%
  mutate(p.Fisher = round(p.Fisher, 3),
         p.adj.Fisher = round(p.adj.Fisher, 3),
         Signif = case_when(p.adj.Fisher < 0.001 ~ "**",
                            p.adj.Fisher <= 0.05 ~ "*",
                            TRUE ~ ""),
         p.adj.Fisher = as.character(p.adj.Fisher),
         p.adj.Fisher = case_when(Signif == "**" ~ "<0.001",
                                 Signif != "**" ~ p.adj.Fisher)) %>% 
  rename("Signif. code" = Signif)

table <-
  pairwise_result %>%
  flextable::regulartable() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "sans", part = "all") %>%
  flextable::autofit() %>% 
  flextable::bold(part = "header") %>% 
  flextable::align(align = "right", part = "header") %>% 
  flextable::align(align = "right", part = "body") %>% 
  flextable::bg(bg = 'white', part = "all")
table

# output
save_as_image(table, path = "output/summary_fisher_fast_growth_cohort.png")


## 3-8. GLM: early_mature ~ fl * sex -------------------------------------------

# filtering
d2 <- mr_data %>%
  filter(between(month, 9,11),
         year == cohort + 1,
         sex %in% c("male","male?","female","female?")) %>% 
  mutate(cohort = as.factor(cohort),
         early_mature = as.numeric(early_mature),
         sex2 = case_when(sex == "male?" ~ "male",
                          sex == "female?" ~ "female",
                          TRUE ~ sex),
         sex2 = as.factor(sex2))

# check data
count(d2, sex2, early_mature)

### 3-8-1. model selection -----------------------------------------------------
# ANOVA
# analysis of deviance

# null model
model.null <- glm(early_mature ~ 1, 
                  data=d2, family=binomial(link="logit"))

# full model
model.a <- glm(early_mature ~ fl * sex2, 
               data=d2, family=binomial(link="logit"))

anova(model.a, test = "Chisq")

# remove co-efficient
model.b <- glm(early_mature ~ fl + sex2, 
               data=d2, family=binomial(link="logit"))

anova(model.b, test = "Chisq")

# compare models
anova(model.a, model.b, test = "Chisq")

# best model
result <- glm(early_mature ~ fl * sex2, 
              data=d2, family=binomial(link="logit"))

# summary to data frame
summary_table <- as.data.frame(summary(result)[[12]])
summary_table <- summary_table %>%
  mutate(Variable = rownames(summary_table)) %>% 
  dplyr::select(., "Variable", "Estimate", SE = "Std. Error", z = "z value",  p = "Pr(>|z|)") %>% 
  mutate(Estimate = round(Estimate, 3),
         SE = round(SE, 3),
         z = round(z, 3),
         p = round(p, 3),
         Signif = case_when(p < 0.001 ~ "**",
                            p <= 0.05 ~ "*",
                            TRUE ~ ""),
         p = as.character(p),
         p = case_when(Signif == "**" ~ "<0.001",
                       Signif != "**" ~ p)) %>% 
  rename("z value" = z,
         "p value" = p,
         "Signif. code" = Signif) %>% 
mutate(Variable = gsub(Variable, pattern = "sex2", replacement = ""))

# make table
table <-
  summary_table %>%
  flextable::regulartable() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "sans", part = "all") %>%
  flextable::autofit() %>% 
  flextable::bold(part = "header") %>% 
  flextable::align(align = "right", part = "header") %>% 
  flextable::align(align = "right", part = "body") %>% 
  flextable::bg(bg = 'white', part = "all")
table

# output
save_as_image(table, path = "output/summary_glm_early_mature-fl*sex.png")
write.csv(summary_table, "output/summary_glm_early_mature-fl*sex.csv")

### 3-8-2. plot ----------------------------------------------------------------

## model
model <- glm(early_mature ~ fl * sex2, 
             data=d2,
             family=binomial(link="logit"))

# setting
cols = brewer.pal(4, "RdYlBu")
cols = c("male"=cols[4], "female"=cols[1])

# prediction
preds <- ggpredict(model, terms = c("fl[all]","sex2"), ci_level = 0.95)

# plot prediction

preds <- tibble(
  fl = as.numeric(preds$x),
  sex2 = as.character(preds$group),
  predicted = preds$predicted,
  conf.low = preds$conf.low,
  conf.high = preds$conf.high
)

# 50%  probability of maturing
len50.female <- -(model$coefficients["(Intercept)"])/model$coefficients["fl"]
len50.male <- -(model$coefficients["(Intercept)"] + model$coefficients["sex2male"])/(model$coefficients["fl"] + model$coefficients["fl:sex2male"])

# plot
g <- 
  ggplot(NULL) +
  geom_jitter(data = d2, aes(x = fl, y = early_mature, colour = sex2), 
              size = 2, height=0, width = 0.2, alpha=0.3) +
  # prediction
  geom_ribbon(data = preds, aes(x = fl, y = predicted,
                                xmin = fl, xmax = fl, 
                                ymin = conf.low, ymax = conf.high, 
                                fill = sex2), alpha=0.2, 
              show.legend = FALSE) +
  geom_line(data = preds, aes(x = fl, y = predicted, colour = sex2), linewidth = 2) +
  
  # lines of 50% probability of maturing
  geom_vline(xintercept=len50.male, linewidth=0.8, linetype="solid", color = cols[1]) +
  geom_vline(xintercept=len50.female, linewidth=0.8, linetype="solid", color = cols[2]) +
  
  # add test
  annotate("text", x=100, y=0.75, 
           label=paste0("male = ", as.character(round(len50.male)), " mm"),
           size = 6, family="sans") +
  annotate("text", x=100, y=0.65, 
           label=paste0("female = ", as.character(round(len50.female)), " mm"),
           size = 6, family="sans") +
  
  scale_color_manual('sex: ', values = cols) +
  scale_x_continuous(limits = c(min(d2$fl)-10, max(d2$fl)+10)) +
  scale_y_continuous(breaks=seq(0,1.5,0.5)) +
  theme_bw() +
  theme(axis.title=element_text(size=30, family="sans", face="bold"),
        axis.text=element_text(size=20, family="sans"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = 'top') +
  labs(x = "Fork length (mm)", y = "Prob. of maturation")
g

# save
ggsave(file = "./output/figure_logistic_early_mature-fl*sex.png", 
       plot = g, dpi = 200, width = 7, height = 5)


## 3-9. GLM: early_mature ~ growth_cluster * sex -------------------------------

# input data
data <- d1 %>% 
  filter(early_mature != "nd",
         growth_cluster != "nd",
         cohort != "2017",
         !(cohort == 2016 & treatment == "early"),
         !(cohort == 2016 & treatment == "late"),
         sex %in% c("male","male?","female","female?")) %>% 
  mutate(early_mature = as.numeric(early_mature),
         growth_cluster = as.factor(growth_cluster),
         cohort = as.factor(cohort),
         sex2 = case_when(sex == "male?" ~ "male",
                          sex == "female?" ~ "female",
                          TRUE ~ sex),
         sex2 = as.factor(sex2))

# check data
count(data, sex2, early_mature)

### 3-9-1. model selection -----------------------------------------------------
# ANOVA
# analysis of deviance

# null model
model.null <- glm(early_mature ~ 1, 
                  data=data, family=binomial(link="logit"))

# full model
model.a <- glm(early_mature ~ growth_cluster * sex2, 
               data=data, family=binomial(link="logit"))

anova(model.a, test = "Chisq")

# remove co-efficient
model.b <- glm(early_mature ~ growth_cluster + sex2, 
               data=data, family=binomial(link="logit"))

anova(model.b, test = "Chisq")

# compare models
anova(model.a, model.b, test = "Chisq")

# best model
result <- glm(early_mature ~ growth_cluster * sex2, 
              data=data, family=binomial(link="logit"))
summary(result)

# summary to data frame
summary_table <- as.data.frame(summary(result)[[12]])
summary_table <- summary_table %>%
  mutate(Variable = rownames(summary_table)) %>% 
  dplyr::select(., "Variable", "Estimate", SE = "Std. Error", z = "z value",  p = "Pr(>|z|)") %>% 
  mutate(Estimate = round(Estimate, 3),
         SE = round(SE, 3),
         z = round(z, 3),
         p = round(p, 3),
         Signif = case_when(p < 0.001 ~ "**",
                            p <= 0.05 ~ "*",
                            TRUE ~ ""),
         p = as.character(p),
         p = case_when(Signif == "**" ~ "<0.001",
                       Signif != "**" ~ p)) %>% 
  rename("z value" = z,
         "p value" = p,
         "Signif. code" = Signif) %>% 
  mutate(Variable = gsub(Variable, pattern = "sex2", replacement = ""))

# make table
table <-
  summary_table %>%
  flextable::regulartable() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "sans", part = "all") %>%
  flextable::autofit() %>% 
  flextable::bold(part = "header") %>% 
  flextable::align(align = "right", part = "header") %>% 
  flextable::align(align = "right", part = "body") %>% 
  flextable::bg(bg = 'white', part = "all")
table

# output
save_as_image(table, path = "output/summary_glm_early_mature-growth_cluster*sex.png")
write.csv(summary_table, "output/summary_glm_early_mature-growth_cluster*sex.csv")


## 3-10. GLM: age1mature ~ growth_cluster * treatment---------------------------

# input data
data <- d1 %>% 
  filter(age1_mature == 0 | age1_mature == 1,
         growth_cluster != "nd",
         cohort == 2016) %>% 
  mutate(treatment = factor(treatment, levels = c("control", "early", "late")))

count(data,age1_mature, growth_cluster, treatment)

### 3-10-1. model selection ----------------------------------------------------

full.model <- glm(age1_mature ~ growth_cluster * treatment, 
              data=data, 
              family=binomial(link="logit"))
summary(full.model)
anova(full.model, test = "Chisq")

model.list <- dredge(full.model, rank="AIC")
model.list

### apartition
indata <- data %>% 
  mutate(subgr = case_when(treatment == "control" ~ 1,
                           treatment == "early" ~ 2,
                           treatment == "late" ~ 3),
  subgr = as.factor(subgr),
  y = as.numeric(age1_mature),
  x = as.factor(growth_cluster))

# check combination
listParts(length(unique(indata$subgr)))
level <- levels(indata$subgr)

# read model for apartition file
source("code/model_for_apartition_binominal_glm_subgr+x.R")

subsets <- listParts(length(level))
partseq <- 2:length(subsets)
result <- purrr::map(partseq, glm4apartition)

# Best model
bestresult <- result[[which.min(purrr::map(result, function(ap) ap$model$aic))]]
best_subset <- subsets[bestresult$partition]
best_subset
bestmodel <- bestresult$model
summary(bestmodel)

sorted_indices <- order(purrr::map_dbl(result, ~.x$model$aic))

# Second model
second_result <- result[[sorted_indices[2]]]
second_subset <- subsets[second_result$partition]
second_model <- second_result$model

# AIC
null <- glm(y ~ subgr + x,
            family = binomial(link="logit"),
            data = indata)
est = c(AIC(null, bestmodel, second_model)$AIC)
models_data_sub <- data.frame(
  "lank" = c("NULL",
             "1",
             "2"),
  "group" = c("NULL", as.character(best_subset), as.character(second_subset)),
  "model" = "growth_cluster+treatment",
  "AIC" = c(est),
  "deltaAIC" = AIC(bestmodel)-c(AIC(null, bestmodel, second_model)$AIC)
)
models_data_sub

# Best cluster: control + late vs early
# but dAIC is 1.996


# edit input data
data <- d1 %>% 
  filter(age1_mature == 0 | age1_mature == 1,
         growth_cluster != "nd",
         cohort == 2016) %>% 
  mutate(treatment2 = case_when(treatment == "control" ~ "control,late",
                                treatment == "late" ~ "control,late",
                                TRUE ~ "early"),
         treatment2 = factor(treatment2, levels = c("control,late", "early")))

count(data,age1_mature, growth_cluster, treatment2)

model <- glm(age1_mature ~ growth_cluster + treatment2, 
                  data=data, 
                  family=binomial(link="logit"))
summary(model)
anova(model, test = "Chisq")

# summary to data frame
summary_table <- as.data.frame(summary(model)[[12]])
summary_table <- summary_table %>%
  mutate(Variable = rownames(summary_table)) %>% 
  dplyr::select(., "Variable", "Estimate", SE = "Std. Error", z = "z value",  p = "Pr(>|z|)") %>% 
  mutate(Estimate = round(Estimate, 3),
         SE = round(SE, 3),
         z = round(z, 3),
         p = round(p, 3),
         Signif = case_when(p < 0.001 ~ "**",
                            p <= 0.05 ~ "*",
                            TRUE ~ ""),
         p = as.character(p),
         p = case_when(Signif == "**" ~ "<0.001",
                       Signif != "**" ~ p)) %>% 
  rename("z value" = z,
         "p value" = p,
         "Signif. code" = Signif) %>% 
  mutate(Variable = gsub(Variable, pattern = "growth_cluster", replacement = "growth cluster-"),
         Variable = gsub(Variable, pattern = "treatment2", replacement = ""))

# make table
table <-
  summary_table %>%
  flextable::regulartable() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "sans", part = "all") %>%
  flextable::autofit() %>% 
  flextable::bold(part = "header") %>% 
  flextable::align(align = "right", part = "header") %>% 
  flextable::align(align = "right", part = "body") %>% 
  flextable::bg(bg = 'white', part = "all")
table

# output
save_as_image(table, path = "output/summary_glm_age1mature-growth_cluster+treatment.png")

### 3-10-2. plot ---------------------------------------------------------------

# setting
cols = c("dimgray","green4")

# prediction
preds <- ggpredict(model, terms = c("growth_cluster","treatment2"), ci_level = 0.95, interval = "confidence")

# plot prediction
preds <- tibble(
  growth_cluster = as.numeric(preds$x),
  treatment = as.character(preds$group),
  predicted = preds$predicted,
  sd = preds$std.error,
  conf.low = preds$conf.low,
  conf.high = preds$conf.high
)

# convert as numeric
data$growth_cluster2 <- as.numeric(data$growth_cluster)

g <- 
  ggplot(NULL) +
  geom_jitter(data = data, aes(x = growth_cluster2, y = age1_mature, colour = treatment2), 
              size = 2, height=0, width = 0.3, alpha=0.3) +
  # prediction
  geom_errorbar(data = preds, 
                aes(x = growth_cluster, y = predicted, 
                    ymin = conf.low, ymax = conf.high,
                    colour = treatment),
                width = 0, linewidth = 1.5, position = position_dodge(0.5)) +
  geom_point(data = preds, aes(x = growth_cluster, y = predicted, colour = treatment),
             size = 5, position = position_dodge(0.5)) +
  
  scale_color_manual('treatment: ', values = cols) +
  scale_x_continuous(limits = c(0.5, 4.5), breaks=seq(0,5,1)) +
  scale_y_continuous(breaks=seq(0,1.5,0.5)) +
  theme_bw() +
  theme(axis.title=element_text(size=30, family="sans", face="bold"),
        axis.text=element_text(size=20, family="sans"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = 'top') +
  labs(x = "Growth cluster", y = "Prob. of maturation")
g

# save
ggsave(file = "./output/figure_logistic_age1mature-growth_cluster+treatment.png", 
       plot = g, dpi = 200, width = 9, height = 8)
write.csv(preds, "./output/predict_prob_age1maturity.csv")

## 3-11. GLM: age2mature ~ growth_cluster * treatment---------------------------

# input data
data <- d1 %>% 
  filter(age2_mature == 0 | age2_mature == 1,
         growth_cluster != "nd",
         cohort == 2015) %>% 
  mutate(treatment = factor(treatment, levels = c("control", "early", "late")),
         age2_mature = as.numeric(age2_mature))

count(data,age2_mature, growth_cluster, treatment)

# Solution for complete separation
# marge growth cluster 1 and 2
data2 <- data %>% 
  mutate(growth_cluster = case_when(growth_cluster == "1" ~ "2",
                                 TRUE ~ growth_cluster))
count(data2,age2_mature, growth_cluster, treatment)

### 3-11-1. model selection ----------------------------------------------------

full.model <- glm(age2_mature ~ growth_cluster * treatment, 
                  data=data2, 
                  family=binomial(link="logit"))
summary(full.model)
anova(full.model, test = "Chisq")

model.list <- dredge(full.model, rank="AIC")
model.list

### apartition
indata <- data %>% 
  mutate(subgr = case_when(treatment == "control" ~ 1,
                           treatment == "early" ~ 2,
                           treatment == "late" ~ 3),
         subgr = as.factor(subgr),
         y = as.numeric(age2_mature),
         x = as.factor(growth_cluster))

# check combination
listParts(length(unique(indata$subgr)))
level <- levels(indata$subgr)

# read model for apartition file
source("code/model_for_apartition_binominal_glm_subgr+x.R")

subsets <- listParts(length(level))
partseq <- 2:length(subsets)
result <- purrr::map(partseq, glm4apartition)

# Best model
bestresult <- result[[which.min(purrr::map(result, function(ap) ap$model$aic))]]
best_subset <- subsets[bestresult$partition]
best_subset
bestmodel <- bestresult$model
summary(bestmodel)

sorted_indices <- order(purrr::map_dbl(result, ~.x$model$aic))

# Second model
second_result <- result[[sorted_indices[2]]]
second_subset <- subsets[second_result$partition]
second_model <- second_result$model

# AIC
null <- glm(y ~ subgr + x,
            family = binomial(link="logit"),
            data = indata)
est = c(AIC(null, bestmodel, second_model)$AIC)
models_data_sub <- data.frame(
  "lank" = c("NULL",
             "1",
             "2"),
  "group" = c("NULL", as.character(best_subset), as.character(second_subset)),
  "model" = "growth_cluster+treatment",
  "AIC" = c(est),
  "deltaAIC" = AIC(bestmodel)-c(AIC(null, bestmodel, second_model)$AIC)
)
models_data_sub

# remove treatment variable
# use 2015 ~ 2021 cohort data
# no data of growth cluster 1 & age 2 immature -> remove
data <- d1 %>% 
  filter(age2_mature == 0 | age2_mature == 1,
         growth_cluster != "nd",
         between(cohort,2015,2021),
         growth_cluster != "1") %>% 
  mutate(age2_mature = as.numeric(age2_mature))

count(data,age2_mature, growth_cluster)

model <- glm(age2_mature ~ growth_cluster, 
             data=data, 
             family=binomial(link="logit"))
summary(model)
anova(model, test = "Chisq")

# summary to data frame
summary_table <- as.data.frame(summary(model)[[12]])
summary_table <- summary_table %>%
  mutate(Variable = rownames(summary_table)) %>% 
  dplyr::select(., "Variable", "Estimate", SE = "Std. Error", z = "z value",  p = "Pr(>|z|)") %>% 
  mutate(Estimate = round(Estimate, 3),
         SE = round(SE, 3),
         z = round(z, 3),
         p = round(p, 3),
         Signif = case_when(p < 0.001 ~ "**",
                            p <= 0.05 ~ "*",
                            TRUE ~ ""),
         p = as.character(p),
         p = case_when(Signif == "**" ~ "<0.001",
                       Signif != "**" ~ p)) %>% 
  rename("z value" = z,
         "p value" = p,
         "Signif. code" = Signif) %>% 
  mutate(Variable = gsub(Variable, pattern = "growth_cluster", replacement = "GC-"),
         Variable = gsub(Variable, pattern = "treatment2", replacement = ""))

# make table
table <-
  summary_table %>%
  flextable::regulartable() %>%
  flextable::theme_booktabs() %>%
  flextable::font(fontname = "sans", part = "all") %>%
  flextable::autofit() %>% 
  flextable::bold(part = "header") %>% 
  flextable::align(align = "right", part = "header") %>% 
  flextable::align(align = "right", part = "body") %>% 
  flextable::bg(bg = 'white', part = "all")
table

# output
save_as_image(table, path = "output/summary_glm_age2mature-growth_cluster.png")

### 3-11-2. plot ---------------------------------------------------------------

# prediction
preds <- ggpredict(model, terms = c("growth_cluster"), ci_level = 0.95, interval = "confidence")

# plot prediction
preds <- tibble(
  growth_cluster = as.numeric(as.character(preds$x)),
  treatment = as.character(preds$group),
  predicted = preds$predicted,
  sd = preds$std.error,
  conf.low = preds$conf.low,
  conf.high = preds$conf.high
) %>% 
  rbind(c(1,1,1,0,0,0)) # growth cluster 1

# convert as numeric
data <- d1 %>% 
  filter(age2_mature == 0 | age2_mature == 1,
         growth_cluster != "nd",
         between(cohort,2015,2021)) %>% 
  mutate(age2_mature = as.numeric(age2_mature),
         growth_cluster2 = as.numeric(growth_cluster))

g <- 
  ggplot(NULL) +
  geom_jitter(data = data, aes(x = growth_cluster2, y = age2_mature), 
              size = 2, height=0, width = 0.3, alpha=0.3) +
  # prediction
  geom_errorbar(data = preds, 
                aes(x = growth_cluster, y = predicted, 
                    ymin = conf.low, ymax = conf.high),
                width = 0, linewidth = 1.5, position = position_dodge(0.5)) +
  geom_point(data = preds, aes(x = growth_cluster, y = predicted),
             size = 5, position = position_dodge(0.5)) +
  
  scale_x_continuous(limits = c(0.5, 4.5), breaks=seq(0,5,1)) +
  scale_y_continuous(breaks=seq(0,1.5,0.5)) +
  theme_bw() +
  theme(axis.title=element_text(size=30, family="sans", face="bold"),
        axis.text=element_text(size=20, family="sans"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = 'top') +
  labs(x = "Growth cluster", y = "Prob. of maturation")
g

# save
ggsave(file = "./output/figure_logistic_age2ature-growth_cluster.png", 
       plot = g, dpi = 200, width = 9, height = 8)
write.csv(preds, "./output/predict_prob_age2maturity.csv")





# end --------------------------------------------------------------------------