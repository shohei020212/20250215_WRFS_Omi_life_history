
# checked by Terui on 2024.07.22
# Modified by Noda on 2025.01.17

# 1. setup ---------------------------------------------------------------------

# package
pacman::p_load(tidyverse, runjags, loo)

# reset environment
rm(list=ls(all.names=T))

# road a model
m <- read.jagsfile("code/model_cjs_r_ver2.R")

# read raw data
d0 <- read.csv("./output/rawdata_filtered.csv") 
gc <- read.csv("./output/k4_data.csv")  

# 2. May.2016 -> Oct.2017 ------------------------------------------------------
## 2-1. prepare input data -----------------------------------------------------

# edit data
data <- d0 %>% 
  left_join(gc, by="id") %>% 
  mutate(growth_cluster = as.numeric(ocluster_rename),
         distance = case_when(between(year, 2016, 2019) ~ (new_landmark - 400) * 10,
                              between(year, 2020, 2023) ~ new_landmark * 10),
         age_at_mature2 = case_when(age_at_mature == "0" ~ "0",
                                    age_at_mature == "1" ~ "1",
                                    age_at_mature == "2" ~ "2",
                                    age_at_mature == "3" ~ "2",
                                    age_at_mature == ">1" ~ "2",
                                    age_at_mature == ">2" ~ "2",
                                    age_at_mature == ">3" ~ "2",
                                    TRUE ~ "nd")) %>% 
  dplyr::select(id, year, month, ym, date, Julian, newsection, treatment, new_landmark, 
                cohort, age_at_mature, age_at_mature2, growth_cluster, distance) %>%
  filter(date >= as.Date('2016-5-01') & date < as.Date('2017-11-01')) %>%
  filter(distance > 0,
         treatment == "early" | treatment == "late" | treatment == "control",
         age_at_mature2 != "0",
         !is.na(growth_cluster),
         cohort == 2016) %>% # extract the data of the autumn survey
  mutate(SecUL = ceiling(distance * 0.1) * 10, # upstream landmark
         SecDL = SecUL - 10, # downstream landmark
         Occasion = paste0("oc", as.numeric(factor(ym)))) %>% # add occasion column
  filter(id %in% names(which(table(id) > 1)) | ym != '2017-10') %>% # remove individuals only captured on last occasion
  filter(id %in% names(which(table(unlist(tapply(id, treatment, unique))) == 1))) %>%　# remove individuals captured in other sections
  distinct(id,year,month, .keep_all = TRUE)

# check cohort error ids
eid <- data %>% 
  filter(year < cohort) %>% 
  pull(id)

# check y-m
occasion_list <- tibble(
  ym = sort(unique(data$ym)),
  oc = sort(unique((data$Occasion)))
)
occasion_list

## 2-2. run CJS model ----------------------------------------------------------
## create history matrices
## - Variable transformation
## - `Occasion` - `year`
## - `CL` - `growth` 
## - `CO` - `cohort` 
## Note: An error will occur if there is data with duplicate ID and captured date.

## CH: Capture history - binary information of recaptured or not
CH <- data %>%
  mutate(capture = 1) %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = capture,
              id_cols = c(id, growth_cluster, treatment),
              values_fill = list(capture = 0)) %>%
  rename(CL = growth_cluster) %>%
  mutate(CL = as.numeric(factor(CL),
                         levels = c("1",
                                    "2",
                                    "3",
                                    "4"))) %>%
  rename(TR = treatment) %>%
  mutate(TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id) # arrange by id

# check number of samples
count <- count(CH, CL, TR)
count

# select sorted occasion data
ch <- CH %>%
  select(sort(unique(data$Occasion)))

# put NA until the first capture
for(i in 1:nrow(ch)) {
  x <- min(which(ch[i,] == 1)) - 1
  if(x != 0) ch[i, 1:x] <- NA
}

CH[, which(colnames(CH) %in% sort(unique(data$Occasion)))] <- ch

print(CH)

## DH: Location history - distance from the downstream end to the upstream landmark of a (re)capture subsection
DH <- data %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = SecUL,
              id_cols = c(id, growth_cluster, treatment)) %>%
  rename(CL = growth_cluster) %>%
  mutate(CL = as.numeric(factor(CL),
                         levels = c("1",
                                    "2",
                                    "3",
                                    "4"))) %>%
  rename(TR = treatment) %>%
  mutate(TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(DH)

## JH: Julian history - Julian date of capture for each occasion
JH <- data %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = Julian,
              id_cols = c(id, growth_cluster, treatment)) %>%
  rename(CL = growth_cluster) %>%
  mutate(CL = as.numeric(factor(CL),
                         levels = c("1",
                                    "2",
                                    "3",
                                    "4"))) %>%
  rename(TR = treatment) %>%
  mutate(TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(JH)


J <- as.matrix(JH[,which(colnames(JH) == "oc1"):ncol(JH)]) # make matrix of JH

### mcmc setting ---------------------------------------------------------------

n.ad <- 100
n.iter <- 2E+4
n.thin <- max(3, ceiling(n.iter/500))
burn <- ceiling(max(10, n.iter/2))
Sample <- ceiling(n.iter/n.thin)


# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

### bayesian inference ---------------------------------------------------------

## Data for JAGS
### Response
X <- as.matrix(DH[,which(colnames(DH) == "oc1"):ncol(DH)] - 5)
Y <- z <- as.matrix(CH[,which(colnames(CH) == "oc1"):ncol(CH)])
z[Y == 0] <- NA # Y data replaced 0 to NA
ObsF <- apply(Y, 1, function(x) min(which(is.na(x) == 0))) # first observed occasion 
ObsL <- apply(Y, 1, function(x) max(which(x == 1))) # last observed occasion 

# put 1 first to last observed occasion 
for(i in 1:length(ObsF)){
  z[i,ObsF[i]:ObsL[i]] <- 1
}

### Explanatory
no_cl <- 1:max(CH$CL)
no_tr <- 1:max(CH$TR)
no_co <- 1:(max(CH$CL)*max(CH$TR))

Nday <- diff(apply(J, 2, median, na.rm = T)) # interval between occasions
Grid <- expand.grid(CL = no_cl, TR = no_tr) # Create unique combinations
Grid$G <- no_co # make new no. column
Grid$id <- as.numeric(paste0(Grid$CL, Grid$TR) ) # make new id column
Grid
# marge Grid and Capture history data (Grouping)
groupdat <- data.frame(id = as.numeric(paste0(CH$CL, CH$TR))) %>%
  left_join(Grid, by = "id")
G <- groupdat$G

# make input data
Djags <- list(X = X,
              Y = Y,
              Nday = Nday,
              z = z,
              Nind = nrow(Y), # Number of sample
              Nt = ncol(Y), # Number of occasion
              Ng = length(unique(G)), # Number of unique combinations for growth and treatment
              G = G, # grouping data
              ObsF = ObsF
)
# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

# Repeat an equation 3 times
inits <- replicate(3,
                   list(logit.p = matrix(0, nrow = Djags$Ng, ncol = Djags$Nt-1),
                        .RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA ),
                   simplify = F)

for(k in 1:3) inits[[k]]$.RNG.seed <- k


### mcmc using runjags ---------------------------------------------------------
# Parameter estimation by MCMC using runjags
post <- run.jags(m$model,
                 monitor = para,
                 data = Djags,
                 n.chains = 3,
                 inits = inits,
                 method = "parallel",
                 burnin = burn,
                 sample = Sample,
                 adapt = n.ad,
                 thin = n.thin,
                 n.sims = 3,
                 modules = "glm")
# out result
bpost <- MCMCvis::MCMCsummary(post$mcmc,probs = c(0.05,0.5,0.95))

# evaluation of posterior distribution
print(max(na.omit(bpost[, "Rhat"])))

# trace plot
x1 <- data.frame(post$mcmc[[1]])$phi.1.8
x2 <- data.frame(post$mcmc[[2]])$phi.1.8
x3 <- data.frame(post$mcmc[[3]])$phi.1.8

par()
plot(x1, type="l", col = 'red', ylim=c(0,1))
par(new=T)
plot(x2, type="l", col = 'green', ylim=c(0,1), ann=F)
par(new=T)
plot(x3, type="l", col = 'blue', ylim=c(0,1), ann=F)


### output ---------------------------------------------------------------------

bpost <- bpost %>% 
  as_tibble(rownames = "param") %>% 
  rename(lower = '5%',
         median = '50%',
         upper = '95%')

## convert date
date <- apply(J, 2, median, na.rm = T) %>% 
  as.Date(origin = as.Date("1970-01-01"))

dat_date <- tibble(occasion = 1:length(date),
                   date = date)
dat_date

# make new data
# phi age0 Oct,2016 -> age1 Oct,2017
result <- bpost %>% 
  filter(str_detect(param, "phi")) %>% 
  mutate(x = str_extract(param, pattern = "\\[.{1,}\\]"),
         y = str_remove_all(x, pattern = "\\[|\\]")) %>%
  separate(y, into = c("group", "occasion"), sep = ",") %>% 
  mutate(group = as.numeric(group),
         occasion = as.numeric(occasion)) %>% 
  left_join(Grid, by = c("group" = "G")) %>% 
  left_join(dat_date, by = "occasion") %>% 
  filter(occasion == length(dat_date$date)) %>% 
  mutate(treatment = case_when(TR == 1 ~ "control",
                               TR == 2 ~ "early",
                               TR == 3 ~ "late"),
         treatment = fct_relevel(treatment, c("early", "late")),
         growth = as.character(CL)) %>% 
  left_join(count, by = c("TR","CL")) %>% 
  mutate(n = case_when(treatment == "early" & CL == 1 ~ paste0('n = ', as.character(n)),
                         TRUE ~ as.character(n)))

# save to csv
write.csv(result, "output/phi_age0_1_2016cohort.csv")

# phi age0 May,2016 -> age0 Oct,2016
result2 <- bpost %>% 
  filter(str_detect(param, "phi")) %>% 
  mutate(x = str_extract(param, pattern = "\\[.{1,}\\]"),
         y = str_remove_all(x, pattern = "\\[|\\]")) %>%
  separate(y, into = c("group", "occasion"), sep = ",") %>% 
  mutate(group = as.numeric(group),
         occasion = as.numeric(occasion)) %>% 
  left_join(Grid, by = c("group" = "G")) %>% 
  left_join(dat_date, by = "occasion") %>% 
  filter(occasion == 4) %>% 
  mutate(treatment = case_when(TR == 1 ~ "control",
                               TR == 2 ~ "early",
                               TR == 3 ~ "late"),
         treatment = fct_relevel(treatment, c("early", "late")),
         growth = as.character(CL)) %>% 
  left_join(count, by = c("TR","CL")) %>% 
  mutate(n = case_when(treatment == "early" & CL == 1 ~ paste0('n = ', as.character(n)),
                       TRUE ~ as.character(n)))

# save to csv
write.csv(result2, "output/phi_age0_2016cohort.csv")

## 2-3. plot -------------------------------------------------------------------

library(RColorBrewer)
# color list
color_list = brewer.pal(4, "RdYlBu")

# phi age0 may -> age1 oct
g <- ggplot(result) +
  geom_errorbar(aes(x = growth, y = mean, ymin = lower, ymax = upper, colour = growth), 
                width = 0, linewidth = 1) +
  geom_point(aes(x = growth, y = mean, colour = growth), size = 4) +
  facet_wrap(facets = ~ treatment, ncol = 3) + 
  ylab("Survival probability") +
  xlab("Growth cluster") +
  theme_bw() +
  scale_color_manual('glowth cluster: ', values = color_list) +
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0,1.5,0.5)) +
  geom_text(aes(x = CL, y = 1.1, label = n), size = 5) +
  theme(strip.text=element_text(size=25, family="sans", face="bold"),
        axis.title=element_text(size=28, family="sans", face="bold"),
        axis.text=element_text(size=25, family="sans"),
        strip.background = element_rect(color="black", fill="white", linewidth=0.6, linetype="solid"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = 'none')
g

# save png format
ggsave(g,
       file = "output/figure_phi_age0_1_2016cohort.png",
       width = 10,
       height = 4)



# phi age0 may -> age0 oct
g <- ggplot(result2) +
  geom_errorbar(aes(x = growth, y = mean, ymin = lower, ymax = upper, colour = growth), 
                width = 0, linewidth = 1) +
  geom_point(aes(x = growth, y = mean, colour = growth), size = 4) +
  facet_wrap(facets = ~ treatment, ncol = 3) + 
  ylab("Survival probability") +
  xlab("Growth cluster") +
  theme_bw() +
  scale_color_manual('glowth cluster: ', values = color_list) +
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0,1.5,0.5)) +
  geom_text(aes(x = CL, y = 1.1, label = n), size = 5) +
  theme(strip.text=element_text(size=25, family="sans", face="bold"),
        axis.title=element_text(size=28, family="sans", face="bold"),
        axis.text=element_text(size=25, family="sans"),
        strip.background = element_rect(color="black", fill="white", linewidth=0.6, linetype="solid"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = 'none')
g

# save png format
ggsave(g,
       file = "output/figure_phi_age0_2016cohort.png",
       width = 10,
       height = 4)

# 3. Oct.2016 -> Oct.2017 ------------------------------------------------------
## 3-1. prepare input data -----------------------------------------------------

# edit data
data <- d0 %>% 
  left_join(gc, by="id") %>% 
  mutate(growth_cluster = as.numeric(ocluster_rename),
         distance = case_when(between(year, 2016, 2019) ~ (new_landmark - 400) * 10,
                              between(year, 2020, 2023) ~ new_landmark * 10),
         age_at_mature2 = case_when(age_at_mature == "0" ~ "0",
                                    age_at_mature == "1" ~ "1",
                                    age_at_mature == "2" ~ "2",
                                    age_at_mature == "3" ~ "2",
                                    age_at_mature == ">1" ~ "2",
                                    age_at_mature == ">2" ~ "2",
                                    age_at_mature == ">3" ~ "2",
                                    TRUE ~ "nd")) %>% 
  dplyr::select(id, year, month, ym, date, Julian, newsection, treatment, new_landmark, 
                cohort, age_at_mature, age_at_mature2, growth_cluster, distance) %>% 
  filter(date >= as.Date('2016-10-01') & date < as.Date('2017-11-01')) %>%
  filter(distance > 0,
         treatment == "early" | treatment == "late" | treatment == "control",
         age_at_mature2 != "0",
         !is.na(growth_cluster),
         cohort == 2016) %>% # extract the data of the autumn survey
  mutate(SecUL = ceiling(distance * 0.1) * 10, # upstream landmark
         SecDL = SecUL - 10, # downstream landmark
         Occasion = paste0("oc", as.numeric(factor(ym)))) %>% # add occasion column
  filter(id %in% names(which(table(id) > 1)) | ym != '2017-10') %>% # remove individuals only captured on last occasion
  filter(id %in% names(which(table(unlist(tapply(id, treatment, unique))) == 1))) %>%　# remove individuals captured in other sections
  distinct(id,year,month, .keep_all = TRUE)

# check cohort error ids
eid <- data %>% 
  filter(year < cohort) %>% 
  pull(id)

# check y-m
occasion_list <- tibble(
  ym = sort(unique(data$ym)),
  oc = sort(unique((data$Occasion)))
)
occasion_list

## 3-2. run CJS model ----------------------------------------------------------
## create history matrices
## - Variable transformation
## - `Occasion` - `year`
## - `CL` - `growth` 
## - `CO` - `cohort` 
## Note: An error will occur if there is data with duplicate ID and captured date.

## CH: Capture history - binary information of recaptured or not
CH <- data %>%
  mutate(capture = 1) %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = capture,
              id_cols = c(id, growth_cluster, treatment),
              values_fill = list(capture = 0)) %>%
  rename(CL = growth_cluster) %>%
  mutate(CL = as.numeric(factor(CL),
                         levels = c("1",
                                    "2",
                                    "3",
                                    "4"))) %>%
  rename(TR = treatment) %>%
  mutate(TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id) # arrange by id

# check number of samples
count <- count(CH, CL, TR)
count

# select sorted occasion data
ch <- CH %>%
  select(sort(unique(data$Occasion)))

# put NA until the first capture
for(i in 1:nrow(ch)) {
  x <- min(which(ch[i,] == 1)) - 1
  if(x != 0) ch[i, 1:x] <- NA
}

CH[, which(colnames(CH) %in% sort(unique(data$Occasion)))] <- ch

print(CH)

## DH: Location history - distance from the downstream end to the upstream landmark of a (re)capture subsection
DH <- data %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = SecUL,
              id_cols = c(id, growth_cluster, treatment)) %>%
  rename(CL = growth_cluster) %>%
  mutate(CL = as.numeric(factor(CL),
                         levels = c("1",
                                    "2",
                                    "3",
                                    "4"))) %>%
  rename(TR = treatment) %>%
  mutate(TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(DH)

## JH: Julian history - Julian date of capture for each occasion
JH <- data %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = Julian,
              id_cols = c(id, growth_cluster, treatment)) %>%
  rename(CL = growth_cluster) %>%
  mutate(CL = as.numeric(factor(CL),
                         levels = c("1",
                                    "2",
                                    "3",
                                    "4"))) %>%
  rename(TR = treatment) %>%
  mutate(TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(JH)

J <- as.matrix(JH[,which(colnames(JH) == "oc1"):ncol(JH)]) # make matrix of JH

### mcmc setting ---------------------------------------------------------------

n.ad <- 100
n.iter <- 2E+4
n.thin <- max(3, ceiling(n.iter/500))
burn <- ceiling(max(10, n.iter/2))
Sample <- ceiling(n.iter/n.thin)

# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

### bayesian inference ---------------------------------------------------------

## Data for JAGS
### Response
X <- as.matrix(DH[,which(colnames(DH) == "oc1"):ncol(DH)] - 5)
Y <- z <- as.matrix(CH[,which(colnames(CH) == "oc1"):ncol(CH)])
z[Y == 0] <- NA # Y data replaced 0 to NA
ObsF <- apply(Y, 1, function(x) min(which(is.na(x) == 0))) # first observed occasion 
ObsL <- apply(Y, 1, function(x) max(which(x == 1))) # last observed occasion 

# put 1 first to last observed occasion 
for(i in 1:length(ObsF)){
  z[i,ObsF[i]:ObsL[i]] <- 1
}

### Explanatory
no_cl <- 1:max(CH$CL)
no_tr <- 1:max(CH$TR)
no_co <- 1:(max(CH$CL)*max(CH$TR))

Nday <- diff(apply(J, 2, median, na.rm = T)) # interval between occasions
Grid <- expand.grid(CL = no_cl, TR = no_tr) # Create unique combinations
Grid$G <- no_co # make new no. column
Grid$id <- as.numeric(paste0(Grid$CL, Grid$TR) ) # make new id column
Grid
# marge Grid and Capture history data (Grouping)
groupdat <- data.frame(id = as.numeric(paste0(CH$CL, CH$TR))) %>%
  left_join(Grid, by = "id")
G <- groupdat$G

# make input data
Djags <- list(X = X,
              Y = Y,
              Nday = Nday,
              z = z,
              Nind = nrow(Y), # Number of sample
              Nt = ncol(Y), # Number of occasion
              Ng = length(unique(G)), # Number of unique combinations for growth and treatment
              G = G, # grouping data
              ObsF = ObsF
)
# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

# Repeat an equation 3 times
inits <- replicate(3,
                   list(logit.p = matrix(0, nrow = Djags$Ng, ncol = Djags$Nt-1),
                        .RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA ),
                   simplify = F)

for(k in 1:3) inits[[k]]$.RNG.seed <- k

### mcmc using runjags ---------------------------------------------------------
# Parameter estimation by MCMC using runjags
post <- run.jags(m$model,
                 monitor = para,
                 data = Djags,
                 n.chains = 3,
                 inits = inits,
                 method = "parallel",
                 burnin = burn,
                 sample = Sample,
                 adapt = n.ad,
                 thin = n.thin,
                 n.sims = 3,
                 modules = "glm")
# out result
bpost <- MCMCvis::MCMCsummary(post$mcmc,probs = c(0.05,0.5,0.95))

# evaluation of posterior distribution
print(max(na.omit(bpost[, "Rhat"])))

# trace plot
x1 <- data.frame(post$mcmc[[1]])$phi.1.5
x2 <- data.frame(post$mcmc[[2]])$phi.1.5
x3 <- data.frame(post$mcmc[[3]])$phi.1.5

par()
plot(x1, type="l", col = 'red', ylim=c(0,1))
par(new=T)
plot(x2, type="l", col = 'green', ylim=c(0,1), ann=F)
par(new=T)
plot(x3, type="l", col = 'blue', ylim=c(0,1), ann=F)


### output ---------------------------------------------------------------------

bpost <- bpost %>% 
  as_tibble(rownames = "param") %>% 
  rename(lower = '5%',
         median = '50%',
         upper = '95%')

## convert date
date <- apply(J, 2, median, na.rm = T) %>% 
  as.Date(origin = as.Date("1970-01-01"))

dat_date <- tibble(occasion = 1:length(date),
                   date = date)
dat_date

# make new data
# phi age0 Oct,2016 -> age1 Oct,2017
result <- bpost %>% 
  filter(str_detect(param, "phi")) %>% 
  mutate(x = str_extract(param, pattern = "\\[.{1,}\\]"),
         y = str_remove_all(x, pattern = "\\[|\\]")) %>%
  separate(y, into = c("group", "occasion"), sep = ",") %>% 
  mutate(group = as.numeric(group),
         occasion = as.numeric(occasion)) %>% 
  left_join(Grid, by = c("group" = "G")) %>% 
  left_join(dat_date, by = "occasion") %>% 
  filter(occasion == length(dat_date$date)) %>% 
  mutate(treatment = case_when(TR == 1 ~ "control",
                               TR == 2 ~ "early",
                               TR == 3 ~ "late"),
         treatment = fct_relevel(treatment, c("early", "late")),
         growth = as.character(CL)) %>% 
  left_join(count, by = c("TR","CL")) %>% 
  mutate(n = case_when(treatment == "early" & CL == 1 ~ paste0('n = ', as.character(n)),
                       TRUE ~ as.character(n)))

# save to csv
write.csv(result, "output/phi_age0_1_oct_2016cohort.csv")

## 3-3. plot -------------------------------------------------------------------

library(RColorBrewer)
# color list
color_list = brewer.pal(4, "RdYlBu")

# phi age0 oct. -> age1 oct.
g <- ggplot(result) +
  geom_errorbar(aes(x = growth, y = mean, ymin = lower, ymax = upper, colour = growth), 
                width = 0, linewidth = 1) +
  geom_point(aes(x = growth, y = mean, colour = growth), size = 4) +
  facet_wrap(facets = ~ treatment, ncol = 3) + 
  ylab("Survival probability") +
  xlab("Growth cluster") +
  theme_bw() +
  scale_color_manual('glowth cluster: ', values = color_list) +
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0,1.5,0.5)) +
  geom_text(aes(x = CL, y = 1.1, label = n), size = 5) +
  theme(strip.text=element_text(size=25, family="sans", face="bold"),
        axis.title=element_text(size=28, family="sans", face="bold"),
        axis.text=element_text(size=25, family="sans"),
        strip.background = element_rect(color="black", fill="white", linewidth=0.6, linetype="solid"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = 'none')
g

# save png format
ggsave(g,
       file = "output/figure_phi_age0_1_oct_2016cohort.png",
       width = 10,
       height = 4)

# 4. Oct.2017 -> Aug.2018 ------------------------------------------------------
## 4-1. prepare input data -----------------------------------------------------

# edit data
data <- d0 %>% 
  left_join(gc, by="id") %>% 
  mutate(growth_cluster = as.numeric(ocluster_rename),
         distance = case_when(between(year, 2016, 2019) ~ (new_landmark - 400) * 10,
                              between(year, 2020, 2023) ~ new_landmark * 10),
         age_at_mature2 = case_when(age_at_mature == "0" ~ "0",
                                    age_at_mature == "1" ~ "1",
                                    age_at_mature == "2" ~ "2",
                                    age_at_mature == "3" ~ "2",
                                    age_at_mature == ">1" ~ "2",
                                    age_at_mature == ">2" ~ "2",
                                    age_at_mature == ">3" ~ "2",
                                    TRUE ~ "nd")) %>% 
  filter(date >= as.Date('2017-10-01') & date < as.Date('2019-10-01')) %>%
  filter(age1_mature == 0 | age1_mature == 1,
         distance > 0,
         treatment == "early" | treatment == "late" | treatment == "control",
         !is.na(growth_cluster),
         cohort == 2016) %>% # extract the data of the autumn survey
  mutate(SecUL = ceiling(distance * 0.1) * 10, # upstream landmark
         SecDL = SecUL - 10, # downstream landmark
         Occasion = paste0("oc", as.numeric(factor(ym)))) %>% # add occasion column
  filter(id %in% names(which(table(id) > 1)) | ym != '2018-8') %>% # remove individuals only captured on last occasion
  filter(id %in% names(which(table(unlist(tapply(id, treatment, unique))) == 1))) %>%　# remove individuals captured in other sections
  distinct(id,year,month, .keep_all = TRUE)

# check cohort error ids
eid <- data %>% 
  filter(year < cohort) %>% 
  pull(id)

# check y-m
occasion_list <- tibble(
  ym = sort(unique(data$ym)),
  oc = sort(unique((data$Occasion)))
)
occasion_list

## 4-2. run CJS model ----------------------------------------------------------
## create history matrices
## - Variable transformation
## - `Occasion` - `year`
## - `CL` - `growth` 
## - `TR` - `treatment` 
## - `MT` - `age at maturity` 
## Note: An error will occur if there is data with duplicate ID and captured date.

## CH: Capture history - binary information of recaptured or not
CH <- data %>%
  mutate(capture = 1) %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = capture,
              id_cols = c(id, growth_cluster, treatment, age1_mature),
              values_fill = list(capture = 0)) %>%
  rename(CL = growth_cluster,
         TR = treatment,
         MT = age1_mature) %>%
  mutate(CL = as.numeric(factor(CL),
                         levels = c("1",
                                    "2",
                                    "3",
                                    "4")),
         TR = as.numeric(factor(TR, 
                                levels = c("control",
                                           "early",
                                           "late"))),
         MT = as.factor(MT)) %>%
  arrange(id) # arrange by id

# check number of samples
count <- count(CH, CL, TR, MT)
count

# select sorted occasion data
ch <- CH %>%
  select(sort(unique(data$Occasion)))

# put NA until the first capture
for(i in 1:nrow(ch)) {
  x <- min(which(ch[i,] == 1)) - 1
  if(x != 0) ch[i, 1:x] <- NA
}

CH[, which(colnames(CH) %in% sort(unique(data$Occasion)))] <- ch

print(CH)

## DH: Location history - distance from the downstream end to the upstream landmark of a (re)capture subsection
DH <- data %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = SecUL,
              id_cols = c(id, growth_cluster, treatment, age1_mature)) %>%
  rename(CL = growth_cluster,
         TR = treatment,
         MT = age1_mature) %>%
  mutate(CL = as.numeric(factor(CL),
                         levels = c("1",
                                    "2",
                                    "3",
                                    "4")),
         TR = as.numeric(factor(TR, 
                                levels = c("control",
                                           "early",
                                           "late"))),
         MT = as.factor(MT)) %>%
  arrange(id) # arrange by id

print(DH)

## JH: Julian history - Julian date of capture for each occasion
JH <- data %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = Julian,
              id_cols = c(id, growth_cluster, treatment, age1_mature)) %>%
  rename(CL = growth_cluster,
         TR = treatment,
         MT = age1_mature) %>%
  mutate(CL = as.numeric(factor(CL),
                         levels = c("1",
                                    "2",
                                    "3",
                                    "4")),
         TR = as.numeric(factor(TR, 
                                levels = c("control",
                                           "early",
                                           "late"))),
         MT = as.factor(MT)) %>%
  arrange(id) # arrange by id

print(JH)

J <- as.matrix(JH[,which(colnames(JH) == "oc1"):ncol(JH)]) # make matrix of JH

### mcmc setting ---------------------------------------------------------------

n.ad <- 100
n.iter <- 2E+4
n.thin <- max(3, ceiling(n.iter/500))
burn <- ceiling(max(10, n.iter/2))
Sample <- ceiling(n.iter/n.thin)

# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

### bayesian inference ---------------------------------------------------------

## Data for JAGS
### Response
X <- as.matrix(DH[,which(colnames(DH) == "oc1"):ncol(DH)] - 5)
Y <- z <- as.matrix(CH[,which(colnames(CH) == "oc1"):ncol(CH)])
z[Y == 0] <- NA # Y data replaced 0 to NA
ObsF <- apply(Y, 1, function(x) min(which(is.na(x) == 0))) # first observed occasion 
ObsL <- apply(Y, 1, function(x) max(which(x == 1))) # last observed occasion 

# put 1 first to last observed occasion 
for(i in 1:length(ObsF)){
  z[i,ObsF[i]:ObsL[i]] <- 1
}

### Explanatory
Nday <- diff(apply(J, 2, median, na.rm = T)) # interval between occasions
Grid <- count %>% dplyr::select(CL,TR,MT) # Create unique combinations
Grid$G <- 1:length(Grid$CL) # make new no. column
Grid$id <- as.numeric(paste0(Grid$CL,Grid$TR,Grid$MT)) # make new id column

# marge Grid and Capture history data (Grouping)
groupdat <- tibble(id = as.numeric(paste0(CH$CL, CH$TR, CH$MT))) %>%
  left_join(Grid, by = "id")
G <- groupdat$G

# make input data
Djags <- list(X = X,
              Y = Y,
              Nday = Nday,
              z = z,
              Nind = nrow(Y), # Number of sample
              Nt = ncol(Y), # Number of occasion
              Ng = length(unique(G)), # Number of unique combinations for growth and treatment
              G = G, # grouping data
              ObsF = ObsF
)
# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

# Repeat an equation 3 times
inits <- replicate(3,
                   list(logit.p = matrix(0, nrow = Djags$Ng, ncol = Djags$Nt-1),
                        .RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA ),
                   simplify = F)

for(k in 1:3) inits[[k]]$.RNG.seed <- k

### mcmc using runjags ---------------------------------------------------------
# Parameter estimation by MCMC using runjags
post <- run.jags(m$model,
                 monitor = para,
                 data = Djags,
                 n.chains = 3,
                 inits = inits,
                 method = "parallel",
                 burnin = burn,
                 sample = Sample,
                 adapt = n.ad,
                 thin = n.thin,
                 n.sims = 3,
                 modules = "glm")
# out result
bpost <- MCMCvis::MCMCsummary(post$mcmc,probs = c(0.05,0.5,0.95))

# evaluation of posterior distribution
print(max(na.omit(bpost[, "Rhat"])))

# trace plot
x1 <- data.frame(post$mcmc[[1]])$phi.1.5
x2 <- data.frame(post$mcmc[[2]])$phi.1.5
x3 <- data.frame(post$mcmc[[3]])$phi.1.5

par()
plot(x1, type="l", col = 'red', ylim=c(0,1))
par(new=T)
plot(x2, type="l", col = 'green', ylim=c(0,1), ann=F)
par(new=T)
plot(x3, type="l", col = 'blue', ylim=c(0,1), ann=F)


### output ---------------------------------------------------------------------

bpost <- bpost %>% 
  as_tibble(rownames = "param") %>% 
  rename(lower = '5%',
         median = '50%',
         upper = '95%')

## convert date
date <- apply(J, 2, median, na.rm = T) %>% 
  as.Date(origin = as.Date("1970-01-01"))

dat_date <- tibble(occasion = 1:length(date),
                   date = date)
dat_date

# make new data
# phi age1 May,2017 -> age2 July,2018
result <- bpost %>% 
  filter(str_detect(param, "phi")) %>% 
  mutate(x = str_extract(param, pattern = "\\[.{1,}\\]"),
         y = str_remove_all(x, pattern = "\\[|\\]")) %>%
  separate(y, into = c("group", "occasion"), sep = ",") %>% 
  mutate(group = as.numeric(group),
         occasion = as.numeric(occasion)) %>% 
  left_join(Grid, by = c("group" = "G")) %>% 
  left_join(dat_date, by = "occasion") %>% 
  filter(occasion == 4) %>% # 2018-07 oc4  
  mutate(treatment = case_when(TR == 1 ~ "control",
                               TR == 2 ~ "early",
                               TR == 3 ~ "late"),
         treatment = fct_relevel(treatment, c("early", "late")),
         growth = as.character(CL),
         age1_mature = as.character(MT)) %>% 
  left_join(count, by = c("TR","CL","MT")) %>% 
  mutate(n = case_when(treatment == "early" & CL == 1 ~ paste0('n = ', as.character(n)),
                       TRUE ~ as.character(n)))

# save to csv
write.csv(result, "output/phi_age1_2_oct_2016cohort.csv")

## 4-3. plot -------------------------------------------------------------------

library(RColorBrewer)
# color list
color_list = brewer.pal(4, "RdYlBu")

# phi age1 oct. -> age2 aug.
g <- ggplot(result, aes(x = growth, y = mean, ymin = lower, ymax = upper, 
                         colour = growth, group = age1_mature, shape = age1_mature)) +
  geom_errorbar(width = 0, linewidth = 1, position = position_dodge(.7)) +
  geom_point(size = 4, position = position_dodge(.7)) +
  facet_wrap(facets = ~ treatment, ncol = 3) + 
  ylab("Survival probability") +
  xlab("Growth cluster") +
  theme_bw() +
  scale_color_manual('glowth cluster: ', values = color_list) +
  scale_shape_manual('age 1 maturity: ', values = c(16,17)) +
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0,1.5,0.5)) +
  geom_text(aes(x = growth, y = 1.1, label = n), size = 4, position = position_dodge(.7)) +
  theme(strip.text=element_text(size=25, family="sans", face="bold"),
        axis.title=element_text(size=28, family="sans", face="bold"),
        axis.text=element_text(size=25, family="sans"),
        strip.background = element_rect(color="black", fill="white", linewidth=0.6, linetype="solid"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = 'top')
g

# save png format
ggsave(g,
       file = "output/figure_phi_age1_2_oct_2016cohort.png",
       width = 12,
       height = 5)


# 5. May.2016 -> Oct.2017, no growth cluster -----------------------------------
## 2-1. prepare input data -----------------------------------------------------

# edit data
data <- d0 %>% 
  left_join(gc, by="id") %>% 
  mutate(distance = case_when(between(year, 2016, 2019) ~ (new_landmark - 400) * 10,
                              between(year, 2020, 2023) ~ new_landmark * 10),
         age_at_mature2 = case_when(age_at_mature == "0" ~ "0",
                                    age_at_mature == "1" ~ "1",
                                    age_at_mature == "2" ~ "2",
                                    age_at_mature == "3" ~ "2",
                                    age_at_mature == ">1" ~ "2",
                                    age_at_mature == ">2" ~ "2",
                                    age_at_mature == ">3" ~ "2",
                                    TRUE ~ "nd")) %>% 
  dplyr::select(id, year, month, ym, date, Julian, newsection, treatment, new_landmark, 
                cohort, age_at_mature, age_at_mature2, distance) %>%
  filter(date >= as.Date('2016-5-01') & date < as.Date('2017-11-01')) %>%
  filter(distance > 0,
         treatment == "early" | treatment == "late" | treatment == "control",
         age_at_mature2 != "0",
         cohort == 2016) %>% # extract the data of the autumn survey
  mutate(SecUL = ceiling(distance * 0.1) * 10, # upstream landmark
         SecDL = SecUL - 10, # downstream landmark
         Occasion = paste0("oc", as.numeric(factor(ym)))) %>% # add occasion column
  filter(id %in% names(which(table(id) > 1)) | ym != '2017-10') %>% # remove individuals only captured on last occasion
  filter(id %in% names(which(table(unlist(tapply(id, treatment, unique))) == 1))) %>%　# remove individuals captured in other sections
  distinct(id,year,month, .keep_all = TRUE)

# check cohort error ids
eid <- data %>% 
  filter(year < cohort) %>% 
  pull(id)

# check y-m
occasion_list <- tibble(
  ym = sort(unique(data$ym)),
  oc = sort(unique((data$Occasion)))
)
occasion_list

## 2-2. run CJS model ----------------------------------------------------------
## create history matrices
## - Variable transformation
## - `Occasion` - `year`
## - `CL` - `growth` 
## - `CO` - `cohort` 
## Note: An error will occur if there is data with duplicate ID and captured date.

## CH: Capture history - binary information of recaptured or not
CH <- data %>%
  mutate(capture = 1) %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = capture,
              id_cols = c(id, treatment),
              values_fill = list(capture = 0)) %>%
  rename(TR = treatment) %>%
  mutate(TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id) # arrange by id

# check number of samples
count <- count(CH, TR)
count

# select sorted occasion data
ch <- CH %>%
  select(sort(unique(data$Occasion)))

# put NA until the first capture
for(i in 1:nrow(ch)) {
  x <- min(which(ch[i,] == 1)) - 1
  if(x != 0) ch[i, 1:x] <- NA
}

CH[, which(colnames(CH) %in% sort(unique(data$Occasion)))] <- ch

print(CH)

## DH: Location history - distance from the downstream end to the upstream landmark of a (re)capture subsection
DH <- data %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = SecUL,
              id_cols = c(id, treatment)) %>%
  rename(TR = treatment) %>%
  mutate(TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(DH)

## JH: Julian history - Julian date of capture for each occasion
JH <- data %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = Julian,
              id_cols = c(id, treatment)) %>%
  rename(TR = treatment) %>%
  mutate(TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(JH)


J <- as.matrix(JH[,which(colnames(JH) == "oc1"):ncol(JH)]) # make matrix of JH

### mcmc setting ---------------------------------------------------------------

n.ad <- 100
n.iter <- 2E+4
n.thin <- max(3, ceiling(n.iter/500))
burn <- ceiling(max(10, n.iter/2))
Sample <- ceiling(n.iter/n.thin)

# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

### bayesian inference ---------------------------------------------------------

## Data for JAGS
### Response
X <- as.matrix(DH[,which(colnames(DH) == "oc1"):ncol(DH)] - 5)
Y <- z <- as.matrix(CH[,which(colnames(CH) == "oc1"):ncol(CH)])
z[Y == 0] <- NA # Y data replaced 0 to NA
ObsF <- apply(Y, 1, function(x) min(which(is.na(x) == 0))) # first observed occasion 
ObsL <- apply(Y, 1, function(x) max(which(x == 1))) # last observed occasion 

# put 1 first to last observed occasion 
for(i in 1:length(ObsF)){
  z[i,ObsF[i]:ObsL[i]] <- 1
}

### Explanatory
Nday <- diff(apply(J, 2, median, na.rm = T)) # interval between occasions
Grid <- count %>% dplyr::select(TR) # Create unique combinations
Grid$G <- 1:length(Grid$TR) # make new no. column
Grid$id <- as.numeric(paste0(Grid$TR)) # make new id column
# marge Grid and Capture history data (Grouping)
groupdat <- tibble(id = as.numeric(paste0(CH$TR))) %>%
  left_join(Grid, by = "id")
G <- groupdat$G

# make input data
Djags <- list(X = X,
              Y = Y,
              Nday = Nday,
              z = z,
              Nind = nrow(Y), # Number of sample
              Nt = ncol(Y), # Number of occasion
              Ng = length(unique(G)), # Number of unique combinations for growth and treatment
              G = G, # grouping data
              ObsF = ObsF
)
# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

# Repeat an equation 3 times
inits <- replicate(3,
                   list(logit.p = matrix(0, nrow = Djags$Ng, ncol = Djags$Nt-1),
                        .RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA ),
                   simplify = F)

for(k in 1:3) inits[[k]]$.RNG.seed <- k

### mcmc using runjags ---------------------------------------------------------
# Parameter estimation by MCMC using runjags
post <- run.jags(m$model,
                 monitor = para,
                 data = Djags,
                 n.chains = 3,
                 inits = inits,
                 method = "parallel",
                 burnin = burn,
                 sample = Sample,
                 adapt = n.ad,
                 thin = n.thin,
                 n.sims = 3,
                 modules = "glm")
# out result
bpost <- MCMCvis::MCMCsummary(post$mcmc,probs = c(0.05,0.5,0.95))

# evaluation of posterior distribution
print(max(na.omit(bpost[, "Rhat"])))

# trace plot
x1 <- data.frame(post$mcmc[[1]])$phi.1.8
x2 <- data.frame(post$mcmc[[2]])$phi.1.8
x3 <- data.frame(post$mcmc[[3]])$phi.1.8

par()
plot(x1, type="l", col = 'red', ylim=c(0,1))
par(new=T)
plot(x2, type="l", col = 'green', ylim=c(0,1), ann=F)
par(new=T)
plot(x3, type="l", col = 'blue', ylim=c(0,1), ann=F)


### output ---------------------------------------------------------------------

bpost <- bpost %>% 
  as_tibble(rownames = "param") %>% 
  rename(lower = '5%',
         median = '50%',
         upper = '95%')

## convert date
date <- apply(J, 2, median, na.rm = T) %>% 
  as.Date(origin = as.Date("1970-01-01"))

dat_date <- tibble(occasion = 1:length(date),
                   date = date)
dat_date

# make new data
# phi age0 May,2016 -> age0 Oct,2016
result <- bpost %>% 
  filter(str_detect(param, "phi")) %>% 
  mutate(x = str_extract(param, pattern = "\\[.{1,}\\]"),
         y = str_remove_all(x, pattern = "\\[|\\]")) %>%
  separate(y, into = c("group", "occasion"), sep = ",") %>% 
  mutate(group = as.numeric(group),
         occasion = as.numeric(occasion)) %>% 
  left_join(Grid, by = c("group" = "G")) %>% 
  left_join(dat_date, by = "occasion") %>% 
  filter(occasion == length(dat_date$date)) %>% 
  mutate(treatment = case_when(TR == 1 ~ "control",
                               TR == 2 ~ "early",
                               TR == 3 ~ "late"),
         treatment = fct_relevel(treatment, c("early", "late"))) %>% 
  left_join(count, by = c("TR")) %>% 
  mutate(n = case_when(treatment == "early" ~ paste0('n = ', as.character(n)),
                       TRUE ~ as.character(n)))

# save to csv
write.csv(result, "output/phi_age0_1_2016cohort_no_growth.csv")

# phi age0 may -> age0 oct
result2 <- bpost %>% 
  filter(str_detect(param, "phi")) %>% 
  mutate(x = str_extract(param, pattern = "\\[.{1,}\\]"),
         y = str_remove_all(x, pattern = "\\[|\\]")) %>%
  separate(y, into = c("group", "occasion"), sep = ",") %>% 
  mutate(group = as.numeric(group),
         occasion = as.numeric(occasion)) %>% 
  left_join(Grid, by = c("group" = "G")) %>% 
  left_join(dat_date, by = "occasion") %>% 
  filter(occasion == 4) %>% 
  mutate(treatment = case_when(TR == 1 ~ "control",
                               TR == 2 ~ "early",
                               TR == 3 ~ "late"),
         treatment = fct_relevel(treatment, c("early", "late"))) %>% 
  left_join(count, by = c("TR")) %>% 
  mutate(n = case_when(treatment == "early" ~ paste0('n = ', as.character(n)),
                       TRUE ~ as.character(n)))

  # save to csv
write.csv(result2, "output/phi_age0_2016cohort_no_growth.csv")

## 2-3. plot -------------------------------------------------------------------

# phi age0 may -> age1 oct
g <- ggplot(result) +
  geom_errorbar(aes(x = treatment, y = mean, ymin = lower, ymax = upper), 
                width = 0, linewidth = 2) +
  geom_point(aes(x = treatment, y = mean), size = 6) +
  ylab("Survival probability") +
  xlab("Treatment") +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0,1.5,0.5)) +
  geom_text(aes(x = treatment, y = 1.1, label = n), size = 5) +
  theme(strip.text=element_text(size=25, family="sans", face="bold"),
        axis.title=element_text(size=28, family="sans", face="bold"),
        axis.text=element_text(size=25, family="sans"),
        strip.background = element_rect(color="black", fill="white", linewidth=0.6, linetype="solid"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = 'none')
g

# save png format
ggsave(g,
       file = "output/figure_phi_age0_1_2016cohort_no_growth.png",
       width = 8,
       height = 6)

# phi age0 Oct,2016 -> age1 Oct,2017
g <- ggplot(result2) +
  geom_errorbar(aes(x = treatment, y = mean, ymin = lower, ymax = upper), 
                width = 0, linewidth = 2) +
  geom_point(aes(x = treatment, y = mean), size = 6) +
  ylab("Survival probability") +
  xlab("Treatment") +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0,1.5,0.5)) +
  geom_text(aes(x = treatment, y = 1.1, label = n), size = 5) +
  theme(strip.text=element_text(size=25, family="sans", face="bold"),
        axis.title=element_text(size=28, family="sans", face="bold"),
        axis.text=element_text(size=25, family="sans"),
        strip.background = element_rect(color="black", fill="white", linewidth=0.6, linetype="solid"),
        panel.background = element_rect("white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = 'none')
g

# save png format
ggsave(g,
       file = "output/figure_phi_age0_2016cohort_no_growth.png",
       width = 8,
       height = 6)




# end --------------------------------------------------------------------------