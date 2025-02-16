
# checked by Terui on 2024.07.22
# Modified by Noda on 2025.01.17

# v2: we use June, Aug. and Oct data in 2016 - 2018

# 1. setup ---------------------------------------------------------------------

# package
pacman::p_load(tidyverse, runjags, loo)
library(RColorBrewer)

# reset environment
rm(list=ls(all.names=T))

# road a model
m <- read.jagsfile("code/model_cjs_r_ver3.R")

# 2. prepare data --------------------------------------------------------------

# read raw data
d0 <- read.csv("./output/rawdata_filtered.csv") 
gc <- read.csv("./output/k4_data.csv")  

# edit data
data <- d0 %>% 
  left_join(gc, by="id") %>% 
  mutate(growth_cluster = as.numeric(ocluster_rename),
         distance = case_when(between(year, 2016, 2019) ~ (new_landmark - 400) * 10,
                              between(year, 2020, 2024) ~ new_landmark * 10),
         age_at_mature2 = case_when(age_at_mature == "0" ~ "0",
                                    age_at_mature == "1" ~ "1",
                                    age_at_mature == "2" ~ "2",
                                    age_at_mature == "3" ~ "2",
                                    age_at_mature == ">1" ~ "2",
                                    age_at_mature == ">2" ~ "2",
                                    age_at_mature == ">3" ~ "2",
                                    TRUE ~ "nd"),
         month2 = case_when(month == 7 ~ 8,
                            month == 9 ~ 10,
                            month == 11 ~ 10,
                            TRUE ~ month),
         occasion = format(as.Date(paste(year, month2, day, sep = "-")), format="%Y-%m")) %>% 
  dplyr::select(id, year, month, ym, occasion, Julian, newsection, treatment, new_landmark, 
         cohort, age_at_mature, age_at_mature2, age2_mature, growth_cluster, distance) %>% 
  filter(distance > 0,
         treatment == "early" | treatment == "late" | treatment == "control",
         age_at_mature2 != "0",
         between(cohort, 2016, 2021)) %>%
  mutate(SecUL = ceiling(distance * 0.1) * 10, # upstream landmark
         SecDL = SecUL - 10, # downstream landmark
         Occasion = as.numeric(factor(occasion)), # add occasion column
         growth = case_when(growth_cluster == 1 ~ 1,
                            growth_cluster == 2 ~ 2,
                            growth_cluster == 3 ~ 3,
                            growth_cluster == 4 ~ 4)) %>%
  filter(id %in% names(which(table(id) > 1)) | ym != '2024-10') %>% # remove individuals only captured on last occasion
  filter(id %in% names(which(table(unlist(tapply(id, treatment, unique))) == 1))) %>%  # remove individuals captured in other treatment
  distinct(id,year,month, .keep_all = TRUE)

# check cohort error ids
eid <- data %>% 
  filter(year < cohort) %>% 
  pull(id)

# check y-m
occasion_list <- tibble(
  ym = sort(unique(data$occasion)),
  oc = sort(unique((data$Occasion)))
)
occasion_list

# remove duplicates for id and occasion
# check
dups <- data %>% 
  group_by(id, occasion) %>%  filter(n() > 1)
data2 <- data %>% 
  distinct(id, occasion, .keep_all = TRUE)


# 3. run CJS model -------------------------------------------------------------
## 3-1 Age 0 -> 1 --------------------------------------------------------------
### 3-1-1 GC and CO ------------------------------------------------------------

out_name <- "output/phi_age0_1_growth_v2.csv"

dat1 <- data2 %>% 
  filter(!is.na(growth_cluster))

## create history matrices
## - Variable transformation
## - `Occasion` - `year`
## - `CL` - `growth` 
## - `CO` - `cohort` 
## Note: An error will occur if there is data with duplicate ID and captured date.

## CH: Capture history - binary information of recaptured or not
CH <- dat1 %>%
  mutate(capture = 1) %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = capture,
              id_cols = c(id, cohort, growth),
              values_fill = list(capture = 0)) %>%
  rename(CO = cohort,
         CL = growth) %>%
  mutate(CO = as.numeric(factor(CO))) %>%
  arrange(id) # arrange by id

# check number of samples
count <- count(CH, CO, CL)
count

# select sorted occasion data
ch <- CH %>%
  dplyr::select(as.factor(sort(unique(as.numeric(dat1$Occasion)))))

# put NA until the first capture
for(i in 1:nrow(ch)) {
  x <- min(which(ch[i,] == 1)) - 1
  if(x != 0) ch[i, 1:x] <- NA
}

CH[, which(colnames(CH) %in% sort(unique(dat1$Occasion)))] <- ch

print(CH)

## DH: Location history - distance from the downstream end to the upstream landmark of a (re)capture subsection
DH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = SecUL,
              id_cols = c(id, cohort, growth)) %>%
  rename(CO = cohort,
         CL = growth) %>%
  mutate(CO = as.numeric(factor(CO))) %>%
  arrange(id)

print(DH)

## JH: Julian history - Julian date of capture for each occasion
JH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = Julian,
              id_cols = c(id, cohort, growth)) %>%
  rename(CO = cohort,
         CL = growth) %>%
  mutate(CO = as.numeric(factor(CO))) %>%
  arrange(id)

print(JH)


J <- as.matrix(JH[,which(colnames(JH) == "1"):ncol(JH)]) # make matrix of JH

#### mcmc setting --------------------------------------------------------------

n.ad <- 100
n.iter <- 2E+4
n.thin <- max(3, ceiling(n.iter/500))
burn <- ceiling(max(10, n.iter/2))
Sample <- ceiling(n.iter/n.thin)

# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

#### bayesian inference --------------------------------------------------------

## Data for JAGS
### Response
X <- as.matrix(DH[,which(colnames(DH) == "1"):ncol(DH)] - 5)
Y <- z <- as.matrix(CH[,which(colnames(CH) == "1"):ncol(CH)])
z[Y == 0] <- NA # Y data replaced 0 to NA
ObsF <- apply(Y, 1, function(x) min(which(is.na(x) == 0))) # first observed occasion 
ObsL <- apply(Y, 1, function(x) max(which(x == 1))) # last observed occasion 

# put 1 first to last observed occasion 
for(i in 1:length(ObsF)){
  z[i,ObsF[i]:ObsL[i]] <- 1
}

### Explanatory
no_coh <- 1:max(CH$CO)
no_cl <- 1:max(CH$CL)
no_g <- 1:(max(CH$CO)*max(CH$CL))

Nday <- diff(apply(J, 2, median, na.rm = T)) # interval between occasions
Grid <- expand.grid(CO = no_coh, CL = no_cl) # Create unique combinations
Grid$G <- no_g # make new no. column
Grid$id <- as.numeric(paste0(Grid$CO,Grid$CL) ) # make new id column
Grid

# marge Grid and Capture history data (Grouping)
groupdat <- tibble(id = as.numeric(paste0(CH$CO, CH$CL))) %>%
  left_join(Grid, by = "id")

G <- groupdat$G
first_occasion <-rep(c(3,7,10,11,12,13), times = 4)

# make input data
Djags <- list(X = X,
              Y = Y,
              Nday = Nday,
              z = z,
              Nind = nrow(Y), # Number of sample
              Nt = ncol(Y), # Number of occasion
              Ng = length(unique(G)), # Number of unique combinations for growth and treatment
              G = G, # grouping data
              Fo = first_occasion,
              ObsF = ObsF
)


# Repeat an equation 3 times
inits <- replicate(3,
                   list(logit.p = matrix(0, nrow = Djags$Ng, ncol = Djags$Nt-1),
                        .RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA ),
                   simplify = F)

for(k in 1:3) inits[[k]]$.RNG.seed <- k

#### mcmc using runjags --------------------------------------------------------
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
x1 <- data.frame(post$mcmc[[1]])$phi.1.2
x2 <- data.frame(post$mcmc[[2]])$phi.1.2
x3 <- data.frame(post$mcmc[[3]])$phi.1.2

par()
plot(x1, type="l", col = 'red', ylim=c(0,1))
par(new=T)
plot(x2, type="l", col = 'green', ylim=c(0,1), ann=F)
par(new=T)
plot(x3, type="l", col = 'blue', ylim=c(0,1), ann=F)


#### output --------------------------------------------------------------------

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
result <- bpost %>% 
  filter(str_detect(param, "phi")) %>% 
  mutate(x = str_extract(param, pattern = "\\[.{1,}\\]"),
         y = str_remove_all(x, pattern = "\\[|\\]")) %>%
  separate(y, into = c("group", "occasion"), sep = ",") %>% 
  mutate(group = as.numeric(group),
         occasion = as.numeric(occasion)) %>% 
  left_join(Grid, by = c("group" = "G")) %>% 
  left_join(dat_date, by = "occasion") %>% 
  mutate(cohort = case_when(CO == 1 ~ "2016",
                            CO == 2 ~ "2017",
                            CO == 3 ~ "2018",
                            CO == 4 ~ "2019",
                            CO == 5 ~ "2020",
                            CO == 6 ~ "2021"),
         year = year(date),
         month = month(date),
         age = year - as.numeric(cohort)) %>% 
  left_join(count, by = c("CO","CL")) %>% 
  mutate(n = case_when(cohort == "2016" & CL == 1 ~ paste0('n = ', as.character(n)),
                       TRUE ~ as.character(n)))

# save to csv
write.csv(result, out_name)

#### make figure ---------------------------------------------------------------

# filtering
dat_sub <- result %>% 
  filter(age == 1,
         month >= 10 | (cohort == 2017 & month == 7)) %>% 
  mutate(growth = as.character(CL))

# color list
color_list = brewer.pal(4, "RdYlBu")

g <- ggplot(dat_sub) +
  geom_errorbar(aes(x = growth, y = mean, ymin = lower, ymax = upper, colour = growth), 
                width = 0, linewidth = 1) +
  geom_point(aes(x = growth, y = mean, colour = growth), size = 4) +
  facet_grid(~cohort) +
  ylab("Survival probability") +
  xlab("Glowth cluster") +
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
       file = paste("output/figure_phi_age0_1_growth_cohort_v2.png", sep = ""),
       dpi = 200, width = 14, height = 4)

### 3-2 GC, CO and TR ----------------------------------------------------------

out_name <- "output/phi_age0_1_growth_treatment_v2.csv"

# filtering by growth cluster
dat1 <- data2 %>% 
  filter(!is.na(growth)) %>% 
  mutate(treatment2 = case_when(cohort >= 2019 ~ "control",
                                TRUE ~ treatment))

## create history matrices
## - Variable transformation
## - `Occasion` - `year`
## - `CL` - `growth` 
## - `CO` - `cohort` 
## - `TR` - `treatment` 
## Note: An error will occur if there is data with duplicate ID and captured date.

## CH: Capture history - binary information of recaptured or not
CH <- dat1 %>%
  mutate(capture = 1) %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = capture,
              id_cols = c(id, cohort, growth, treatment2),
              values_fill = list(capture = 0)) %>%
  rename(CO = cohort,
         CL = growth,
         TR = treatment2) %>%
  mutate(CO = as.numeric(factor(CO)),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id) # arrange by id

# check number of samples
count <- count(CH, CO, TR, CL)
count

# select sorted occasion data
ch <- CH %>%
  dplyr::select(as.factor(sort(unique(as.numeric(dat1$Occasion)))))

# put NA until the first capture
for(i in 1:nrow(ch)) {
  x <- min(which(ch[i,] == 1)) - 1
  if(x != 0) ch[i, 1:x] <- NA
}

CH[, which(colnames(CH) %in% sort(unique(dat1$Occasion)))] <- ch

print(CH)

## DH: Location history - distance from the downstream end to the upstream landmark of a (re)capture subsection
DH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = SecUL,
              id_cols = c(id, cohort, growth, treatment2)) %>%
  rename(CO = cohort,
         CL = growth,
         TR = treatment2) %>%
  mutate(CO = as.numeric(factor(CO)),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(DH)

## JH: Julian history - Julian date of capture for each occasion
JH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = Julian,
              id_cols = c(id, cohort, growth, treatment2)) %>%
  rename(CO = cohort,
         CL = growth,
         TR = treatment2) %>%
  mutate(CO = as.numeric(factor(CO)),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(JH)

J <- as.matrix(JH[,which(colnames(JH) == "1"):ncol(JH)]) # make matrix of JH

#### mcmc setting --------------------------------------------------------------

n.ad <- 100
n.iter <- 2E+4
n.thin <- max(3, ceiling(n.iter/500))
burn <- ceiling(max(10, n.iter/2))
Sample <- ceiling(n.iter/n.thin)

# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

#### bayesian inference --------------------------------------------------------

## Data for JAGS
### Response
X <- as.matrix(DH[,which(colnames(DH) == "1"):ncol(DH)] - 5)
Y <- z <- as.matrix(CH[,which(colnames(CH) == "1"):ncol(CH)])
z[Y == 0] <- NA # Y data replaced 0 to NA
ObsF <- apply(Y, 1, function(x) min(which(is.na(x) == 0))) # first observed occasion 
ObsL <- apply(Y, 1, function(x) max(which(x == 1))) # last observed occasion 

# put 1 first to last observed occasion 
for(i in 1:length(ObsF)){
  z[i,ObsF[i]:ObsL[i]] <- 1
}

### Explanatory
Nday <- diff(apply(J, 2, median, na.rm = T)) # interval between occasions
Grid <- count %>% dplyr::select(CO,CL,TR) # Create unique combinations
Grid$G <- 1:length(Grid$CO) # make new no. column
Grid$id <- as.numeric(paste0(Grid$CO,Grid$CL,Grid$TR)) # make new id column

# marge Grid and Capture history data (Grouping)
groupdat <- tibble(id = as.numeric(paste0(CH$CO, CH$CL, CH$TR))) %>%
  left_join(Grid, by = "id")

G <- groupdat$G
first_occasion <- Grid %>% 
  mutate(Fo = case_when(CO == 1 ~ 3,
                        CO == 2 ~ 7,
                        CO == 3 ~ 10,
                        CO == 4 ~ 11,
                        CO == 5 ~ 12,
                        CO == 6 ~ 13))

# make input data
Djags <- list(X = X,
              Y = Y,
              Nday = Nday,
              z = z,
              Nind = nrow(Y), # Number of sample
              Nt = ncol(Y), # Number of occasion
              Ng = length(unique(G)), # Number of unique combinations for groups
              G = G, # grouping data
              Fo = first_occasion$Fo,
              ObsF = ObsF
)


# Repeat an equation 3 times
inits <- replicate(3,
                   list(logit.p = matrix(0, nrow = Djags$Ng, ncol = Djags$Nt-1),
                        .RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA ),
                   simplify = F)

for(k in 1:3) inits[[k]]$.RNG.seed <- k

#### mcmc using runjags --------------------------------------------------------
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

#### output --------------------------------------------------------------------

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
result <- bpost %>% 
  filter(str_detect(param, "phi")) %>% 
  mutate(x = str_extract(param, pattern = "\\[.{1,}\\]"),
         y = str_remove_all(x, pattern = "\\[|\\]")) %>%
  separate(y, into = c("group", "occasion"), sep = ",") %>% 
  mutate(group = as.numeric(group),
         occasion = as.numeric(occasion)) %>% 
  left_join(Grid, by = c("group" = "G")) %>% 
  left_join(dat_date, by = "occasion") %>% 
  mutate(cohort = case_when(CO == 1 ~ "2016",
                            CO == 2 ~ "2017",
                            CO == 3 ~ "2018",
                            CO == 4 ~ "2019",
                            CO == 5 ~ "2020",
                            CO == 6 ~ "2021"),
         treatment = case_when(TR == 1 ~ "control",
                               TR == 2 ~ "early",
                               TR == 3 ~ "late"),
         treatment = fct_relevel(treatment, c("early", "late")),
         year = year(date),
         month = month(date),
         age = year - as.numeric(cohort)) %>% 
  left_join(count, by = c("CO","CL","TR")) %>% 
  mutate(n = case_when(cohort == "2016" & CL == 1 & treatment == "early" ~ paste0('n = ', as.character(n)),
                       TRUE ~ as.character(n)))

# save to csv
write.csv(result, out_name)

#### make figure ---------------------------------------------------------------

## 2016 ~ 2021 cohort

# filtering
dat_sub <- result %>% 
  filter(age == 1,
         month >= 10 | (cohort == 2017 & month == 7)) %>% 
  mutate(growth = as.character(CL))

# color list
color_list = brewer.pal(4, "RdYlBu")

g <- ggplot(dat_sub) +
  geom_errorbar(aes(x = growth, y = mean, ymin = lower, ymax = upper, colour = growth), 
                width = 0, linewidth = 1) +
  geom_point(aes(x = growth, y = mean, colour = growth), size = 4) +
  facet_grid(treatment~cohort) +
  ylab("Survival probability") +
  xlab("Glowth cluster") +
  theme_bw() +
  scale_color_manual('glowth cluster: ', values = color_list) +
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0,1.5,0.5)) +
  geom_text(aes(x = CL, y = 1.1, label = n), size = 4) +
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
       file = paste("output/figure_phi_age0_1_growth_cohort_treatment_v2.png", sep = ""),
       dpi = 200, width = 14, height = 8)

## 2016 cohort
# filtering
dat_sub <- result %>% 
  filter(age == 1,
         month == 10,
         cohort == 2016) %>% 
  mutate(growth = as.character(CL))

# color list

g <- ggplot(dat_sub) +
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
       file = "output/figure_phi_age0_1_oct_2016cohort_ver2.png",
       width = 10,
       height = 4)

### 3-3 CO and TR --------------------------------------------------------------

out_name <- "output/phi_age0_1_treatment_v2.csv"

dat1 <- data2 %>% 
  mutate(treatment2 = case_when(cohort >= 2019 ~ "control",
                                TRUE ~ treatment))

## create history matrices
## - Variable transformation
## - `Occasion` - `year`
## - `CO` - `cohort` 
## - `TR` - `treatment` 
## Note: An error will occur if there is data with duplicate ID and captured date.

## CH: Capture history - binary information of recaptured or not
CH <- dat1 %>%
  mutate(capture = 1) %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = capture,
              id_cols = c(id, cohort, treatment2),
              values_fill = list(capture = 0)) %>%
  rename(CO = cohort,
         TR = treatment2) %>%
  mutate(CO = as.numeric(factor(CO)),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id) # arrange by id

# check number of samples
count <- count(CH, CO, TR)
count

# select sorted occasion data
ch <- CH %>%
  dplyr::select(as.factor(sort(unique(as.numeric(dat1$Occasion)))))

# put NA until the first capture
for(i in 1:nrow(ch)) {
  x <- min(which(ch[i,] == 1)) - 1
  if(x != 0) ch[i, 1:x] <- NA
}

CH[, which(colnames(CH) %in% sort(unique(dat1$Occasion)))] <- ch

print(CH)

## DH: Location history - distance from the downstream end to the upstream landmark of a (re)capture subsection
DH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = SecUL,
              id_cols = c(id, cohort, treatment2)) %>%
  rename(CO = cohort,
         TR = treatment2) %>%
  mutate(CO = as.numeric(factor(CO)),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(DH)

## JH: Julian history - Julian date of capture for each occasion
JH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = Julian,
              id_cols = c(id, cohort, treatment2)) %>%
  rename(CO = cohort,
         TR = treatment2) %>%
  mutate(CO = as.numeric(factor(CO)),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(JH)

J <- as.matrix(JH[,which(colnames(JH) == "1"):ncol(JH)]) # make matrix of JH

#### mcmc setting --------------------------------------------------------------

n.ad <- 100
n.iter <- 2E+4
n.thin <- max(3, ceiling(n.iter/500))
burn <- ceiling(max(10, n.iter/2))
Sample <- ceiling(n.iter/n.thin)

# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

#### bayesian inference --------------------------------------------------------

## Data for JAGS
### Response
X <- as.matrix(DH[,which(colnames(DH) == "1"):ncol(DH)] - 5)
Y <- z <- as.matrix(CH[,which(colnames(CH) == "1"):ncol(CH)])
z[Y == 0] <- NA # Y data replaced 0 to NA
ObsF <- apply(Y, 1, function(x) min(which(is.na(x) == 0))) # first observed occasion 
ObsL <- apply(Y, 1, function(x) max(which(x == 1))) # last observed occasion 

# put 1 first to last observed occasion 
for(i in 1:length(ObsF)){
  z[i,ObsF[i]:ObsL[i]] <- 1
}

### Explanatory
Nday <- diff(apply(J, 2, median, na.rm = T)) # interval between occasions
Grid <- count %>% dplyr::select(CO,TR) # Create unique combinations
Grid$G <- 1:length(Grid$CO) # make new no. column
Grid$id <- as.numeric(paste0(Grid$CO,Grid$TR)) # make new id column

# marge Grid and Capture history data (Grouping)
groupdat <- tibble(id = as.numeric(paste0(CH$CO, CH$TR))) %>%
  left_join(Grid, by = "id")

G <- groupdat$G
first_occasion <- Grid %>% 
  mutate(Fo = case_when(CO == 1 ~ 3,
                        CO == 2 ~ 7,
                        CO == 3 ~ 10,
                        CO == 4 ~ 11,
                        CO == 5 ~ 12,
                        CO == 6 ~ 13))

# make input data
Djags <- list(X = X,
              Y = Y,
              Nday = Nday,
              z = z,
              Nind = nrow(Y), # Number of sample
              Nt = ncol(Y), # Number of occasion
              Ng = length(unique(G)), # Number of unique combinations for groups
              G = G, # grouping data
              Fo = first_occasion$Fo,
              ObsF = ObsF
)

# Repeat an equation 3 times
inits <- replicate(3,
                   list(logit.p = matrix(0, nrow = Djags$Ng, ncol = Djags$Nt-1),
                        .RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA ),
                   simplify = F)

for(k in 1:3) inits[[k]]$.RNG.seed <- k

#### mcmc using runjags --------------------------------------------------------
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

#### output --------------------------------------------------------------------

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
result <- bpost %>% 
  filter(str_detect(param, "phi")) %>% 
  mutate(x = str_extract(param, pattern = "\\[.{1,}\\]"),
         y = str_remove_all(x, pattern = "\\[|\\]")) %>%
  separate(y, into = c("group", "occasion"), sep = ",") %>% 
  mutate(group = as.numeric(group),
         occasion = as.numeric(occasion)) %>% 
  left_join(Grid, by = c("group" = "G")) %>% 
  left_join(dat_date, by = "occasion") %>% 
  mutate(cohort = case_when(CO == 1 ~ "2016",
                            CO == 2 ~ "2017",
                            CO == 3 ~ "2018",
                            CO == 4 ~ "2019",
                            CO == 5 ~ "2020",
                            CO == 6 ~ "2021"),
         treatment = case_when(TR == 1 ~ "control",
                               TR == 2 ~ "early",
                               TR == 3 ~ "late"),
         treatment = fct_relevel(treatment, c("early", "late")),
         year = year(date),
         month = month(date),
         age = year - as.numeric(cohort)) %>% 
  left_join(count, by = c("CO","TR")) %>% 
  mutate(n = case_when(cohort == "2016" & treatment == "early" ~ paste0('n = ', as.character(n)),
                       TRUE ~ as.character(n)))

# save to csv
write.csv(result, out_name)

#### make figure ---------------------------------------------------------------

## 2016 ~ 2021 cohort

# filtering
dat_sub <- result %>% 
  filter(age == 1,
         month >= 10 | (cohort == 2017 & month == 7))

g <- ggplot(dat_sub) +
  geom_errorbar(aes(x = treatment, y = mean, ymin = lower, ymax = upper), 
                width = 0, linewidth = 1) +
  geom_point(aes(x = treatment, y = mean), size = 4) +
  facet_grid(~cohort) +
  ylab("Survival probability") +
  xlab("Treatment") +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0,1.5,0.5)) +
  geom_text(aes(x = treatment, y = 1.1, label = n), size = 4) +
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
       file = paste("output/figure_phi_age0_1_cohort_treatment_v2.png", sep = ""),
       dpi = 200, width = 20, height = 4)

# ...----

## 3-2 Age 1 -> 2 --------------------------------------------------------------
### 3-2-1 GC, CO, TR and MT ----------------------------------------------------
## Grouped by maturity at age 1 or not
## use annual captured data

out_name <- "output/phi_age1_2_growth_treatment_age1mature_v2.csv"

# filtering by growth cluster
dat1 <- data2 %>% 
  filter(!is.na(growth),
         age_at_mature2 != "nd",
         between(month, 9, 11) | 
           (year == 2018 & between(month, 7, 8))) %>% 
  mutate(treatment2 = case_when(cohort >= 2019 ~ "control",
                                TRUE ~ treatment),
         Occasion = as.numeric(factor(year)))
unique(dat1$cohort)

# check y-m
occasion_list <- tibble(
  ym = sort(unique(dat1$occasion)),
  oc = sort(unique((dat1$Occasion)))
)
occasion_list

## create history matrices
## - Variable transformation
## - `Occasion` - `year`
## - `CL` - `growth` 
## - `CO` - `cohort` 
## - `TR` - `treatment` 
## - `MT` - `age at maturity` 
## Note: An error will occur if there is data with duplicate ID and captured date.

## CH: Capture history - binary information of recaptured or not
CH <- dat1 %>%
  mutate(capture = 1) %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = capture,
              id_cols = c(id, cohort, growth, treatment2, age_at_mature2),
              values_fill = list(capture = 0)) %>%
  rename(CO = cohort,
         CL = growth,
         TR = treatment2,
         MT = age_at_mature2) %>%
  mutate(CO = as.numeric(factor(CO)),
         MT = as.factor(MT),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id) # arrange by id

# check number of samples
count <- count(CH, CO, TR, CL, MT)
count

tmp <- count(CH, CO, TR, CL, MT, `3`) %>% 
  filter(CO == 1)

# select sorted occasion data
ch <- CH %>%
  dplyr::select(as.factor(sort(unique(as.numeric(dat1$Occasion)))))

# put NA until the first capture
for(i in 1:nrow(ch)) {
  x <- min(which(ch[i,] == 1)) - 1
  if(x != 0) ch[i, 1:x] <- NA
}

CH[, which(colnames(CH) %in% sort(unique(dat1$Occasion)))] <- ch

print(CH)

## DH: Location history - distance from the downstream end to the upstream landmark of a (re)capture subsection
DH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = SecUL,
              id_cols = c(id, cohort, growth, treatment2, age_at_mature2)) %>%
  rename(CO = cohort,
         CL = growth,
         TR = treatment2,
         MT = age_at_mature2) %>%
  mutate(CO = as.numeric(factor(CO)),
         MT = as.factor(MT),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(DH)

## JH: Julian history - Julian date of capture for each occasion
JH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = Julian,
              id_cols = c(id, cohort, growth, treatment2, age_at_mature2)) %>%
  rename(CO = cohort,
         CL = growth,
         TR = treatment2,
         MT = age_at_mature2) %>%
  mutate(CO = as.numeric(factor(CO)),
         MT = as.factor(MT),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(JH)

J <- as.matrix(JH[,which(colnames(JH) == "1"):ncol(JH)]) # make matrix of JH

#### mcmc setting --------------------------------------------------------------

n.ad <- 100
n.iter <- 2E+4
n.thin <- max(3, ceiling(n.iter/500))
burn <- ceiling(max(10, n.iter/2))
Sample <- ceiling(n.iter/n.thin)

# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

#### bayesian inference --------------------------------------------------------

## Data for JAGS
### Response
X <- as.matrix(DH[,which(colnames(DH) == "1"):ncol(DH)] - 5)
Y <- z <- as.matrix(CH[,which(colnames(CH) == "1"):ncol(CH)])
z[Y == 0] <- NA # Y data replaced 0 to NA
ObsF <- apply(Y, 1, function(x) min(which(is.na(x) == 0))) # first observed occasion 
ObsL <- apply(Y, 1, function(x) max(which(x == 1))) # last observed occasion 

# put 1 first to last observed occasion 
for(i in 1:length(ObsF)){
  z[i,ObsF[i]:ObsL[i]] <- 1
}

### Explanatory
Nday <- diff(apply(J, 2, median, na.rm = T)) # interval between occasions
Grid <- count %>% dplyr::select(CO,CL,TR,MT) # Create unique combinations
Grid$G <- 1:length(Grid$CO) # make new no. column
Grid$id <- as.numeric(paste0(Grid$CO,Grid$CL,Grid$TR,Grid$MT)) # make new id column

# marge Grid and Capture history data (Grouping)
groupdat <- tibble(id = as.numeric(paste0(CH$CO, CH$CL, CH$TR, CH$MT))) %>%
  left_join(Grid, by = "id")

G <- groupdat$G
first_occasion <- Grid %>% 
  mutate(Fo = case_when(CO == 1 ~ 2,
                        CO == 2 ~ 4,
                        CO == 3 ~ 5,
                        CO == 4 ~ 6,
                        CO == 5 ~ 7))

# make input data
Djags <- list(X = X,
              Y = Y,
              Nday = Nday,
              z = z,
              Nind = nrow(Y), # Number of sample
              Nt = ncol(Y), # Number of occasion
              Ng = length(unique(G)), # Number of unique combinations for groups
              G = G, # grouping data
              Fo = first_occasion$Fo,
              ObsF = ObsF
)


# Repeat an equation 3 times
inits <- replicate(3,
                   list(logit.p = matrix(0, nrow = Djags$Ng, ncol = Djags$Nt-1),
                        .RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA ),
                   simplify = F)

for(k in 1:3) inits[[k]]$.RNG.seed <- k

#### mcmc using runjags --------------------------------------------------------
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

#### output --------------------------------------------------------------------

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
result <- bpost %>% 
  filter(str_detect(param, "phi")) %>% 
  mutate(x = str_extract(param, pattern = "\\[.{1,}\\]"),
         y = str_remove_all(x, pattern = "\\[|\\]")) %>%
  separate(y, into = c("group", "occasion"), sep = ",") %>% 
  mutate(group = as.numeric(group),
         occasion = as.numeric(occasion)) %>% 
  left_join(Grid, by = c("group" = "G")) %>% 
  left_join(dat_date, by = "occasion") %>% 
  mutate(cohort = case_when(CO == 1 ~ "2016",
                            CO == 2 ~ "2018",
                            CO == 3 ~ "2019",
                            CO == 4 ~ "2020",
                            CO == 5 ~ "2021"),
         treatment = case_when(TR == 1 ~ "control",
                               TR == 2 ~ "early",
                               TR == 3 ~ "late"),
         treatment = fct_relevel(treatment, c("early", "late")),
         age_at_maturity = MT,
         year = year(date),
         month = month(date),
         age = year - as.numeric(cohort)) %>% 
  left_join(count, by = c("CO","CL","TR","MT")) %>% 
  mutate(n = case_when(cohort == "2016" & CL == 1 & treatment == "early" ~ paste0('n = ', as.character(n)),
                       TRUE ~ as.character(n)))

# save to csv
write.csv(result, out_name)

#### make figure ---------------------------------------------------------------

result <- read_csv(out_name)
# filtering
dat_sub <- result %>% 
  filter(age == 2,
         month >= 7) %>% 
  mutate(treatment = fct_relevel(treatment, c("early", "late")),
         growth = as.character(CL),
         age_at_maturity = as.character(age_at_maturity))

# color list
color_list = brewer.pal(4, "RdYlBu")

## 2016 ~ 2021 cohort
g <- ggplot(dat_sub, aes(x = growth, y = mean, ymin = lower, ymax = upper, 
                         colour = growth, group = age_at_maturity, shape = age_at_maturity)) +
  geom_errorbar(width = 0, linewidth = 1, position = position_dodge(.7)) +
  geom_point(size = 4, position = position_dodge(.7)) +
  facet_grid(treatment~cohort) +
  ylab("Survival probability") +
  xlab("Glowth cluster") +
  theme_bw() +
  scale_color_manual('glowth cluster: ', values = color_list) +
  scale_shape_manual('age at maturity: ', values = c(16,17)) +
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0,1.5,0.5)) +
  geom_text(aes(x = growth, y = 1.1, label = n), size = 4, color = "black",
            position = position_dodge(.7)) +
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
       file = paste("output/figure_phi_age1_2_growth_cohort_treatment_age1mature_v2.png", sep = ""),
       dpi = 200, width = 18, height = 10)


## 2016 cohort
# filtering
dat_sub2 <- dat_sub %>% 
  filter(cohort == 2016)

g <- ggplot(dat_sub2, aes(x = growth, y = mean, ymin = lower, ymax = upper, 
                         colour = growth, group = age_at_maturity, shape = age_at_maturity)) +
  geom_errorbar(width = 0, linewidth = 1, position = position_dodge(.7)) +
  geom_point(size = 4, position = position_dodge(.7)) +
  facet_wrap(facets = ~ treatment, ncol = 3) +
  ylab("Survival probability") +
  xlab("Glowth cluster") +
  theme_bw() +
  scale_color_manual('glowth cluster: ', values = color_list) +
  scale_shape_manual('age at maturity: ', values = c(16,17)) +
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0,1.5,0.5)) +
  geom_text(aes(x = growth, y = 1.1, label = n), size = 4, color = "black",
            position = position_dodge(.7)) +
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
       file = paste("output/figure_phi_age1_2_growth_2016cohort_treatment_age1mature_v2.png", sep = ""),
       dpi = 200, width = 12, height = 5)



### 3-2-2 CO, TR and MT --------------------------------------------------------
## Grouped by maturity at age 1 or not
## use annual captured data

out_name <- "output/phi_age1_2_treatment_age1mature_v2.csv"

dat1 <- data2 %>% 
  filter(!is.na(growth),
         age_at_mature2 != "nd",
         between(month, 9, 11) | 
           (year == 2018 & between(month, 7, 8))) %>% 
  mutate(treatment2 = case_when(cohort >= 2019 ~ "control",
                                TRUE ~ treatment),
         Occasion = as.numeric(factor(year)))
unique(dat1$cohort)

# check y-m
occasion_list <- tibble(
  ym = sort(unique(dat1$occasion)),
  oc = sort(unique((dat1$Occasion)))
)
occasion_list

## create history matrices
## - Variable transformation
## - `Occasion` - `year`
## - `CO` - `cohort` 
## - `TR` - `treatment` 
## - `MT` - `age at maturity` 
## Note: An error will occur if there is data with duplicate ID and captured date.

## CH: Capture history - binary information of recaptured or not
CH <- dat1 %>%
  mutate(capture = 1) %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = capture,
              id_cols = c(id, cohort, treatment2, age_at_mature2),
              values_fill = list(capture = 0)) %>%
  rename(CO = cohort,
         TR = treatment2,
         MT = age_at_mature2) %>%
  mutate(CO = as.numeric(factor(CO)),
         MT = as.factor(MT),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id) # arrange by id

# check number of samples
count <- count(CH, CO, TR, MT)
count

# select sorted occasion data
ch <- CH %>%
  dplyr::select(as.factor(sort(unique(as.numeric(dat1$Occasion)))))

# put NA until the first capture
for(i in 1:nrow(ch)) {
  x <- min(which(ch[i,] == 1)) - 1
  if(x != 0) ch[i, 1:x] <- NA
}

CH[, which(colnames(CH) %in% sort(unique(dat1$Occasion)))] <- ch

print(CH)

## DH: Location history - distance from the downstream end to the upstream landmark of a (re)capture subsection
DH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = SecUL,
              id_cols = c(id, cohort, treatment2, age_at_mature2)) %>%
  rename(CO = cohort,
         TR = treatment2,
         MT = age_at_mature2) %>%
  mutate(CO = as.numeric(factor(CO)),
         MT = as.factor(MT),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(DH)

## JH: Julian history - Julian date of capture for each occasion
JH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = Julian,
              id_cols = c(id, cohort, treatment2, age_at_mature2)) %>%
  rename(CO = cohort,
         TR = treatment2,
         MT = age_at_mature2) %>%
  mutate(CO = as.numeric(factor(CO)),
         MT = as.factor(MT),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(JH)

J <- as.matrix(JH[,which(colnames(JH) == "1"):ncol(JH)]) # make matrix of JH

#### mcmc setting --------------------------------------------------------------

n.ad <- 100
n.iter <- 2E+4
n.thin <- max(3, ceiling(n.iter/500))
burn <- ceiling(max(10, n.iter/2))
Sample <- ceiling(n.iter/n.thin)

# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

#### bayesian inference --------------------------------------------------------

## Data for JAGS
### Response
X <- as.matrix(DH[,which(colnames(DH) == "1"):ncol(DH)] - 5)
Y <- z <- as.matrix(CH[,which(colnames(CH) == "1"):ncol(CH)])
z[Y == 0] <- NA # Y data replaced 0 to NA
ObsF <- apply(Y, 1, function(x) min(which(is.na(x) == 0))) # first observed occasion 
ObsL <- apply(Y, 1, function(x) max(which(x == 1))) # last observed occasion 

# put 1 first to last observed occasion 
for(i in 1:length(ObsF)){
  z[i,ObsF[i]:ObsL[i]] <- 1
}

### Explanatory
Nday <- diff(apply(J, 2, median, na.rm = T)) # interval between occasions
Grid <- count %>% dplyr::select(CO,TR,MT) # Create unique combinations
Grid$G <- 1:length(Grid$CO) # make new no. column
Grid$id <- as.numeric(paste0(Grid$CO,Grid$TR,Grid$MT)) # make new id column

# marge Grid and Capture history data (Grouping)
groupdat <- tibble(id = as.numeric(paste0(CH$CO, CH$TR, CH$MT))) %>%
  left_join(Grid, by = "id")

G <- groupdat$G
first_occasion <- Grid %>% 
  mutate(Fo = case_when(CO == 1 ~ 2,
                        CO == 2 ~ 4,
                        CO == 3 ~ 5,
                        CO == 4 ~ 6,
                        CO == 5 ~ 7))

# make input data
Djags <- list(X = X,
              Y = Y,
              Nday = Nday,
              z = z,
              Nind = nrow(Y), # Number of sample
              Nt = ncol(Y), # Number of occasion
              Ng = length(unique(G)), # Number of unique combinations for groups
              G = G, # grouping data
              Fo = first_occasion$Fo,
              ObsF = ObsF
)


# Repeat an equation 3 times
inits <- replicate(3,
                   list(logit.p = matrix(0, nrow = Djags$Ng, ncol = Djags$Nt-1),
                        .RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA ),
                   simplify = F)

for(k in 1:3) inits[[k]]$.RNG.seed <- k

#### mcmc using runjags --------------------------------------------------------
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

#### output --------------------------------------------------------------------

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
result <- bpost %>% 
  filter(str_detect(param, "phi")) %>% 
  mutate(x = str_extract(param, pattern = "\\[.{1,}\\]"),
         y = str_remove_all(x, pattern = "\\[|\\]")) %>%
  separate(y, into = c("group", "occasion"), sep = ",") %>% 
  mutate(group = as.numeric(group),
         occasion = as.numeric(occasion)) %>% 
  left_join(Grid, by = c("group" = "G")) %>% 
  left_join(dat_date, by = "occasion") %>% 
  mutate(cohort = case_when(CO == 1 ~ "2016",
                            CO == 2 ~ "2018",
                            CO == 3 ~ "2019",
                            CO == 4 ~ "2020",
                            CO == 5 ~ "2021"),
         treatment = case_when(TR == 1 ~ "control",
                               TR == 2 ~ "early",
                               TR == 3 ~ "late"),
         treatment = fct_relevel(treatment, c("early", "late")),
         age_at_maturity = MT,
         year = year(date),
         month = month(date),
         age = year - as.numeric(cohort)) %>% 
  left_join(count, by = c("CO","TR","MT")) %>% 
  mutate(n = case_when(cohort == "2016" & treatment == "early" & age_at_maturity == 1 ~ paste0('n = ', as.character(n)),
                       TRUE ~ as.character(n)))

# save to csv
write.csv(result, out_name)

#### make figure ---------------------------------------------------------------

## 2016 ~ 2021 cohort

# filtering
dat_sub <- result %>% 
  filter(age == 2,
         month >= 7)

g <- ggplot(dat_sub, aes(x = treatment, y = mean, ymin = lower, ymax = upper, 
                        group = age_at_maturity, shape = age_at_maturity)) +
  geom_errorbar(width = 0, linewidth = 1, position = position_dodge(.7)) +
  geom_point(size = 4, position = position_dodge(.7)) +
  facet_grid(~cohort) +
  ylab("Survival probability") +
  xlab("Treatment") +
  theme_bw() +
  scale_shape_manual('age at maturity: ', values = c(16,17)) +
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0,1.5,0.5)) +
  geom_text(aes(x = treatment, y = 1.1, label = n), size = 4, position = position_dodge(.7)) +
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
       file = paste("output/figure_phi_age1_2_cohort_treatment_age1mature_v2.png", sep = ""),
       dpi = 200, width = 18, height = 5)

### 3-2-3 GC, CO and TR --------------------------------------------------------
## Grouped by maturity at age 1 or not
## use annual captured data

out_name <- "output/phi_age1_2_growth_treatment_v2.csv"

# filtering by growth cluster
dat1 <- data2 %>% 
  filter(!is.na(growth),
         age_at_mature2 != "nd",
         between(month, 9, 11) | 
           (year == 2018 & between(month, 7, 8))) %>% 
  mutate(treatment2 = case_when(cohort >= 2019 ~ "control",
                                TRUE ~ treatment),
         Occasion = as.numeric(factor(year)))
unique(dat1$cohort)

# check y-m
occasion_list <- tibble(
  ym = sort(unique(dat1$occasion)),
  oc = sort(unique((dat1$Occasion)))
)
occasion_list

## create history matrices
## - Variable transformation
## - `Occasion` - `year`
## - `CL` - `growth` 
## - `CO` - `cohort` 
## - `TR` - `treatment` 
## Note: An error will occur if there is data with duplicate ID and captured date.

## CH: Capture history - binary information of recaptured or not
CH <- dat1 %>%
  mutate(capture = 1) %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = capture,
              id_cols = c(id, cohort, growth, treatment2),
              values_fill = list(capture = 0)) %>%
  rename(CO = cohort,
         CL = growth,
         TR = treatment2) %>%
  mutate(CO = as.numeric(factor(CO)),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id) # arrange by id

# check number of samples
count <- count(CH, CO, TR, CL)
count

# select sorted occasion data
ch <- CH %>%
  dplyr::select(as.factor(sort(unique(as.numeric(dat1$Occasion)))))

# put NA until the first capture
for(i in 1:nrow(ch)) {
  x <- min(which(ch[i,] == 1)) - 1
  if(x != 0) ch[i, 1:x] <- NA
}

CH[, which(colnames(CH) %in% sort(unique(dat1$Occasion)))] <- ch

print(CH)

## DH: Location history - distance from the downstream end to the upstream landmark of a (re)capture subsection
DH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = SecUL,
              id_cols = c(id, cohort, growth, treatment2)) %>%
  rename(CO = cohort,
         CL = growth,
         TR = treatment2) %>%
  mutate(CO = as.numeric(factor(CO)),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(DH)

## JH: Julian history - Julian date of capture for each occasion
JH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = Julian,
              id_cols = c(id, cohort, growth, treatment2)) %>%
  rename(CO = cohort,
         CL = growth,
         TR = treatment2) %>%
  mutate(CO = as.numeric(factor(CO)),
         TR = as.numeric(factor(TR,
                                levels = c("control",
                                           "early",
                                           "late")))) %>%
  arrange(id)

print(JH)

J <- as.matrix(JH[,which(colnames(JH) == "1"):ncol(JH)]) # make matrix of JH

#### mcmc setting --------------------------------------------------------------

n.ad <- 100
n.iter <- 2E+4
n.thin <- max(3, ceiling(n.iter/500))
burn <- ceiling(max(10, n.iter/2))
Sample <- ceiling(n.iter/n.thin)

# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

#### bayesian inference --------------------------------------------------------

## Data for JAGS
### Response
X <- as.matrix(DH[,which(colnames(DH) == "1"):ncol(DH)] - 5)
Y <- z <- as.matrix(CH[,which(colnames(CH) == "1"):ncol(CH)])
z[Y == 0] <- NA # Y data replaced 0 to NA
ObsF <- apply(Y, 1, function(x) min(which(is.na(x) == 0))) # first observed occasion 
ObsL <- apply(Y, 1, function(x) max(which(x == 1))) # last observed occasion 

# put 1 first to last observed occasion 
for(i in 1:length(ObsF)){
  z[i,ObsF[i]:ObsL[i]] <- 1
}

### Explanatory
Nday <- diff(apply(J, 2, median, na.rm = T)) # interval between occasions
Grid <- count %>% dplyr::select(CO,CL,TR) # Create unique combinations
Grid$G <- 1:length(Grid$CO) # make new no. column
Grid$id <- as.numeric(paste0(Grid$CO,Grid$CL,Grid$TR)) # make new id column

# marge Grid and Capture history data (Grouping)
groupdat <- tibble(id = as.numeric(paste0(CH$CO, CH$CL, CH$TR))) %>%
  left_join(Grid, by = "id")

G <- groupdat$G
first_occasion <- Grid %>% 
  mutate(Fo = case_when(CO == 1 ~ 2,
                        CO == 2 ~ 4,
                        CO == 3 ~ 5,
                        CO == 4 ~ 6,
                        CO == 5 ~ 7))

# make input data
Djags <- list(X = X,
              Y = Y,
              Nday = Nday,
              z = z,
              Nind = nrow(Y), # Number of sample
              Nt = ncol(Y), # Number of occasion
              Ng = length(unique(G)), # Number of unique combinations for groups
              G = G, # grouping data
              Fo = first_occasion$Fo,
              ObsF = ObsF
)


# Repeat an equation 3 times
inits <- replicate(3,
                   list(logit.p = matrix(0, nrow = Djags$Ng, ncol = Djags$Nt-1),
                        .RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA ),
                   simplify = F)

for(k in 1:3) inits[[k]]$.RNG.seed <- k

#### mcmc using runjags --------------------------------------------------------
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

#### output --------------------------------------------------------------------

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
result <- bpost %>% 
  filter(str_detect(param, "phi")) %>% 
  mutate(x = str_extract(param, pattern = "\\[.{1,}\\]"),
         y = str_remove_all(x, pattern = "\\[|\\]")) %>%
  separate(y, into = c("group", "occasion"), sep = ",") %>% 
  mutate(group = as.numeric(group),
         occasion = as.numeric(occasion)) %>% 
  left_join(Grid, by = c("group" = "G")) %>% 
  left_join(dat_date, by = "occasion") %>% 
  mutate(cohort = case_when(CO == 1 ~ "2016",
                            CO == 2 ~ "2018",
                            CO == 3 ~ "2019",
                            CO == 4 ~ "2020",
                            CO == 5 ~ "2021"),
         treatment = case_when(TR == 1 ~ "control",
                               TR == 2 ~ "early",
                               TR == 3 ~ "late"),
         treatment = fct_relevel(treatment, c("early", "late")),
         year = year(date),
         month = month(date),
         age = year - as.numeric(cohort)) %>% 
  left_join(count, by = c("CO","CL","TR")) %>% 
  mutate(n = case_when(cohort == "2016" & CL == 1 & treatment == "early" ~ paste0('n = ', as.character(n)),
                       TRUE ~ as.character(n)))

# save to csv
write.csv(result, out_name)

#### make figure ---------------------------------------------------------------

## 2016 ~ 2021 cohort

# filtering
dat_sub <- result %>% 
  filter(age == 2,
         month >= 7) %>% 
  mutate(growth = as.character(CL))

# color list
color_list = brewer.pal(4, "RdYlBu")

g <- ggplot(dat_sub) +
  geom_errorbar(aes(x = growth, y = mean, ymin = lower, ymax = upper, colour = growth), 
                width = 0, linewidth = 1) +
  geom_point(aes(x = growth, y = mean, colour = growth), size = 4) +
  facet_grid(treatment~cohort) +
  ylab("Survival probability") +
  xlab("Glowth cluster") +
  theme_bw() +
  scale_color_manual('glowth cluster: ', values = color_list) +
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0,1.5,0.5)) +
  geom_text(aes(x = CL, y = 1.1, label = n), size = 4) +
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
       file = paste("output/figure_phi_age1_2_growth_cohort_treatment_v2.png", sep = ""),
       dpi = 200, width = 14, height = 8)

## 3-3 Age 2 -> 3 --------------------------------------------------------------
### 3-3-1 GC, CO and MT --------------------------------------------------------
## Grouped by maturity at age 1 or not
## use annual captured data
## use 2018 - 2021 cohorts

out_name <- "output/phi_age2_3_growth_age2mature_v2.csv"

# filtering by growth cluster
dat1 <- data2 %>% 
  filter(!is.na(growth),
         age2_mature != "nd",
         cohort >= 2018) %>% 
  mutate(Occasion = as.numeric(factor(year))) # edit occasion

unique(dat1$cohort)

# check y-m
occasion_list <- tibble(
  ym = sort(unique(dat1$occasion)),
  oc = sort(unique((dat1$Occasion)))
)
occasion_list

## create history matrices
## - Variable transformation
## - `Occasion` - `year`
## - `CL` - `growth` 
## - `CO` - `cohort`
## - `MT` - `age at maturity` 
## Note: An error will occur if there is data with duplicate ID and captured date.

## CH: Capture history - binary information of recaptured or not
CH <- dat1 %>%
  mutate(capture = 1) %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = capture,
              id_cols = c(id, cohort, growth, age2_mature),
              values_fill = list(capture = 0)) %>%
  rename(CO = cohort,
         CL = growth,
         MT = age2_mature) %>%
  mutate(CO = as.numeric(factor(CO)),
         MT = as.factor(MT)) %>%
  arrange(id) # arrange by id

# check number of samples
count <- count(CH, CO, CL, MT)
count

# select sorted occasion data
ch <- CH %>%
  dplyr::select(as.factor(sort(unique(as.numeric(dat1$Occasion)))))

# put NA until the first capture
for(i in 1:nrow(ch)) {
  x <- min(which(ch[i,] == 1)) - 1
  if(x != 0) ch[i, 1:x] <- NA
}

CH[, which(colnames(CH) %in% sort(unique(dat1$Occasion)))] <- ch

print(CH)

## DH: Location history - distance from the downstream end to the upstream landmark of a (re)capture subsection
DH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = SecUL,
              id_cols = c(id, cohort, growth, age2_mature)) %>%
  rename(CO = cohort,
         CL = growth,
         MT = age2_mature) %>%
  mutate(CO = as.numeric(factor(CO)),
         MT = as.factor(MT)) %>%
  arrange(id)

print(DH)

## JH: Julian history - Julian date of capture for each occasion
JH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = Julian,
              id_cols = c(id, cohort, growth, age2_mature)) %>%
  rename(CO = cohort,
         CL = growth,
         MT = age2_mature) %>%
  mutate(CO = as.numeric(factor(CO)),
         MT = as.factor(MT)) %>%
  arrange(id)

print(JH)

J <- as.matrix(JH[,which(colnames(JH) == "1"):ncol(JH)]) # make matrix of JH

#### mcmc setting --------------------------------------------------------------

n.ad <- 100
n.iter <- 2E+4
n.thin <- max(3, ceiling(n.iter/500))
burn <- ceiling(max(10, n.iter/2))
Sample <- ceiling(n.iter/n.thin)

# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

#### bayesian inference --------------------------------------------------------

## Data for JAGS
### Response
X <- as.matrix(DH[,which(colnames(DH) == "1"):ncol(DH)] - 5)
Y <- z <- as.matrix(CH[,which(colnames(CH) == "1"):ncol(CH)])
z[Y == 0] <- NA # Y data replaced 0 to NA
ObsF <- apply(Y, 1, function(x) min(which(is.na(x) == 0))) # first observed occasion 
ObsL <- apply(Y, 1, function(x) max(which(x == 1))) # last observed occasion 

# put 1 first to last observed occasion 
for(i in 1:length(ObsF)){
  z[i,ObsF[i]:ObsL[i]] <- 1
}

### Explanatory
Nday <- diff(apply(J, 2, median, na.rm = T)) # interval between occasions
Grid <- count %>% dplyr::select(CO,CL,MT) # Create unique combinations
Grid$G <- 1:length(Grid$CO) # make new no. column
Grid$id <- as.numeric(paste0(Grid$CO,Grid$CL,Grid$MT)) # make new id column

# marge Grid and Capture history data (Grouping)
groupdat <- tibble(id = as.numeric(paste0(CH$CO, CH$CL, CH$MT))) %>%
  left_join(Grid, by = "id")

G <- groupdat$G
first_occasion <- Grid %>% 
  mutate(Fo = case_when(CO == 1 ~ 2,
                        CO == 2 ~ 3,
                        CO == 3 ~ 4,
                        CO == 4 ~ 5))

# make input data
Djags <- list(X = X,
              Y = Y,
              Nday = Nday,
              z = z,
              Nind = nrow(Y), # Number of sample
              Nt = ncol(Y), # Number of occasion
              Ng = length(unique(G)), # Number of unique combinations for groups
              G = G, # grouping data
              Fo = first_occasion$Fo,
              ObsF = ObsF
)


# Repeat an equation 3 times
inits <- replicate(3,
                   list(logit.p = matrix(0, nrow = Djags$Ng, ncol = Djags$Nt-1),
                        .RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA ),
                   simplify = F)

for(k in 1:3) inits[[k]]$.RNG.seed <- k

#### mcmc using runjags --------------------------------------------------------
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

#### output --------------------------------------------------------------------

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
result <- bpost %>% 
  filter(str_detect(param, "phi")) %>% 
  mutate(x = str_extract(param, pattern = "\\[.{1,}\\]"),
         y = str_remove_all(x, pattern = "\\[|\\]")) %>%
  separate(y, into = c("group", "occasion"), sep = ",") %>% 
  mutate(group = as.numeric(group),
         occasion = as.numeric(occasion)) %>% 
  left_join(Grid, by = c("group" = "G")) %>% 
  left_join(dat_date, by = "occasion") %>% 
  mutate(cohort = case_when(CO == 1 ~ "2018",
                            CO == 2 ~ "2019",
                            CO == 3 ~ "2020",
                            CO == 4 ~ "2021"),
         age2_maturity = MT,
         year = year(date),
         month = month(date),
         age = year - as.numeric(cohort)) %>% 
  left_join(count, by = c("CO","CL","MT"))

# save to csv
write.csv(result, out_name)

#### make figure ---------------------------------------------------------------

## 2018 ~ 2021 cohort

# filtering
dat_sub <- result %>% 
  filter(age == 3) %>% 
  mutate(growth = as.character(CL))

# color list
color_list = brewer.pal(4, "RdYlBu")

g <- ggplot(dat_sub, aes(x = growth, y = mean, ymin = lower, ymax = upper, 
                         colour = growth, group = age2_maturity, shape = age2_maturity)) +
  geom_errorbar(width = 0, linewidth = 1, position = position_dodge(.7)) +
  geom_point(size = 4, position = position_dodge(.7)) +
  facet_grid(~cohort) +
  ylab("Survival probability") +
  xlab("Glowth cluster") +
  theme_bw() +
  scale_color_manual('glowth cluster: ', values = color_list) +
  scale_shape_manual('age-2 maturity: ', values = c(16,17)) +
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
       file = paste("output/figure_phi_age2_3_growth_cohort_age2mature_v2.png", sep = ""),
       dpi = 200, width = 13, height = 5)


### 3-3-2 CO and MT --------------------------------------------------------
## Grouped by maturity at age 1 or not
## use annual captured data
## use 2018 - 2021 cohorts

out_name <- "output/phi_age2_3_age2mature_v2.csv"

# filtering by growth cluster
dat1 <- data2 %>% 
  filter(age2_mature != "nd",
         cohort >= 2018) %>% 
  mutate(Occasion = as.numeric(factor(year))) # edit occasion

unique(dat1$cohort)

# check y-m
occasion_list <- tibble(
  ym = sort(unique(dat1$occasion)),
  oc = sort(unique((dat1$Occasion)))
)
occasion_list

## create history matrices
## - Variable transformation
## - `Occasion` - `year`
## - `CL` - `growth` 
## - `CO` - `cohort`
## - `MT` - `age at maturity` 
## Note: An error will occur if there is data with duplicate ID and captured date.

## CH: Capture history - binary information of recaptured or not
CH <- dat1 %>%
  mutate(capture = 1) %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = capture,
              id_cols = c(id, cohort, age2_mature),
              values_fill = list(capture = 0)) %>%
  rename(CO = cohort,
         MT = age2_mature) %>%
  mutate(CO = as.numeric(factor(CO)),
         MT = as.factor(MT)) %>%
  arrange(id) # arrange by id

# check number of samples
count <- count(CH, CO, MT)
count

# select sorted occasion data
ch <- CH %>%
  dplyr::select(as.factor(sort(unique(as.numeric(dat1$Occasion)))))

# put NA until the first capture
for(i in 1:nrow(ch)) {
  x <- min(which(ch[i,] == 1)) - 1
  if(x != 0) ch[i, 1:x] <- NA
}

CH[, which(colnames(CH) %in% sort(unique(dat1$Occasion)))] <- ch

print(CH)

## DH: Location history - distance from the downstream end to the upstream landmark of a (re)capture subsection
DH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = SecUL,
              id_cols = c(id, cohort, age2_mature)) %>%
  rename(CO = cohort,
         MT = age2_mature) %>%
  mutate(CO = as.numeric(factor(CO)),
         MT = as.factor(MT)) %>%
  arrange(id)

print(DH)

## JH: Julian history - Julian date of capture for each occasion
JH <- dat1 %>%
  arrange(Occasion) %>%
  pivot_wider(names_from = Occasion,
              values_from = Julian,
              id_cols = c(id, cohort, age2_mature)) %>%
  rename(CO = cohort,
         MT = age2_mature) %>%
  mutate(CO = as.numeric(factor(CO)),
         MT = as.factor(MT)) %>%
  arrange(id)

print(JH)

J <- as.matrix(JH[,which(colnames(JH) == "1"):ncol(JH)]) # make matrix of JH

#### mcmc setting --------------------------------------------------------------

n.ad <- 100
n.iter <- 2E+4
n.thin <- max(3, ceiling(n.iter/500))
burn <- ceiling(max(10, n.iter/2))
Sample <- ceiling(n.iter/n.thin)

# parameter list
para <- c("xi", "mu.p", "sigma.p", "p", "pi", "phi", "alpha")

#### bayesian inference --------------------------------------------------------

## Data for JAGS
### Response
X <- as.matrix(DH[,which(colnames(DH) == "1"):ncol(DH)] - 5)
Y <- z <- as.matrix(CH[,which(colnames(CH) == "1"):ncol(CH)])
z[Y == 0] <- NA # Y data replaced 0 to NA
ObsF <- apply(Y, 1, function(x) min(which(is.na(x) == 0))) # first observed occasion 
ObsL <- apply(Y, 1, function(x) max(which(x == 1))) # last observed occasion 

# put 1 first to last observed occasion 
for(i in 1:length(ObsF)){
  z[i,ObsF[i]:ObsL[i]] <- 1
}

### Explanatory
Nday <- diff(apply(J, 2, median, na.rm = T)) # interval between occasions
Grid <- count %>% dplyr::select(CO,MT) # Create unique combinations
Grid$G <- 1:length(Grid$CO) # make new no. column
Grid$id <- as.numeric(paste0(Grid$CO,Grid$MT)) # make new id column

# marge Grid and Capture history data (Grouping)
groupdat <- tibble(id = as.numeric(paste0(CH$CO, CH$MT))) %>%
  left_join(Grid, by = "id")

G <- groupdat$G
first_occasion <- Grid %>% 
  mutate(Fo = case_when(CO == 1 ~ 2,
                        CO == 2 ~ 3,
                        CO == 3 ~ 4,
                        CO == 4 ~ 5))

# make input data
Djags <- list(X = X,
              Y = Y,
              Nday = Nday,
              z = z,
              Nind = nrow(Y), # Number of sample
              Nt = ncol(Y), # Number of occasion
              Ng = length(unique(G)), # Number of unique combinations for groups
              G = G, # grouping data
              Fo = first_occasion$Fo,
              ObsF = ObsF
)


# Repeat an equation 3 times
inits <- replicate(3,
                   list(logit.p = matrix(0, nrow = Djags$Ng, ncol = Djags$Nt-1),
                        .RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA ),
                   simplify = F)

for(k in 1:3) inits[[k]]$.RNG.seed <- k

#### mcmc using runjags --------------------------------------------------------
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

#### output --------------------------------------------------------------------

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
result <- bpost %>% 
  filter(str_detect(param, "phi")) %>% 
  mutate(x = str_extract(param, pattern = "\\[.{1,}\\]"),
         y = str_remove_all(x, pattern = "\\[|\\]")) %>%
  separate(y, into = c("group", "occasion"), sep = ",") %>% 
  mutate(group = as.numeric(group),
         occasion = as.numeric(occasion)) %>% 
  left_join(Grid, by = c("group" = "G")) %>% 
  left_join(dat_date, by = "occasion") %>% 
  mutate(cohort = case_when(CO == 1 ~ "2018",
                            CO == 2 ~ "2019",
                            CO == 3 ~ "2020",
                            CO == 4 ~ "2021"),
         age2_maturity = MT,
         year = year(date),
         month = month(date),
         age = year - as.numeric(cohort)) %>% 
  left_join(count, by = c("CO","MT"))  %>% 
  mutate(n = case_when(cohort == "2018" & age2_maturity == 0 ~ paste0('n = ', as.character(n)),
                                                            TRUE ~ as.character(n)))

# save to csv
write.csv(result, out_name)

#### make figure ---------------------------------------------------------------

## 2018 ~ 2021 cohort

# filtering
dat_sub <- result %>% 
  filter(age == 3)

g <- ggplot(dat_sub, aes(x = age2_maturity, y = mean, ymin = lower, ymax = upper)) +
  geom_errorbar(width = 0, linewidth = 1, position = position_dodge(.7)) +
  geom_point(size = 4, position = position_dodge(.7)) +
  facet_grid(~cohort) +
  ylab("Survival probability") +
  xlab("Age 2 maturity") +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1.2), breaks=seq(0,1.5,0.5)) +
  geom_text(aes(x = age2_maturity, y = 1.1, label = n), size = 4, position = position_dodge(.7)) +
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
       file = paste("output/figure_phi_age2_3_cohort_age2mature_v2.png", sep = ""),
       dpi = 200, width = 13, height = 4)



# ...----


# end --------------------------------------------------------------------------