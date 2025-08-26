
# ============================================================================ #
# ========       LONGITUDINAL CLUSTERING METHODS COMPARISON         ========== #
# ========              Analysis with simulated data                 ========= #
# ============================================================================ #

rm(list = ls())

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary packages:
library(longitudinalData)
library(ggplot2)
library(factoextra)
library(knitr)
library(papeR)
library(gridExtra)
library(ggpubr)
library(lcmm)
library(fpc)
library(RColorBrewer)
library(tidyr)
library(ggvenn)
library(FactoMineR)
library(factoextra)
library(kableExtra)
library(dplyr)
library(openxlsx)
library(qgraph)
library(networkD3)
library(dplyr)
library(splines)
library(lme4) # for ranef(), extract residuals of random factors
library(iccTraj)
library(tidyverse)
library(psych)
library(gtools) 

#Load my functions:
source("3_simulations_functions.R")

# -------------------------------------- #
# 1. DATA SIMULATION SPECIFICATIONS ---- 
# -------------------------------------- #

# Respon: al final del seguiment, la puntuació PANSS ha baixat 10 punts o més
# No respon: al final del seguiment, la puntuació PANSS ha baixat menys de 10 punts
# Simulo 240 individus com en la mostra original, però igual proporció de n en tots els grups

n_ind <- 240  # Número total de individuos

## 1.1. Tres grups ----
prop_grupos <- c(NonResp = 0.3333333, Fast_Resp = 0.3333333, Slow_Resp = 0.3333333)
n_grupos <- round(n_ind * prop_grupos)


### 3 separats ----
set.seed(123)
sim_data_3_sep <- bind_rows(
  simular_grupo(beta0=30, slopes= c(4, 3, 0, 0), sd_b =1 , sd_epsilon=1, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "NonResp", n=n_grupos["NonResp"]),
  simular_grupo(beta0=30, slopes= c(-1, -6, -7, -9), sd_b =1 , sd_epsilon=1, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Fast_Resp", n=n_grupos["Fast_Resp"]),
  simular_grupo(beta0=30, slopes= c(-12, -14, -15, -19), sd_b =1 , sd_epsilon=1, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Slow_Resp", n=n_grupos["Slow_Resp"]))
sim_data_3_sep$group <- factor(sim_data_3_sep$group, levels = c("NonResp", "Fast_Resp", "Slow_Resp"))
icc_3_sep <- calculate_icc(sim_data_3_sep)

### 1 separat i 2 solapats ----
set.seed(123)
sim_data_3_sep_sol <- bind_rows(
  simular_grupo(beta0=30, slopes= c(4, 3, 0, 0), sd_b =1 , sd_epsilon=2, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "NonResp", n=n_grupos["NonResp"]),
  simular_grupo(beta0=30, slopes= c(-6, -12, -13, -15), sd_b =2 , sd_epsilon=3, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Fast_Resp", n=n_grupos["Fast_Resp"]),
  simular_grupo(beta0=30, slopes= c(-12, -13, -15, -19), sd_b =2 , sd_epsilon=3, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Slow_Resp", n=n_grupos["Slow_Resp"]))
sim_data_3_sep_sol$group <- factor(sim_data_3_sep_sol$group, levels = c("NonResp", "Fast_Resp", "Slow_Resp"))
icc_3_sep_sol <- calculate_icc(sim_data_3_sep_sol)

### 3 solapats ----
set.seed(123)
sim_data_3_sol <- bind_rows(
  simular_grupo(beta0=30, slopes= c(-5, -5, -6, -6), sd_b =2 , sd_epsilon=4, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "NonResp", n=n_grupos["NonResp"]),
  simular_grupo(beta0=30, slopes= c(-6, -10, -11, -13), sd_b =2 , sd_epsilon=4, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Fast_Resp", n=n_grupos["Fast_Resp"]),
  simular_grupo(beta0=30, slopes= c(-12, -14, -15, -15), sd_b =2 , sd_epsilon=4, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Slow_Resp", n=n_grupos["Slow_Resp"]))
sim_data_3_sol$group <- factor(sim_data_3_sol$group, levels = c("NonResp", "Fast_Resp", "Slow_Resp"))
icc_3_sol <- calculate_icc(sim_data_3_sol)




## 1.2. Quatre grups ----
prop_grupos <- c(Severe_NonResp = 0.25, Moderate_NonResp = 0.25, Fast_Resp = 0.25, Slow_Resp = 0.25)
n_grupos <- round(n_ind * prop_grupos)

### 4 separats ----
set.seed(123)
sim_data_4_sep <- bind_rows(
  simular_grupo(beta0=30, slopes= c(4, 3, 3, 2.5), sd_b =1 , sd_epsilon=1, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Severe_NonResp", n=n_grupos["Severe_NonResp"]),
  simular_grupo(beta0=30, slopes= c(-3.5, -9, -9.5, -10), sd_b =1 , sd_epsilon=1, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Fast_Resp", n=n_grupos["Fast_Resp"]),
  simular_grupo(beta0=30, slopes= c(-12, -14, -15, -19), sd_b =1 , sd_epsilon=1, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Slow_Resp", n=n_grupos["Slow_Resp"]),
  simular_grupo(beta0=30, slopes= c(-3, -3.5, -3.1, -3.6), sd_b =1 , sd_epsilon=1, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Moderate_NonResp", n=n_grupos["Moderate_NonResp"]))
sim_data_4_sep$group <- factor(sim_data_4_sep$group, levels = c("Severe_NonResp", "Fast_Resp", "Slow_Resp", "Moderate_NonResp"))
icc_4_sep <- calculate_icc(sim_data_4_sep)


### 2 separats i 2 solapats ----
set.seed(123)
sim_data_4_sep_sol <- bind_rows(
  simular_grupo(beta0=30, slopes= c(4, 3, 3, 2.5), sd_b =2 , sd_epsilon=3, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Severe_NonResp", n=n_grupos["Severe_NonResp"]),
  simular_grupo(beta0=30, slopes= c(-7, -13, -14, -16), sd_b =2 , sd_epsilon=3, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Fast_Resp", n=n_grupos["Fast_Resp"]),
  simular_grupo(beta0=30, slopes= c(-13, -14, -16, -20), sd_b =2 , sd_epsilon=3, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Slow_Resp", n=n_grupos["Slow_Resp"]),
  simular_grupo(beta0=30, slopes= c(0, -0.5, -0.1, -0.6), sd_b =2 , sd_epsilon=3, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Moderate_NonResp", n=n_grupos["Moderate_NonResp"]))
sim_data_4_sep_sol$group <- factor(sim_data_4_sep_sol$group, levels = c("Severe_NonResp", "Fast_Resp", "Slow_Resp", "Moderate_NonResp"))
icc_4_sep_sol <- calculate_icc(sim_data_4_sep_sol)



### 4 solapats ----
set.seed(123)
sim_data_4_sol <- bind_rows(
  simular_grupo(beta0=30, slopes= c(-4, -4, -4.5, -4.5), sd_b =2 , sd_epsilon=4, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Severe_NonResp", n=n_grupos["Severe_NonResp"]),
  simular_grupo(beta0=30, slopes= c(-6, -7.5, -7.5, -7), sd_b =2 , sd_epsilon=4, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Moderate_NonResp", n=n_grupos["Moderate_NonResp"]),
  simular_grupo(beta0=30, slopes= c(-6, -12, -12, -13), sd_b =2 , sd_epsilon=4, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Fast_Resp", n=n_grupos["Fast_Resp"]),
  simular_grupo(beta0=30, slopes= c(-11, -13, -15, -16), sd_b =2 , sd_epsilon=4, ymin=7, ymax=49,
                time= c(0,2,6,12), class_name= "Slow_Resp", n=n_grupos["Slow_Resp"]))
sim_data_4_sol$group <- factor(sim_data_4_sol$group, levels = c("Severe_NonResp", "Fast_Resp", "Slow_Resp", "Moderate_NonResp"))
icc_4_sol <- calculate_icc(sim_data_4_sol)


# ------------- #
# 2. PLOTS ----
# ------------- #

## 2.1. Tres grups  ----
p_3_sep <- ggplot(sim_data_3_sep, aes(x = time, y = PANSS, color = group, group = id)) +
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = group), geom = "line", fun = mean, linewidth = 1.5) +
  scale_x_continuous(breaks = c(0,2,6,12)) +
  labs(title = "Three groups", subtitle = paste0("All separated (ICC=", round(icc_3_sep$r,2),")"),
       x = "Time (months)", y = "Score") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set3") + 
  scale_y_continuous(limits = c(5,40),breaks = seq(10, 40, 10)) 

p_3_sep_sol <- ggplot(sim_data_3_sep_sol, aes(x = time, y = PANSS, color = group, group = id)) +
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = group), geom = "line", fun = mean, linewidth = 1.5) +
  scale_x_continuous(breaks = c(0,2,6,12)) +
  labs(title = " ", subtitle = paste0("One separated and two overlapped (ICC=", round(icc_3_sep_sol$r,2),")"),
       x = "Time (months)", y = "Score") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set3") + 
  scale_y_continuous(limits = c(5,40),breaks = seq(10, 40, 10)) 

p_3_sol <- ggplot(sim_data_3_sol, aes(x = time, y = PANSS, color = group, group = id)) +
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = group), geom = "line", fun = mean, linewidth = 1.5) +
  scale_x_continuous(breaks = c(0,2,6,12)) +
  labs(title = " ", subtitle = paste0("All overlapped (ICC=", round(icc_3_sol$r,2),")"),
       x = "Time (months)", y = "Score") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set3") + 
  scale_y_continuous(limits = c(5,40),breaks = seq(10, 40, 10)) 

## 2.2. Quatre grups ----
p_4_sep <- ggplot(sim_data_4_sep, aes(x = time, y = PANSS, color = group, group = id)) +
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = group), geom = "line", fun = mean, linewidth = 1.5) +
  scale_x_continuous(breaks = c(0,2,6,12)) +
  labs(title = "Four groups", subtitle = paste0("All separated (ICC=", round(icc_4_sep$r,2),")"),
       x = "Time (months)", y = "Score") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set3") + 
  scale_y_continuous(limits = c(5,40),breaks = seq(10, 40, 10)) 

p_4_sep_sol <- ggplot(sim_data_4_sep_sol, aes(x = time, y = PANSS, color = group, group = id)) +
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = group), geom = "line", fun = mean, linewidth = 1.5) +
  scale_x_continuous(breaks = c(0,2,6,12)) +
  labs(title = " ", subtitle = paste0("Two separated and two overlapped (ICC=", round(icc_4_sep_sol$r,2),")"),
       x = "Time (months)", y = "Score") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set3") + 
  scale_y_continuous(limits = c(5,40),breaks = seq(10, 40, 10)) 

p_4_sol <- ggplot(sim_data_4_sol, aes(x = time, y = PANSS, color = group, group = id)) +
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = group), geom = "line", fun = mean, linewidth = 1.5) +
  scale_x_continuous(breaks = c(0,2,6,12)) +
  labs(title = " ", subtitle = paste0("All overlapped (ICC=", round(icc_4_sol$r,2),")"),
       x = "Time (months)", y = "Score") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set3") + 
  scale_y_continuous(limits = c(5,40),breaks = seq(10, 40, 10)) 

# Save plots:
pdf(paste0(path, "/Plots_simulations.pdf"), pointsize = 22, width = 19, height = 10)
grid.arrange(p_3_sep, p_3_sep_sol, p_3_sol, ncol=3)
grid.arrange(p_4_sep, p_4_sep_sol, p_4_sol, ncol=3)
dev.off()



# ------------------- #
# 3. Kappa index ----
# ------------------- #

## 3.1. Tres grups ----
k=3
### 3.1.1. 3 separats ----
#### Clustering techniques
sim_data_3_sep_wide <- sim_data_3_sep %>% select(id, group, time, PANSS) %>%
  pivot_wider(names_from = time, values_from = PANSS, names_prefix = "PANSS_") %>% 
  distinct() %>% mutate(group_num = as.numeric(group))
D <- dist(sim_data_3_sep_wide[,3:6])
set.seed(123)
clust.hier <- hclust(D, method = "ward.D2")
sim_data_3_sep_wide$group.hier <- cutree(clust.hier, k)
sim_data_3_sep_wide$group.kmeans <- kmeans(D, centers = k)$cluster
sim_data_3_sep_wide$group.pam <- pam(D, k = k, diss = TRUE)$clustering

kappa_clust_3_sep <- calcular_kappa_clusters(sim_data_3_sep_wide, cluster_cols = c("group.hier", "group.kmeans", "group.pam"), k=k)

#### LCMMs
sim_data_3_sep_wide <- transform(sim_data_3_sep_wide, id = as.numeric(factor(id)))
sim_data_3_sep <- transform(sim_data_3_sep, id = as.numeric(factor(id))) %>% mutate(group_num = as.numeric(group))
set.seed(123)
lcmm.linear <- hlme(PANSS ~ time, random = ~1, subject = "id", data = sim_data_3_sep, ng = 1)
lcmm.quad <- hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = sim_data_3_sep, ng = 1)
lcmm.splines <- hlme(PANSS ~ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = sim_data_3_sep, ng = 1)
x <- hlme(PANSS ~ time, random = ~1, subject = "id", data = sim_data_3_sep, mixture = ~ time, ng = k, B = random(lcmm.linear))
sim_data_3_sep_wide <-  merge(sim_data_3_sep_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.lin = class)
x <- hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = sim_data_3_sep, mixture = ~ time + I(time^2), ng = k, B = random(lcmm.quad))
sim_data_3_sep_wide <-  merge(sim_data_3_sep_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.quad= class)
x <- hlme(PANSS ~ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = sim_data_3_sep, mixture = ~ns(time, knots = c(2, 6)), ng = k, B = random(lcmm.splines))
sim_data_3_sep_wide <-  merge(sim_data_3_sep_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.spl = class)

kappa_lcmm_3_sep <- calcular_kappa_clusters(sim_data_3_sep_wide, cluster_cols = c("group.lcmm.lin", "group.lcmm.quad", "group.lcmm.spl"), k=k)


### 3.1.2. 1 separat i 2 solapats ----
#### Clustering techniques
sim_data_3_sep_sol_wide <- sim_data_3_sep_sol %>% select(id, group, time, PANSS) %>%
  pivot_wider(names_from = time, values_from = PANSS, names_prefix = "PANSS_") %>% 
  distinct() %>% mutate(group_num = as.numeric(group))
D <- dist(sim_data_3_sep_sol_wide[,3:6])
set.seed(123)
clust.hier <- hclust(D, method = "ward.D2")
sim_data_3_sep_sol_wide$group.hier <- cutree(clust.hier, k)
sim_data_3_sep_sol_wide$group.kmeans <- kmeans(D, centers = k)$cluster
sim_data_3_sep_sol_wide$group.pam <- pam(D, k = k, diss = TRUE)$clustering

kappa_clust_3_sep_sol <- calcular_kappa_clusters(sim_data_3_sep_sol_wide, cluster_cols = c("group.hier", "group.kmeans", "group.pam"), k=k)

#### LCMMs
sim_data_3_sep_sol_wide <- transform(sim_data_3_sep_sol_wide, id = as.numeric(factor(id)))
sim_data_3_sep_sol <- transform(sim_data_3_sep_sol, id = as.numeric(factor(id))) %>% mutate(group_num = as.numeric(group))
set.seed(123)
lcmm.linear <- hlme(PANSS ~ time, random = ~1, subject = "id", data = sim_data_3_sep_sol, ng = 1)
lcmm.quad <- hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = sim_data_3_sep_sol, ng = 1)
lcmm.splines <- hlme(PANSS ~ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = sim_data_3_sep_sol, ng = 1)
x <- hlme(PANSS ~ time, random = ~1, subject = "id", data = sim_data_3_sep_sol, mixture = ~ time, ng = k, B = random(lcmm.linear))
sim_data_3_sep_sol_wide <-  merge(sim_data_3_sep_sol_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.lin = class)
x <- hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = sim_data_3_sep_sol, mixture = ~ time + I(time^2), ng = k, B = random(lcmm.quad))
sim_data_3_sep_sol_wide <-  merge(sim_data_3_sep_sol_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.quad= class)
x <- hlme(PANSS ~ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = sim_data_3_sep_sol, mixture = ~ns(time, knots = c(2, 6)), ng = k, B = random(lcmm.splines))
sim_data_3_sep_sol_wide <-  merge(sim_data_3_sep_sol_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.spl = class)

kappa_lcmm_3_sep_sol <- calcular_kappa_clusters(sim_data_3_sep_sol_wide, cluster_cols = c("group.lcmm.lin", "group.lcmm.quad", "group.lcmm.spl"), k=k)



### 3.1.3. 3 solapats ----
#### Clustering techniques
sim_data_3_sol_wide <- sim_data_3_sol %>% select(id, group, time, PANSS) %>%
  pivot_wider(names_from = time, values_from = PANSS, names_prefix = "PANSS_") %>% 
  distinct() %>% mutate(group_num = as.numeric(group))
D <- dist(sim_data_3_sol_wide[,3:6])
set.seed(123)
clust.hier <- hclust(D, method = "ward.D2")
sim_data_3_sol_wide$group.hier <- cutree(clust.hier, k)
sim_data_3_sol_wide$group.kmeans <- kmeans(D, centers = k)$cluster
sim_data_3_sol_wide$group.pam <- pam(D, k = k, diss = TRUE)$clustering

kappa_clust_3_sol <- calcular_kappa_clusters(sim_data_3_sol_wide, cluster_cols = c("group.hier", "group.kmeans", "group.pam"), k=k)

#### LCMMs
sim_data_3_sol_wide <- transform(sim_data_3_sol_wide, id = as.numeric(factor(id)))
sim_data_3_sol <- transform(sim_data_3_sol, id = as.numeric(factor(id))) %>% mutate(group_num = as.numeric(group))
set.seed(123)
lcmm.linear <- hlme(PANSS ~ time, random = ~1, subject = "id", data = sim_data_3_sol, ng = 1)
lcmm.quad <- hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = sim_data_3_sol, ng = 1)
lcmm.splines <- hlme(PANSS ~ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = sim_data_3_sol, ng = 1)
x <- hlme(PANSS ~ time, random = ~1, subject = "id", data = sim_data_3_sol, mixture = ~ time, ng = k, B = random(lcmm.linear))
sim_data_3_sol_wide <-  merge(sim_data_3_sol_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.lin = class)
x <- hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = sim_data_3_sol, mixture = ~ time + I(time^2), ng = k, B = random(lcmm.quad))
sim_data_3_sol_wide <-  merge(sim_data_3_sol_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.quad= class)
x <- hlme(PANSS ~ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = sim_data_3_sol, mixture = ~ns(time, knots = c(2, 6)), ng = k, B = random(lcmm.splines))
sim_data_3_sol_wide <-  merge(sim_data_3_sol_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.spl = class)

kappa_lcmm_3_sol <- calcular_kappa_clusters(sim_data_3_sol_wide, cluster_cols = c("group.lcmm.lin", "group.lcmm.quad", "group.lcmm.spl"), k=k)


## 3.2. Quatre grups ----
k=4
### 3.2.1. 4 separats ----
#### Clustering techniques
sim_data_4_sep_wide <- sim_data_4_sep %>% select(id, group, time, PANSS) %>%
  pivot_wider(names_from = time, values_from = PANSS, names_prefix = "PANSS_") %>% 
  distinct() %>% mutate(group_num = as.numeric(group))
D <- dist(sim_data_4_sep_wide[,3:6])
set.seed(123)
clust.hier <- hclust(D, method = "ward.D2")
sim_data_4_sep_wide$group.hier <- cutree(clust.hier, k)
sim_data_4_sep_wide$group.kmeans <- kmeans(D, centers = k)$cluster
sim_data_4_sep_wide$group.pam <- pam(D, k = k, diss = TRUE)$clustering

kappa_clust_4_sep <- calcular_kappa_clusters(sim_data_4_sep_wide, cluster_cols = c("group.hier", "group.kmeans", "group.pam"), k=k)

#### LCMMs
sim_data_4_sep_wide <- transform(sim_data_4_sep_wide, id = as.numeric(factor(id)))
sim_data_4_sep <- transform(sim_data_4_sep, id = as.numeric(factor(id))) %>% mutate(group_num = as.numeric(group))
set.seed(123)
lcmm.linear <- hlme(PANSS ~ time, random = ~1, subject = "id", data = sim_data_4_sep, ng = 1)
lcmm.quad <- hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = sim_data_4_sep, ng = 1)
lcmm.splines <- hlme(PANSS ~ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = sim_data_4_sep, ng = 1)
x <- hlme(PANSS ~ time, random = ~1, subject = "id", data = sim_data_4_sep, mixture = ~ time, ng = k, B = random(lcmm.linear))
sim_data_4_sep_wide <-  merge(sim_data_4_sep_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.lin = class)
x <- hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = sim_data_4_sep, mixture = ~ time + I(time^2), ng = k, B = random(lcmm.quad))
sim_data_4_sep_wide <-  merge(sim_data_4_sep_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.quad= class)
x <- hlme(PANSS ~ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = sim_data_4_sep, mixture = ~ns(time, knots = c(2, 6)), ng = k, B = random(lcmm.splines))
sim_data_4_sep_wide <-  merge(sim_data_4_sep_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.spl = class)

kappa_lcmm_4_sep <- calcular_kappa_clusters(sim_data_4_sep_wide, cluster_cols = c("group.lcmm.lin", "group.lcmm.quad", "group.lcmm.spl"), k=k)



### 3.2.2. 2 separat i 2 solapats ----
#### Clustering techniques
sim_data_4_sep_sol_wide <- sim_data_4_sep_sol %>% select(id, group, time, PANSS) %>%
  pivot_wider(names_from = time, values_from = PANSS, names_prefix = "PANSS_") %>% 
  distinct() %>% mutate(group_num = as.numeric(group))
D <- dist(sim_data_4_sep_sol_wide[,3:6])
set.seed(123)
clust.hier <- hclust(D, method = "ward.D2")
sim_data_4_sep_sol_wide$group.hier <- cutree(clust.hier, k)
sim_data_4_sep_sol_wide$group.kmeans <- kmeans(D, centers = k)$cluster
sim_data_4_sep_sol_wide$group.pam <- pam(D, k = k, diss = TRUE)$clustering

kappa_clust_4_sep_sol <- calcular_kappa_clusters(sim_data_4_sep_sol_wide, cluster_cols = c("group.hier", "group.kmeans", "group.pam"), k=k)

#### LCMMs
sim_data_4_sep_sol_wide <- transform(sim_data_4_sep_sol_wide, id = as.numeric(factor(id)))
sim_data_4_sep_sol <- transform(sim_data_4_sep_sol, id = as.numeric(factor(id))) %>% mutate(group_num = as.numeric(group))
set.seed(123)
lcmm.linear <- hlme(PANSS ~ time, random = ~1, subject = "id", data = sim_data_4_sep_sol, ng = 1)
lcmm.quad <- hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = sim_data_4_sep_sol, ng = 1)
lcmm.splines <- hlme(PANSS ~ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = sim_data_4_sep_sol, ng = 1)
x <- hlme(PANSS ~ time, random = ~1, subject = "id", data = sim_data_4_sep_sol, mixture = ~ time, ng = k, B = random(lcmm.linear))
sim_data_4_sep_sol_wide <-  merge(sim_data_4_sep_sol_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.lin = class)
x <- hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = sim_data_4_sep_sol, mixture = ~ time + I(time^2), ng = k, B = random(lcmm.quad))
sim_data_4_sep_sol_wide <-  merge(sim_data_4_sep_sol_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.quad= class)
x <- hlme(PANSS ~ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = sim_data_4_sep_sol, mixture = ~ns(time, knots = c(2, 6)), ng = k, B = random(lcmm.splines))
sim_data_4_sep_sol_wide <-  merge(sim_data_4_sep_sol_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.spl = class)

kappa_lcmm_4_sep_sol <- calcular_kappa_clusters(sim_data_4_sep_sol_wide, cluster_cols = c("group.lcmm.lin", "group.lcmm.quad", "group.lcmm.spl"), k=k)


### 3.2.3. 4 solapats ----
#### Clustering techniques
sim_data_4_sol_wide <- sim_data_4_sol %>% select(id, group, time, PANSS) %>%
  pivot_wider(names_from = time, values_from = PANSS, names_prefix = "PANSS_") %>% 
  distinct() %>% mutate(group_num = as.numeric(group))
D <- dist(sim_data_4_sol_wide[,3:6])
set.seed(123)
clust.hier <- hclust(D, method = "ward.D2")
sim_data_4_sol_wide$group.hier <- cutree(clust.hier, k)
sim_data_4_sol_wide$group.kmeans <- kmeans(D, centers = k)$cluster
sim_data_4_sol_wide$group.pam <- pam(D, k = k, diss = TRUE)$clustering

kappa_clust_4_sol <- calcular_kappa_clusters(sim_data_4_sol_wide, cluster_cols = c("group.hier", "group.kmeans", "group.pam"), k=k)


#### LCMMs
sim_data_4_sol_wide <- transform(sim_data_4_sol_wide, id = as.numeric(factor(id)))
sim_data_4_sol <- transform(sim_data_4_sol, id = as.numeric(factor(id))) %>% mutate(group_num = as.numeric(group))
set.seed(123)
lcmm.linear <- hlme(PANSS ~ time, random = ~1, subject = "id", data = sim_data_4_sol, ng = 1)
lcmm.quad <- hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = sim_data_4_sol, ng = 1)
lcmm.splines <- hlme(PANSS ~ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = sim_data_4_sol, ng = 1)
x <- hlme(PANSS ~ time, random = ~1, subject = "id", data = sim_data_4_sol, mixture = ~ time, ng = k, B = random(lcmm.linear))
sim_data_4_sol_wide <-  merge(sim_data_4_sol_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.lin = class)
x <- hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = sim_data_4_sol, mixture = ~ time + I(time^2), ng = k, B = random(lcmm.quad))
sim_data_4_sol_wide <-  merge(sim_data_4_sol_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.quad= class)
x <- hlme(PANSS ~ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = sim_data_4_sol, mixture = ~ns(time, knots = c(2, 6)), ng = k, B = random(lcmm.splines))
sim_data_4_sol_wide <-  merge(sim_data_4_sol_wide, x$pprob[,1:2], by = "id") %>% 
  rename(group.lcmm.spl = class)

kappa_lcmm_4_sol <- calcular_kappa_clusters(sim_data_4_sol_wide, cluster_cols = c("group.lcmm.lin", "group.lcmm.quad", "group.lcmm.spl"), k=k)





# ------------------------------- #
# 4. SIMULATIONS AND METRICS ----
# ------------------------------- #
n_ind <- 240

## 4.1. Tres grups ----
prop_grupos <- c(NonResp = 0.3333333, Fast_Resp = 0.3333333, Slow_Resp = 0.3333333)

### 4.1.1. 3 separats ----
#### Clustering techniques
start_time <- Sys.time()
set.seed(123)
results_clust_3_sep <- perform_clustering_simulation(n_ind = n_ind, prop_groups = prop_grupos, beta0_values= c(30, 30, 30),
                                              slopes_values = list(c(4, 3, 0, 0), c(-1, -6, -7, -9), c(-12, -14, -15, -19)),
                                              sd_b=c(1,1,1), sd_epsilon=c(1,1,1),
                                              ymin=c(7,7,7), ymax=c(49,49,49), format= "wide", time= c(0,2,6,12),
                                              num_resamples = 500, k_range = 2:5, trace = TRUE, min.percentage = 5)
end_time <- Sys.time()
time_clust_3_sep <- end_time - start_time


### LCMMs
start_time <- Sys.time()
set.seed(123)
results_lcmm_3_sep <- perform_lcmm_simulation(n_ind = n_ind, prop_groups = prop_grupos, beta0_values= c(30, 30, 30),
                                              slopes_values = list(c(4, 3, 0, 0), c(-1, -6, -7, -9), c(-12, -14, -15, -19)),
                                              sd_b=c(1,1,1), sd_epsilon=c(1,1,1),
                                              ymin=c(7,7,7), ymax=c(49,49,49), format= "long", time= c(0,2,6,12),
                                              num_resamples = 500, k_range = 2:5, trace = TRUE, min.percentage = 5)
end_time <- Sys.time()
time_lcmm_3_sep <- end_time - start_time


### 4.1.2. 1 separat i 2 solapats ----
#### Clustering techniques
start_time <- Sys.time()
set.seed(123)
results_clust_3_sep_sol <- perform_clustering_simulation(n_ind = n_ind, prop_groups = prop_grupos, beta0_values= c(30, 30, 30),
                                                     slopes_values = list(c(4, 3, 0, 0), c(-6, -12, -13, -15), c(-12, -13, -15, -19)),
                                                     sd_b=c(1,2,2), sd_epsilon=c(2,3,3),
                                                     ymin=c(7,7,7), ymax=c(49,49,49), format= "wide", time= c(0,2,6,12),
                                                     num_resamples = 500, k_range = 2:5, trace = TRUE, min.percentage = 5)
end_time <- Sys.time()
time_clust_3_sep_sol <- end_time - start_time

### LCMMs
start_time <- Sys.time()
set.seed(123)
results_lcmm_3_sep_sol <- perform_lcmm_simulation(n_ind = n_ind, prop_groups = prop_grupos, beta0_values= c(30, 30, 30),
                                              slopes_values = list(c(4, 3, 0, 0), c(-6, -12, -13, -15), c(-12, -13, -15, -19)),
                                              sd_b=c(1,2,2), sd_epsilon=c(2,3,3),
                                              ymin=c(7,7,7), ymax=c(49,49,49), format= "long", time= c(0,2,6,12),
                                              num_resamples = 500, k_range = 2:5, trace = TRUE, min.percentage = 5)
end_time <- Sys.time()
time_lcmm_3_sep_sol <- end_time - start_time


### 4.1.3. 3 solapats ----
#### Clustering techniques
start_time <- Sys.time()
set.seed(123)
results_clust_3_sol <- perform_clustering_simulation(n_ind = n_ind, prop_groups = prop_grupos, beta0_values= c(30, 30, 30),
                                                     slopes_values = list(c(-5, -5, -6, -6), c(-6, -10, -11, -13), c(-12, -14, -15, -15)),
                                                     sd_b=c(2,2,2), sd_epsilon=c(4,4,4),
                                                     ymin=c(7,7,7), ymax=c(49,49,49), format= "wide", time= c(0,2,6,12),
                                                     num_resamples = 500, k_range = 2:5, trace = TRUE, min.percentage = 5)
end_time <- Sys.time()
time_clust_3_sol <- end_time - start_time

### LCMMs
start_time <- Sys.time()
set.seed(123)
results_lcmm_3_sol <- perform_lcmm_simulation(n_ind = n_ind, prop_groups = prop_grupos, beta0_values= c(30, 30, 30),
                                              slopes_values = list(c(-5, -5, -6, -6), c(-6, -10, -11, -13), c(-12, -14, -15, -15)),
                                              sd_b=c(2,2,2), sd_epsilon=c(4,4,4),
                                              ymin=c(7,7,7), ymax=c(49,49,49), format= "long", time= c(0,2,6,12),
                                              num_resamples = 500, k_range = 2:5, trace = TRUE, min.percentage = 5)
end_time <- Sys.time()
time_lcmm_3_sol <- end_time - start_time


save.image("Results_simulations_170425_k3.RData")
#load("Results_simulations_170425_k3.RData")

## 4.2. Quatre grups ----
prop_grupos <- c(Severe_NonResp = 0.25, Moderate_NonResp = 0.25, Fast_Resp = 0.25, Slow_Resp = 0.25)

### 4.2.1. 4 separats ----
### Clustering techniques
start_time <- Sys.time()
set.seed(123)
results_clust_4_sep <- perform_clustering_simulation(n_ind = n_ind, prop_groups = prop_grupos, beta0_values= c(30, 30, 30, 30),
                                                     slopes_values = list(c(4, 3, 3, 2.5), c(-3, -3.5, -3.1, -3.6), c(-3.5, -9, -9.5, -10), c(-12, -14, -15, -19)),
                                                     sd_b=c(1,1,1,1), sd_epsilon=c(1,1,1,1),
                                                     ymin=c(7,7,7,7), ymax=c(49,49,49,49), format= "wide", time= c(0,2,6,12),
                                                     num_resamples = 500, k_range = 2:5, trace = TRUE, min.percentage = 5)
end_time <- Sys.time()
time_clust_4_sep <- end_time - start_time


### LCMMs
start_time <- Sys.time()
set.seed(123)
results_lcmm_4_sep <- perform_lcmm_simulation(n_ind = n_ind, prop_groups = prop_grupos, beta0_values= c(30, 30, 30, 30),
                                              slopes_values = list(c(4, 3, 3, 2.5), c(-3, -3.5, -3.1, -3.6), c(-3.5, -9, -9.5, -10), c(-12, -14, -15, -19)),
                                              sd_b=c(1,1,1,1), sd_epsilon=c(1,1,1,1),
                                              ymin=c(7,7,7,7), ymax=c(49,49,49,49), format= "long", time= c(0,2,6,12),
                                              num_resamples = 500, k_range = 2:5, trace = TRUE, min.percentage = 5)
end_time <- Sys.time()
time_lcmm_4_sep <- end_time - start_time


### 4.2.2. 2 separats i 2 solapats ----
#### Clustering techniques
start_time <- Sys.time()
set.seed(123)
results_clust_4_sep_sol <- perform_clustering_simulation(n_ind = n_ind, prop_groups = prop_grupos, beta0_values= c(30, 30, 30, 30),
                                                     slopes_values = list(c(4, 3, 3, 2.5), c(0, -0.5, -0.1, -0.6), c(-7, -13, -14, -16), c(-13, -14, -16, -20)),
                                                     sd_b=c(2,2,2,2), sd_epsilon=c(3,3,3,3),
                                                     ymin=c(7,7,7,7), ymax=c(49,49,49,49), format= "wide", time= c(0,2,6,12),
                                                     num_resamples = 500, k_range = 2:5, trace = TRUE, min.percentage = 5)
end_time <- Sys.time()
time_clust_4_sep_sol <- end_time - start_time


### LCMMs
start_time <- Sys.time()
set.seed(123)
results_lcmm_4_sep_sol <- perform_lcmm_simulation(n_ind = n_ind, prop_groups = prop_grupos, beta0_values= c(30, 30, 30, 30),
                                              slopes_values = list(c(4, 3, 3, 2.5), c(0, -0.5, -0.1, -0.6), c(-7, -13, -14, -16), c(-13, -14, -16, -20)),
                                              sd_b=c(2,2,2,2), sd_epsilon=c(3,3,3,3),
                                              ymin=c(7,7,7,7), ymax=c(49,49,49,49), format= "long", time= c(0,2,6,12),
                                              num_resamples = 500, k_range = 2:5, trace = TRUE, min.percentage = 5)
end_time <- Sys.time()
time_lcmm_4_sep_sol <- end_time - start_time


### 4.2.3. 4 solapats ----
#### Clustering techniques
start_time <- Sys.time()
set.seed(123)
results_clust_4_sol <- perform_clustering_simulation(n_ind = n_ind, prop_groups = prop_grupos, beta0_values= c(30, 30, 30, 30),
                                              slopes_values = list(c(-4, -4, -4.5, -4.5), c(-6, -7.5, -7.5, -7), c(-6, -12, -12, -13), c(-11, -13, -15, -16)),
                                              sd_b=c(2,2,2,2), sd_epsilon=c(4,4,4,4),
                                              ymin=c(7,7,7,7), ymax=c(49,49,49,49), format= "wide", time= c(0,2,6,12),
                                              num_resamples = 500, k_range = 2:5, trace = TRUE, min.percentage = 5)
end_time <- Sys.time()
time_clust_4_sol <- end_time - start_time


### LCMMs
start_time <- Sys.time()
set.seed(123)
results_lcmm_4_sol <- perform_lcmm_simulation(n_ind = n_ind, prop_groups = prop_grupos, beta0_values= c(30, 30, 30, 30),
                                              slopes_values = list(c(-4, -4, -4.5, -4.5), c(-6, -7.5, -7.5, -7), c(-6, -12, -12, -13), c(-11, -13, -15, -16)),
                                              sd_b=c(2,2,2,2), sd_epsilon=c(4,4,4,4),
                                              ymin=c(7,7,7,7), ymax=c(49,49,49,49), format= "long", time= c(0,2,6,12),
                                              num_resamples = 500, k_range = 2:5, trace = TRUE, min.percentage = 5)
end_time <- Sys.time()
time_lcmm_4_sol <- end_time - start_time


save.image("Results_simulations_310525.RData")
