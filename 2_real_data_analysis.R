
# ============================================================================ #
# ========       LONGITUDINAL CLUSTERING METHODS COMPARISON         ========== #
# ========                 Analysis with real data                   ========= #
# ============================================================================ #

rm(list = ls())

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary packages
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
library(patchwork)
library(gtools) 
library(fda)

# Load my functions:
source("1_real_data_functions.R")


# ------------------ #
# 1. READ DATA -----
# ------------------ #
data <- read.csv2("../../Data/PEP_Completa.csv") # 588 pacients


## 1.1. Filter observations ----
# Define NAs
data[data== ""] <- NA # lots of elements are "" instead of NA

# Remove controls
data <- data[-which(data$Tipo_sujeto == "Control"),] # eliminem 253 controls = 335 pacients

# Remove observations with edad_estudio < 16 or NA
sum(data$Edad_estudio < 16, na.rm=T) # # Hi ha 9 pacients amb edat VB < 16
sum(is.na(data$Edad_estudio)) # No hi ha pacients sense edat estudi
data <- data[-which(data$Edad_estudio < 16),] # eliminem 9 pacients = 326 pacients

# Remove observations with psicosis tóxicas
sum(data$PSICOSIS_TOXICAS == 1) # 13 pacients amb psicosis toxiques
data <- data[-which(data$PSICOSIS_TOXICAS == 1),] # 326-13 =313 pacients

# Remove afectius
sum(data$FES_SSD_NoAffective_NoToxicPsych == 0) # 51 pacients no afectius
data <- data[-which(data$FES_SSD_NoAffective_NoToxicPsych == 0),] # 313-51 = 262 pacients

# Remove pacients que no prenen antipsicotics en cap moment
id <- which(data$N_Antipsicotico_VB == 0  & data$N_Antipsicotico_V2M == 0 &
              data$N_Antipsicotico_V6M == 0 & data$N_Antipsicotico12M == 0)
length(id) # 16 pacients no prenen antipsicotics
data <- data[-id,] # 262-16 = 246 pacients
rm(id)


# ---------------------------- #
# 2. CREATE SUB-DATASETS -----
# ---------------------------- #

## 2.1. PANSS Negativos ----
data_neg <- data[,c(1,grep("PANSS_negativos", names(data)))]

### Remove pacients que no tenen valor als 2, 6 y 12 mesos
id <- which(is.na(data_neg$PANSS_negativos_V2M) & is.na(data_neg$PANSS_negativos_V6M) & is.na(data_neg$PANSS_negativos_V12M))
data_neg <- data_neg[-id,] # 246-9 = 237 pacients
rm(id)

### Convert data to long format
data_neg_long <- reshape2::melt(data_neg, 
                                id.vars = "Ident_caso", 
                                variable.name = "Visit", 
                                value.name = "PANSS_negativos")

### Prepare data for the Latent Class Models
IDs <- data.frame("Ident_caso" = unique(data_neg_long$Ident_caso))
IDs$ID <- as.numeric(1:nrow(IDs))
data_neg_long <- merge(data_neg_long, IDs, by="Ident_caso")
data_neg_long$Time <- ifelse(data_neg_long$Visit == "PANSS_negativos_VB", 0, 
                             ifelse(data_neg_long$Visit == "PANSS_negativos_V2M", 2, 
                                    ifelse(data_neg_long$Visit == "PANSS_negativos_V6M", 6, 12)))


## 2.2. PANSS Positivos ----
data_pos <- data[,c(1,grep("PANSS_positivos", names(data)))]

### Remove pacients que no tenen valor als 2, 6 y 12 mesos
id <- which(is.na(data_pos$PANSS_positivos_V2M) & is.na(data_pos$PANSS_positivos_V6M) & is.na(data_pos$PANSS_positivos_V12M))
data_pos <- data_pos[-id,] # 246-9 = 237 pacients
rm(id)

### Convert data to long format
data_pos_long <- reshape2::melt(data_pos, 
                                id.vars = "Ident_caso", 
                                variable.name = "Visit", 
                                value.name = "PANSS_positivos")

### Prepare data for the Latent Class Models
IDs <- data.frame("Ident_caso" = unique(data_pos_long$Ident_caso))
IDs$ID <- as.numeric(1:nrow(IDs))
data_pos_long <- merge(data_pos_long, IDs, by="Ident_caso")
data_pos_long$Time <- ifelse(data_pos_long$Visit == "PANSS_positivos_VB", 0, 
                             ifelse(data_pos_long$Visit == "PANSS_positivos_V2M", 2, 
                                    ifelse(data_pos_long$Visit == "PANSS_positivos_V6M", 6, 12)))


# ----------------------------- #
# 3. DESCRIPTIVE ANALYSIS -----
# ----------------------------- #

## 3.1. Density plot and histogram ----
path <- getwd()
pdf(paste0(path, "/Plots_normality.pdf"), paper="a4r", pointsize = 14, width = 35, height = 26)
p1 <- ggdensity(data_neg_long$PANSS_negativos, xlab = "PANSS_negativos")
p2 <- ggplot(data_neg_long, aes(x=PANSS_negativos))+ 
  geom_histogram(color="darkblue", fill="lightblue")+ theme_minimal()
p3 <- ggqqplot(data_neg_long$PANSS_negativos)
p4 <- ggdensity(data_pos_long$PANSS_positivos, xlab = "PANSS_positivos")
p5 <- ggplot(data_pos_long, aes(x=PANSS_positivos))+ 
  geom_histogram(color="darkblue", fill="lightblue")+ theme_minimal()
p6 <- ggqqplot(data_pos_long$PANSS_positivos)

grid.arrange(p1, p2, p3,p4,p5,p6,ncol=3)

aux.neg <- data_neg_long[complete.cases(data_neg_long),] 
aux.pos <- data_pos_long[complete.cases(data_pos_long),] 
aux.neg$PANSS_negativos_norm <- scale(aux.neg$PANSS_negativos)
aux.neg$PANSS_negativos_log<- log(aux.neg$PANSS_negativos)
aux.pos$PANSS_positivos_norm <- scale(aux.pos$PANSS_positivos)
aux.pos$PANSS_positivos_log<- log(aux.pos$PANSS_positivos)

p1 <- ggdensity(aux.neg$PANSS_negativos_norm, xlab = "PANSS_negativos (norm)")
p2 <- ggplot(aux.neg, aes(x=PANSS_negativos_norm))+ 
  geom_histogram(color="darkblue", fill="lightblue")+ theme_minimal()
p3 <- ggqqplot(aux.neg$PANSS_negativos_norm)
p4 <- ggdensity(aux.pos$PANSS_positivos_norm, xlab = "PANSS_positivos (norm)")
p5 <- ggplot(aux.pos, aes(x=PANSS_positivos_norm))+ 
  geom_histogram(color="darkblue", fill="lightblue")+ theme_minimal()
p6 <- ggqqplot(aux.pos$PANSS_positivos_norm)
grid.arrange(p1, p2, p3,p4,p5,p6,ncol=3)


p1 <- ggdensity(aux.neg$PANSS_negativos_log, xlab = "PANSS_negativos (log)")
p2 <- ggplot(aux.neg, aes(x=PANSS_negativos_log))+ 
  geom_histogram(color="darkblue", fill="lightblue")+ theme_minimal()
p3 <- ggqqplot(aux.neg$PANSS_negativos_log)

p4 <- ggdensity(aux.pos$PANSS_positivos_log, xlab = "PANSS_positivos (log)")
p5 <- ggplot(aux.pos, aes(x=PANSS_positivos_log))+ 
  geom_histogram(color="darkblue", fill="lightblue")+ theme_minimal()
p6 <- ggqqplot(aux.pos$PANSS_positivos_log)

grid.arrange(p1, p2, p3,p4,p5,p6,ncol=3)
dev.off()


## 3.2. Extract residuals of the model ----
### 3.2.1. Negative PANSS ----
#### LINEAR MODELS
mod1 <- hlme(PANSS_negativos ~ Time, random = ~ 1, subject = "ID", ng = 1, data = data_neg_long)
mod2 <- hlme(PANSS_negativos ~ Time, random = ~ 1, mixture = ~ Time, subject = "ID", ng = 2, data = data_neg_long, B = random(mod1))
mod3 <- hlme(PANSS_negativos ~ Time, random = ~ 1, mixture = ~ Time, subject = "ID", ng = 3, data = data_neg_long, B = random(mod1))
mod4 <- hlme(PANSS_negativos ~ Time, random = ~ 1, mixture = ~ Time, subject = "ID", ng = 4, data = data_neg_long, B = random(mod1))
mod5 <- hlme(PANSS_negativos ~ Time, random = ~ 1, mixture = ~ Time, subject = "ID", ng = 5, data = data_neg_long, B = random(mod1))
# c(mod1$BIC, mod2$BIC, mod3$BIC, mod4$BIC, mod5$BIC)
# acf(residuals(mod2), lag.max=10)


pdf(paste0(path, "/Plots_residuals_Negative_linear.pdf"), paper="a4r", pointsize = 14, width = 35, height = 26)
par(mfrow=c(2,2))
plot(residuals(mod1), main="Residuals (G=1)", xlab = "Index", ylab="Residuals")
hist(residuals(mod1), main="Residuals (G=1)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod1), main="Random effects (G=1)", xlab = "Index", ylab="Intercept")
hist(ranef(mod1), main="Random effects (G=1)")

par(mfrow=c(2,2))
plot(residuals(mod2), main="Residuals (G=2)", xlab = "Index", ylab="Residuals")
hist(residuals(mod2), main="Residuals (G=2)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod2), main="Random effects (G=2)", xlab = "Index", ylab="Intercept")
hist(ranef(mod2), main="Random effects (G=2)")

par(mfrow=c(2,2))
plot(residuals(mod3), main="Residuals (G=3)", xlab = "Index", ylab="Residuals")
hist(residuals(mod3), main="Residuals (G=3)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod3), main="Random effects (G=3)", xlab = "Index", ylab="Intercept")
hist(ranef(mod3), main="Random effects (G=3)")

par(mfrow=c(2,2))
plot(residuals(mod4), main="Residuals (G=4)", xlab = "Index", ylab="Residuals")
hist(residuals(mod4), main="Residuals (G=4)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod4), main="Random effects (G=4)", xlab = "Index", ylab="Intercept")
hist(ranef(mod4), main="Random effects (G=4)")


par(mfrow=c(2,2))
plot(residuals(mod5), main="Residuals (G=5)", xlab = "Index", ylab="Residuals")
hist(residuals(mod5), main="Residuals (G=5)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod5), main="Random effects (G=5)", xlab = "Index", ylab="Intercept")
hist(ranef(mod5), main="Random effects (G=5)")
dev.off()

#### QUADRATIC MODELS
mod1 <- hlme(PANSS_negativos ~ Time + I(Time^2), random = ~ 1, subject = "ID", ng = 1, data = data_neg_long)
mod2 <- hlme(PANSS_negativos ~ Time + I(Time^2), random = ~ 1, mixture = ~ Time + I(Time^2), subject = "ID", ng = 2, data = data_neg_long, B = random(mod1))
mod3 <- hlme(PANSS_negativos ~ Time + I(Time^2), random = ~ 1, mixture = ~ Time + I(Time^2), subject = "ID", ng = 3, data = data_neg_long, B = random(mod1))
mod4 <- hlme(PANSS_negativos ~ Time + I(Time^2), random = ~ 1, mixture = ~ Time + I(Time^2), subject = "ID", ng = 4, data = data_neg_long, B = random(mod1))
mod5 <- hlme(PANSS_negativos ~ Time + I(Time^2), random = ~ 1, mixture = ~ Time + I(Time^2), subject = "ID", ng = 5, data = data_neg_long, B = random(mod1))
# c(mod1$BIC, mod2$BIC, mod3$BIC, mod4$BIC, mod5$BIC)

pdf(paste0(path, "/Plots_residuals_Negative_quadratic.pdf"), paper="a4r", pointsize = 14, width = 35, height = 26)
par(mfrow=c(2,2))
plot(residuals(mod1), main="Residuals (G=1)", xlab = "Index", ylab="Residuals")
hist(residuals(mod1), main="Residuals (G=1)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod1), main="Random effects (G=1)", xlab = "Index", ylab="Intercept")
hist(ranef(mod1), main="Random effects (G=1)")

par(mfrow=c(2,2))
plot(residuals(mod2), main="Residuals (G=2)", xlab = "Index", ylab="Residuals")
hist(residuals(mod2), main="Residuals (G=2)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod2), main="Random effects (G=2)", xlab = "Index", ylab="Intercept")
hist(ranef(mod2), main="Random effects (G=2)")

par(mfrow=c(2,2))
plot(residuals(mod3), main="Residuals (G=3)", xlab = "Index", ylab="Residuals")
hist(residuals(mod3), main="Residuals (G=3)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod3), main="Random effects (G=3)", xlab = "Index", ylab="Intercept")
hist(ranef(mod3), main="Random effects (G=3)")

par(mfrow=c(2,2))
plot(residuals(mod4), main="Residuals (G=4)", xlab = "Index", ylab="Residuals")
hist(residuals(mod4), main="Residuals (G=4)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod4), main="Random effects (G=4)", xlab = "Index", ylab="Intercept")
hist(ranef(mod4), main="Random effects (G=4)")


par(mfrow=c(2,2))
plot(residuals(mod5), main="Residuals (G=5)", xlab = "Index", ylab="Residuals")
hist(residuals(mod5), main="Residuals (G=5)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod5), main="Random effects (G=5)", xlab = "Index", ylab="Intercept")
hist(ranef(mod5), main="Random effects (G=5)")
dev.off()

#### SPLINE MODELS
mod1 <- hlme(PANSS_negativos ~ ns(Time, knots = c(2, 6)), random = ~ 1, subject = "ID", ng = 1, data = data_neg_long)
mod2 <- hlme(PANSS_negativos ~ ns(Time, knots = c(2, 6)), random = ~ 1, mixture = ~ ns(Time, knots = c(2, 6)), subject = "ID", ng = 2, data = data_neg_long, B = random(mod1))
mod3 <- hlme(PANSS_negativos ~ ns(Time, knots = c(2, 6)), random = ~ 1, mixture = ~ ns(Time, knots = c(2, 6)), subject = "ID", ng = 3, data = data_neg_long, B = random(mod1))
mod4 <- hlme(PANSS_negativos ~ ns(Time, knots = c(2, 6)), random = ~ 1, mixture = ~ ns(Time, knots = c(2, 6)), subject = "ID", ng = 4, data = data_neg_long, B = random(mod1))
mod5 <- hlme(PANSS_negativos ~ ns(Time, knots = c(2, 6)), random = ~ 1, mixture = ~ ns(Time, knots = c(2, 6)), subject = "ID", ng = 5, data = data_neg_long, B = random(mod1))
# c(mod1$BIC, mod2$BIC, mod3$BIC, mod4$BIC, mod5$BIC)

pdf(paste0(path, "/Plots_residuals_Negative_splines.pdf"), paper="a4r", pointsize = 14, width = 35, height = 26)
par(mfrow=c(2,2))
plot(residuals(mod1), main="Residuals (G=1)", xlab = "Index", ylab="Residuals")
hist(residuals(mod1), main="Residuals (G=1)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod1), main="Random effects (G=1)", xlab = "Index", ylab="Intercept")
hist(ranef(mod1), main="Random effects (G=1)")

par(mfrow=c(2,2))
plot(residuals(mod2), main="Residuals (G=2)", xlab = "Index", ylab="Residuals")
hist(residuals(mod2), main="Residuals (G=2)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod2), main="Random effects (G=2)", xlab = "Index", ylab="Intercept")
hist(ranef(mod2), main="Random effects (G=2)")

par(mfrow=c(2,2))
plot(residuals(mod3), main="Residuals (G=3)", xlab = "Index", ylab="Residuals")
hist(residuals(mod3), main="Residuals (G=3)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod3), main="Random effects (G=3)", xlab = "Index", ylab="Intercept")
hist(ranef(mod3), main="Random effects (G=3)")

par(mfrow=c(2,2))
plot(residuals(mod4), main="Residuals (G=4)", xlab = "Index", ylab="Residuals")
hist(residuals(mod4), main="Residuals (G=4)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod4), main="Random effects (G=4)", xlab = "Index", ylab="Intercept")
hist(ranef(mod4), main="Random effects (G=4)")


par(mfrow=c(2,2))
plot(residuals(mod5), main="Residuals (G=5)", xlab = "Index", ylab="Residuals")
hist(residuals(mod5), main="Residuals (G=5)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod5), main="Random effects (G=5)", xlab = "Index", ylab="Intercept")
hist(ranef(mod5), main="Random effects (G=5)")
dev.off()

### 3.2.2. Positive PANSS ----
#### LINEAR MODELS
mod1 <- hlme(PANSS_positivos ~ Time, random = ~ 1, subject = "ID", ng = 1, data = data_pos_long)
mod2 <- hlme(PANSS_positivos ~ Time, random = ~ 1, mixture = ~ Time, subject = "ID", ng = 2, data = data_pos_long, B = random(mod1))
mod3 <- hlme(PANSS_positivos ~ Time, random = ~ 1, mixture = ~ Time, subject = "ID", ng = 3, data = data_pos_long, B = random(mod1))
mod4 <- hlme(PANSS_positivos ~ Time, random = ~ 1, mixture = ~ Time, subject = "ID", ng = 4, data = data_pos_long, B = random(mod1))
mod5 <- hlme(PANSS_positivos ~ Time, random = ~ 1, mixture = ~ Time, subject = "ID", ng = 5, data = data_pos_long, B = random(mod1))
# c(mod1$BIC, mod2$BIC, mod3$BIC, mod4$BIC, mod5$BIC)

pdf(paste0(path, "/Plots_residuals_Positive_linear.pdf"), paper="a4r", pointsize = 14, width = 35, height = 26)
par(mfrow=c(2,2))
plot(residuals(mod1), main="Residuals (G=1)", xlab = "Index", ylab="Residuals")
hist(residuals(mod1), main="Residuals (G=1)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod1), main="Random effects (G=1)", xlab = "Index", ylab="Intercept")
hist(ranef(mod1), main="Random effects (G=1)")

par(mfrow=c(2,2))
plot(residuals(mod2), main="Residuals (G=2)", xlab = "Index", ylab="Residuals")
hist(residuals(mod2), main="Residuals (G=2)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod2), main="Random effects (G=2)", xlab = "Index", ylab="Intercept")
hist(ranef(mod2), main="Random effects (G=2)")

par(mfrow=c(2,2))
plot(residuals(mod3), main="Residuals (G=3)", xlab = "Index", ylab="Residuals")
hist(residuals(mod3), main="Residuals (G=3)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod3), main="Random effects (G=3)", xlab = "Index", ylab="Intercept")
hist(ranef(mod3), main="Random effects (G=3)")

par(mfrow=c(2,2))
plot(residuals(mod4), main="Residuals (G=4)", xlab = "Index", ylab="Residuals")
hist(residuals(mod4), main="Residuals (G=4)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod4), main="Random effects (G=4)", xlab = "Index", ylab="Intercept")
hist(ranef(mod4), main="Random effects (G=4)")


par(mfrow=c(2,2))
plot(residuals(mod5), main="Residuals (G=5)", xlab = "Index", ylab="Residuals")
hist(residuals(mod5), main="Residuals (G=5)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod5), main="Random effects (G=5)", xlab = "Index", ylab="Intercept")
hist(ranef(mod5), main="Random effects (G=5)")
dev.off()

#### QUADRATIC MODELS
mod1 <- hlme(PANSS_positivos ~ Time + I(Time^2), random = ~ 1, subject = "ID", ng = 1, data = data_pos_long)
mod2 <- hlme(PANSS_positivos ~ Time + I(Time^2), random = ~ 1, mixture = ~ Time + I(Time^2), subject = "ID", ng = 2, data = data_pos_long, B = random(mod1))
mod3 <- hlme(PANSS_positivos ~ Time + I(Time^2), random = ~ 1, mixture = ~ Time + I(Time^2), subject = "ID", ng = 3, data = data_pos_long, B = random(mod1))
mod4 <- hlme(PANSS_positivos ~ Time + I(Time^2), random = ~ 1, mixture = ~ Time + I(Time^2), subject = "ID", ng = 4, data = data_pos_long, B = random(mod1))
mod5 <- hlme(PANSS_positivos ~ Time + I(Time^2), random = ~ 1, mixture = ~ Time + I(Time^2), subject = "ID", ng = 5, data = data_pos_long, B = random(mod1))
# c(mod1$BIC, mod2$BIC, mod3$BIC, mod4$BIC, mod5$BIC)

pdf(paste0(path, "/Plots_residuals_Positive_quadratic.pdf"), paper="a4r", pointsize = 14, width = 35, height = 26)
par(mfrow=c(2,2))
plot(residuals(mod1), main="Residuals (G=1)", xlab = "Index", ylab="Residuals")
hist(residuals(mod1), main="Residuals (G=1)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod1), main="Random effects (G=1)", xlab = "Index", ylab="Intercept")
hist(ranef(mod1), main="Random effects (G=1)")

par(mfrow=c(2,2))
plot(residuals(mod2), main="Residuals (G=2)", xlab = "Index", ylab="Residuals")
hist(residuals(mod2), main="Residuals (G=2)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod2), main="Random effects (G=2)", xlab = "Index", ylab="Intercept")
hist(ranef(mod2), main="Random effects (G=2)")

par(mfrow=c(2,2))
plot(residuals(mod3), main="Residuals (G=3)", xlab = "Index", ylab="Residuals")
hist(residuals(mod3), main="Residuals (G=3)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod3), main="Random effects (G=3)", xlab = "Index", ylab="Intercept")
hist(ranef(mod3), main="Random effects (G=3)")

par(mfrow=c(2,2))
plot(residuals(mod4), main="Residuals (G=4)", xlab = "Index", ylab="Residuals")
hist(residuals(mod4), main="Residuals (G=4)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod4), main="Random effects (G=4)", xlab = "Index", ylab="Intercept")
hist(ranef(mod4), main="Random effects (G=4)")


par(mfrow=c(2,2))
plot(residuals(mod5), main="Residuals (G=5)", xlab = "Index", ylab="Residuals")
hist(residuals(mod5), main="Residuals (G=5)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod5), main="Random effects (G=5)", xlab = "Index", ylab="Intercept")
hist(ranef(mod5), main="Random effects (G=5)")
dev.off()

#### SPLINE MODELS
mod1 <- hlme(PANSS_positivos ~ ns(Time, knots = c(2, 6)), random = ~ 1, subject = "ID", ng = 1, data = data_pos_long)
mod2 <- hlme(PANSS_positivos ~ ns(Time, knots = c(2, 6)), random = ~ 1, mixture = ~ ns(Time, knots = c(2, 6)), subject = "ID", ng = 2, data = data_pos_long, B = random(mod1))
mod3 <- hlme(PANSS_positivos ~ ns(Time, knots = c(2, 6)), random = ~ 1, mixture = ~ ns(Time, knots = c(2, 6)), subject = "ID", ng = 3, data = data_pos_long, B = random(mod1))
mod4 <- hlme(PANSS_positivos ~ ns(Time, knots = c(2, 6)), random = ~ 1, mixture = ~ ns(Time, knots = c(2, 6)), subject = "ID", ng = 4, data = data_pos_long, B = random(mod1))
mod5 <- hlme(PANSS_positivos ~ ns(Time, knots = c(2, 6)), random = ~ 1, mixture = ~ ns(Time, knots = c(2, 6)), subject = "ID", ng = 5, data = data_pos_long, B = random(mod1))
# c(mod1$BIC, mod2$BIC, mod3$BIC, mod4$BIC, mod5$BIC)

pdf(paste0(path, "/Plots_residuals_Positive_splines.pdf"), paper="a4r", pointsize = 14, width = 35, height = 26)
par(mfrow=c(2,2))
plot(residuals(mod1), main="Residuals (G=1)", xlab = "Index", ylab="Residuals")
hist(residuals(mod1), main="Residuals (G=1)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod1), main="Random effects (G=1)", xlab = "Index", ylab="Intercept")
hist(ranef(mod1), main="Random effects (G=1)")

par(mfrow=c(2,2))
plot(residuals(mod2), main="Residuals (G=2)", xlab = "Index", ylab="Residuals")
hist(residuals(mod2), main="Residuals (G=2)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod2), main="Random effects (G=2)", xlab = "Index", ylab="Intercept")
hist(ranef(mod2), main="Random effects (G=2)")

par(mfrow=c(2,2))
plot(residuals(mod3), main="Residuals (G=3)", xlab = "Index", ylab="Residuals")
hist(residuals(mod3), main="Residuals (G=3)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod3), main="Random effects (G=3)", xlab = "Index", ylab="Intercept")
hist(ranef(mod3), main="Random effects (G=3)")

par(mfrow=c(2,2))
plot(residuals(mod4), main="Residuals (G=4)", xlab = "Index", ylab="Residuals")
hist(residuals(mod4), main="Residuals (G=4)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod4), main="Random effects (G=4)", xlab = "Index", ylab="Intercept")
hist(ranef(mod4), main="Random effects (G=4)")


par(mfrow=c(2,2))
plot(residuals(mod5), main="Residuals (G=5)", xlab = "Index", ylab="Residuals")
hist(residuals(mod5), main="Residuals (G=5)", xlab = "Residuals", ylab="Frequency")
plot(ranef(mod5), main="Random effects (G=5)", xlab = "Index", ylab="Intercept")
hist(ranef(mod5), main="Random effects (G=5)")
dev.off()



# ---------------------------------------- #
# 4. RESULTS IN THE ORIGINAL DATASET -----
# ---------------------------------------- #

## 4.1. PANSS NEGATIUS ----
#### Clustering methods
set.seed(123)
B=1
clust.hier.results.neg <- initialize_clust_results(num_resamples = 1, k_range = 2:5)
clust.kmeans.results.neg <- initialize_clust_results(num_resamples = 1, k_range = 2:5)
clust.pam.results.neg <- initialize_clust_results(num_resamples = 1, k_range = 2:5)

D <- dist(data_neg[,-1])
clust.hier <- hclust(D, method = "ward.D2")

for (k in 2:5) {
  # 3.1. Hierarchical clustering with Ward's method
  x <- cutree(clust.hier, k)
  stats <- cluster.stats(D, x, G2 = TRUE)
  clust.hier.results.neg <- update_clust_results(clust.hier.results.neg, B, k, stats)
  
  # 3.2. K-means clustering
  x <- kmeans(D, centers = k)$cluster
  stats <- cluster.stats(D, x, G2 = TRUE)
  clust.kmeans.results.neg <- update_clust_results(clust.kmeans.results.neg, B, k, stats)
  
  # 3.3. PAM clustering
  x <- pam(D, k = k, diss = TRUE)$clustering
  stats <- cluster.stats(D, x, G2 = TRUE)
  clust.pam.results.neg <- update_clust_results(clust.pam.results.neg, B, k, stats)
}


#### LCMMs
B=1; response = "PANSS_negativos"; time= "Time"; knots = "c(2,6)"
lcmm.linear.results.neg <- initialize_lcmm_results(num_resamples =  B, k_range = 2:5)
lcmm.quadratic.results.neg <- initialize_lcmm_results(num_resamples = B, k_range = 2:5)
lcmm.splines.results.neg <- initialize_lcmm_results(num_resamples = B, k_range = 2:5)

formula.linear <- as.formula(paste(response, "~", time))
mixture.linear <- as.formula(paste("~", time))
formula.quad <- as.formula(paste(response, "~", time, "+ I(", time, "^2)"))
mixture.quad <- as.formula(paste("~", time, "+ I(", time, "^2)"))
formula.splines <- as.formula(paste(response, "~ ns(", time, ", knots =", knots, ")"))
mixture.splines <- as.formula(paste("~ ns(", time, ", knots =", knots, ")"))

set.seed(123)
lcmm.linear <- hlme(formula.linear, random=~1, subject = "ID", data = data_neg_long, ng = 1)
lcmm.quad <- hlme(formula.quad, random=~1, subject = "ID", data = data_neg_long, ng = 1)
lcmm.splines <- hlme(formula.splines, random=~1, subject = "ID", data = data_neg_long, ng = 1)


for (k in 2:5) {
  
  # 3.1. Linear
  x <- hlme(formula.linear, random = ~1, subject = "ID", data = data_neg_long,
            mixture = mixture.linear, ng = k, B = random(lcmm.linear))
  lcmm.linear.results.neg <- update_lcmm_results(lcmm.linear.results.neg, B, k, 
                                             as.data.frame(summarytable(x, which = c("G", "loglik", "AIC", "BIC", "SABIC", "entropy", "%class"), display = FALSE)))
  
  # 3.2. Quadratic
  x <- hlme(formula.quad, random = ~1, subject = "ID", data = data_neg_long,
            mixture = mixture.quad, ng = k, B = random(lcmm.quad))
  lcmm.quadratic.results.neg <- update_lcmm_results(lcmm.quadratic.results.neg, B, k,
                                                as.data.frame(summarytable(x, which = c("G", "loglik", "AIC", "BIC", "SABIC", "entropy", "%class"), display = FALSE)))
  
  
  # 3.3. Splines
  x <- hlme(formula.splines, random = ~1, subject = "ID", data = data_neg_long,
            mixture = mixture.splines, ng = k, B = random(lcmm.splines))
  lcmm.splines.results.neg <- update_lcmm_results(lcmm.splines.results.neg, B, k, 
                                              as.data.frame(summarytable(x, which = c("G", "loglik", "AIC", "BIC", "SABIC", "entropy", "%class"), display = FALSE)))
  
}


## 4.2. PANSS POSITIUS ----
#### Clustering methods
set.seed(123)
B=1
clust.hier.results.pos <- initialize_clust_results(num_resamples = 1, k_range = 2:5)
clust.kmeans.results.pos <- initialize_clust_results(num_resamples = 1, k_range = 2:5)
clust.pam.results.pos <- initialize_clust_results(num_resamples = 1, k_range = 2:5)

D <- dist(data_pos[,-1])
clust.hier <- hclust(D, method = "ward.D2")

for (k in 2:5) {
  # 3.1. Hierarchical clustering with Ward's method
  x <- cutree(clust.hier, k)
  stats <- cluster.stats(D, x, G2 = TRUE)
  clust.hier.results.pos <- update_clust_results(clust.hier.results.pos, B, k, stats)
  
  # 3.2. K-means clustering
  x <- kmeans(D, centers = k)$cluster
  stats <- cluster.stats(D, x, G2 = TRUE)
  clust.kmeans.results.pos <- update_clust_results(clust.kmeans.results.pos, B, k, stats)
  
  # 3.3. PAM clustering
  x <- pam(D, k = k, diss = TRUE)$clustering
  stats <- cluster.stats(D, x, G2 = TRUE)
  clust.pam.results.pos <- update_clust_results(clust.pam.results.pos, B, k, stats)
}

#### LCMMs
B=1; response = "PANSS_positivos"; time= "Time"; knots = "c(2,6)"
lcmm.linear.results.pos <- initialize_lcmm_results(num_resamples =  B, k_range = 2:5)
lcmm.quadratic.results.pos <- initialize_lcmm_results(num_resamples = B, k_range = 2:5)
lcmm.splines.results.pos <- initialize_lcmm_results(num_resamples = B, k_range = 2:5)

formula.linear <- as.formula(paste(response, "~", time))
mixture.linear <- as.formula(paste("~", time))
formula.quad <- as.formula(paste(response, "~", time, "+ I(", time, "^2)"))
mixture.quad <- as.formula(paste("~", time, "+ I(", time, "^2)"))
formula.splines <- as.formula(paste(response, "~ ns(", time, ", knots =", knots, ")"))
mixture.splines <- as.formula(paste("~ ns(", time, ", knots =", knots, ")"))

set.seed(123)
lcmm.linear <- hlme(formula.linear, random=~1, subject = "ID", data = data_pos_long, ng = 1)
lcmm.quad <- hlme(formula.quad, random=~1, subject = "ID", data = data_pos_long, ng = 1)
lcmm.splines <- hlme(formula.splines, random=~1, subject = "ID", data = data_pos_long, ng = 1)


for (k in 2:5) {
  
  # 3.1. Linear
  x <- hlme(formula.linear, random = ~1, subject = "ID", data = data_pos_long,
            mixture = mixture.linear, ng = k, B = random(lcmm.linear))
  lcmm.linear.results.pos <- update_lcmm_results(lcmm.linear.results.pos, B, k, 
                                                 as.data.frame(summarytable(x, which = c("G", "loglik", "AIC", "BIC", "SABIC", "entropy", "%class"), display = FALSE)))
  
  # 3.2. Quadratic
  x <- hlme(formula.quad, random = ~1, subject = "ID", data = data_pos_long,
            mixture = mixture.quad, ng = k, B = random(lcmm.quad))
  lcmm.quadratic.results.pos <- update_lcmm_results(lcmm.quadratic.results.pos, B, k,
                                                    as.data.frame(summarytable(x, which = c("G", "loglik", "AIC", "BIC", "SABIC", "entropy", "%class"), display = FALSE)))
  
  
  # 3.3. Splines
  x <- hlme(formula.splines, random = ~1, subject = "ID", data = data_pos_long,
            mixture = mixture.splines, ng = k, B = random(lcmm.splines))
  lcmm.splines.results.pos <- update_lcmm_results(lcmm.splines.results.pos, B, k, 
                                                  as.data.frame(summarytable(x, which = c("G", "loglik", "AIC", "BIC", "SABIC", "entropy", "%class"), display = FALSE)))
  
}


## 4.3. Make final table ----
clust.hier.results.neg <- clust.hier.results.neg %>%
  mutate(Method = "Hierarchical", Type = "N-PANSS") %>%
  select(Method, Type, everything())
clust.kmeans.results.neg <- clust.kmeans.results.neg %>%
  mutate(Method = "K-means", Type = "N-PANSS") %>%
  select(Method, Type, everything())
clust.pam.results.neg <- clust.pam.results.neg %>% 
  mutate(Method = "PAM", Type = "N-PANSS") %>%
  select(Method, Type, everything())
clust.results.neg <- bind_rows(clust.hier.results.neg, clust.kmeans.results.neg, clust.pam.results.neg) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

lcmm.linear.results.neg <- lcmm.linear.results.neg %>%
  mutate(Method = "LCMM Linear", Type = "N-PANSS") %>%
  select(Method, Type, everything())
lcmm.quadratic.results.neg <- lcmm.quadratic.results.neg %>%
  mutate(Method = "LCMM Quadratic", Type = "N-PANSS") %>%
  select(Method, Type, everything())
lcmm.splines.results.neg <- lcmm.splines.results.neg %>%
  mutate(Method = "LCMM Splines", Type = "N-PANSS") %>%
  select(Method, Type, everything())
lcmm.results.neg <- bind_rows(lcmm.linear.results.neg, lcmm.quadratic.results.neg, lcmm.splines.results.neg) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

clust.hier.results.pos <- clust.hier.results.pos %>%
  mutate(Method = "Hierarchical", Type = "P-PANSS") %>%
  select(Method, Type, everything())
clust.kmeans.results.pos <- clust.kmeans.results.pos %>%
  mutate(Method = "K-means", Type = "P-PANSS") %>%
  select(Method, Type, everything())
clust.pam.results.pos <- clust.pam.results.pos %>% 
  mutate(Method = "PAM", Type = "P-PANSS") %>%
  select(Method, Type, everything())
clust.results.pos <- bind_rows(clust.hier.results.pos, clust.kmeans.results.pos, clust.pam.results.pos) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

lcmm.linear.results.pos <- lcmm.linear.results.pos %>%
  mutate(Method = "LCMM Linear", Type = "P-PANSS") %>%
  select(Method, Type, everything())
lcmm.quadratic.results.pos <- lcmm.quadratic.results.pos %>%
  mutate(Method = "LCMM Quadratic", Type = "P-PANSS") %>%
  select(Method, Type, everything())
lcmm.splines.results.pos <- lcmm.splines.results.pos %>%
  mutate(Method = "LCMM Splines", Type = "P-PANSS") %>%
  select(Method, Type, everything())
lcmm.results.pos <- bind_rows(lcmm.linear.results.pos, lcmm.quadratic.results.pos, lcmm.splines.results.pos) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
  
# Export resuts of metrics summary into an Excel file
wb = createWorkbook()
addWorksheet(wb, "Clust_neg"); writeData(wb, "Clust_neg", clust.results.neg)
addWorksheet(wb, "LCMM_neg"); writeData(wb, "LCMM_neg", lcmm.results.neg)
addWorksheet(wb, "Clust_pos"); writeData(wb, "Clust_pos", clust.results.pos)
addWorksheet(wb, "LCMM_pos"); writeData(wb, "LCMM_pos", lcmm.results.pos)
saveWorkbook(wb, "C:/Users/laurajulia/OneDrive - Universitat de Barcelona/Escritorio/IDIBAPS/farmaPRED/Comparison_Methods/1. Real dataset/figures/Metrics_real_data.xlsx", overwrite = TRUE)

# ---------------------------- #
# 5. BOOTSTRAP RESAMPLES -----
# ---------------------------- #

## 5.1. NEGATIUS amb N= 500, k=2:5 ----
### Clustering techniques
start_time <- Sys.time()
set.seed(123)
results.clust_neg <- perform_clustering_analysis(data_neg, 
                                                 num_resamples = 500, 
                                                 k_range = 2:5,
                                                 min.percentage = 5)
end_time <- Sys.time()
time_clust_neg <- end_time - start_time


### LCMMs
start_time <- Sys.time()
set.seed(123)
results.lcmm_neg  <- perform_lcmm_analysis(data= data_neg_long, 
                                           response= "PANSS_negativos",
                                           time= "Time", 
                                           knots = "c(2,6)", 
                                           num_resamples = 500,
                                           k_range = 2:5, 
                                           min.percentage = 5)
end_time <- Sys.time()
time_lcmm_neg <- end_time - start_time

save.image("Results_bootstraps_180425_neg.RData")
#load("Results_bootstraps_180425_neg.RData")

## 5.2. POSITIUS amb N= 500, k=2:5 ----
### Clustering techniques
start_time <- Sys.time()
set.seed(123)
results.clust_pos <- perform_clustering_analysis(data_pos,
                                                 num_resamples = 500,
                                                 k_range = 2:5,
                                                 min.percentage = 5)
end_time <- Sys.time()
time_clust_pos <- end_time - start_time

### LCMMs
start_time <- Sys.time()
set.seed(123)
results.lcmm_pos  <- perform_lcmm_analysis(data= data_pos_long, 
                                           response= "PANSS_positivos",
                                           time= "Time", 
                                           knots = "c(2,6)", 
                                           num_resamples = 500,
                                           k_range = 2:5, 
                                           min.percentage = 5)
end_time <- Sys.time()
time_lcmm_pos <- end_time - start_time




save.image("Results_bootstraps_300525.RData")

# ------------------- #
# 6. INDEX ICC  -----
# ------------------- #

## 6.1. PANSS Negativos ----
data_icc_neg <- data_neg %>% rename(id = Ident_caso,
                                    PANSS_0 = PANSS_negativos_VB,
                                    PANSS_2 = PANSS_negativos_V2M,
                                    PANSS_6 = PANSS_negativos_V6M,
                                    PANSS_12 = PANSS_negativos_V12M)

set.seed(123)
plot_neg_2 <- plot_icc_clusters(data_icc_neg, k = 2)
plot_neg_3 <- plot_icc_clusters(data_icc_neg, k = 3)
plot_neg_4 <- plot_icc_clusters(data_icc_neg, k = 4)
plot_neg_5 <- plot_icc_clusters(data_icc_neg, k = 5)

## 6.2. PANSS Positivos ----
data_icc_pos <- data_pos %>% rename(id = Ident_caso,
                                    PANSS_0 = PANSS_positivos_VB,
                                    PANSS_2 = PANSS_positivos_V2M,
                                    PANSS_6 = PANSS_positivos_V6M,
                                    PANSS_12 = PANSS_positivos_V12M)

set.seed(123)
plot_pos_2 <- plot_icc_clusters(data_icc_pos, k = 2)
plot_pos_3 <- plot_icc_clusters(data_icc_pos, k = 3)
plot_pos_4 <- plot_icc_clusters(data_icc_pos, k = 4)
plot_pos_5 <- plot_icc_clusters(data_icc_pos, k = 5)

