##########################################################
##                        Stats                         ##
##                   Alexandra Hahn                     ##
##########################################################

#this script produces the statistical analysis for
#"Phenotypic Plasticity Drives Seasonal Thermal Tolerance in a Baltic Copepod"

#setwd("~/your_path")

library(tidyverse)
library(olsrr)
library(car)
library(emmeans)
library(multcomp)
library(multcompView)
library(lme4)
library(sjPlot)
library(gridExtra)

#import full data set
assays <- read.csv("data.csv", sep=";", header = TRUE)
assays <- assays[-c(1, 59, 136, 246), ]#remove dead or inactive
assays <-  assays[!(is.na(assays$collection)),]#remove potential NA values created when importing

####data preparation####

#run for prep 
{
  
  # Adding columns based on tank side and mean length:
  assays<- assays %>%
    mutate(tank_side = case_when(
      grepl("R", assays$ID) ~ "right",
      grepl("L", assays$ID) ~ "left"
    ))
  
  assays$length <- (assays$length1+assays$length2+assays$length3)/3  #new column mean length in µm
  
  # set as factor and order levels
  assays$collection <- factor(assays$collection, levels =c("1", "2", "3", "4", "5"))
  assays$treatment <- factor(assays$treatment, levels = c("wild", "cold", "warm"))
  assays$generation <- factor(assays$generation, levels = c("parental", "f1"))
  assays$sex_confirmed <- factor(assays$sex_confirmed, levels = c("f", "m"))
  assays$tank_side <- factor(assays$tank_side, levels = c("left", "right"))
  assays$position <- factor(assays$position, levels = c("side", "top"))
  assays$mean2 <- as.factor(assays$X2.week_mean)
  
  #subsets
  hudsonica <- assays %>%
    filter(species == "hudsonica")%>%
    filter(collection != 3)
  
  hudsonica.2 <- hudsonica[which(hudsonica$collection != 2 | hudsonica$treatment != "cold"),]
  
  wild <- assays[which(assays$treatment == "wild"),]
  wild2 <- wild[which(wild$species != "tonsa"),]
  f1 <- hudsonica[which(hudsonica$generation == "f1"),]
  f1.2 <- f1[which(f1$collection != 2 | f1$treatment != "cold"),]
} 


####CTmax models####

#full CTmax model with random effects
m0 <- lmer(Ctmax~treatment*collection* sex_confirmed*length + (1|vial_number)+ (1|time_assay) + (1|tank_side),
           data = hudsonica)
summary(m0)
Anova(m0, type = 3)

# this visualizes the "random effects"
p <- plot_model(m0, type = "re", facet.grid=FALSE) 
grid.arrange(p[[1]], p[[2]], p[[3]])# no visible strong effects....

#random effects do not improve model fit --> drop random effects

#full CTmax model without random effects
m0_1 <- lm(Ctmax~treatment*collection*sex_confirmed*length, data = hudsonica)
summary(m0_1)
Anova(m0_1, singular.ok = T, type = 3)

#biological CTmax model after model testing and comparing biological relevance
m <- lm(Ctmax~collection*treatment + sex_confirmed + length, data = hudsonica)
summary(m)
Anova(m, singular.ok = T, type = 3)

#model comparison
AIC(m, m0, m0_1)
BIC(m, m0, m0_1)
#biological model is chosen, lowest AIC and BIC

#diagnostics for biological Ctmax model (m)

hist(resid(m))#normal distribution
ols_test_normality(m)#Kolmogorov-Smirnov for observation > 50, not sig = normal distribution
ols_test_correlation(m)#correlation between observed and expected residuals (under normality)

par(mfrow = c(2,2))
plot(m)
par(mfrow = c(1,1))

durbinWatsonTest(m)#no autocorrelation if p > 0.05

#post-hoc pairwise comparison CTmax
m_emmeans <- emmeans(m, pairwise ~ collection:treatment)
summary(m_emmeans$emmeans)
summary(m_emmeans$contrasts)

#model substituting collection with collection temperature
m2 <- lm(Ctmax~treatment*tempf1_parental + sex_confirmed + length, data = hudsonica)
summary(m2)#with treatment and mean temperature same information is conveyed for common garden
Anova(m2, singular.ok = T, type = 3)

hist(resid(m2))

par(mfrow = c(2,2))
plot(m2)
par(mfrow = c(1,1))

AIC(m, m2)#m has lower AIC
#statistical output very similar for both models

#Ctmax models excluding cold col-2
m.2 <- lm(Ctmax~collection*treatment + sex_confirmed + length, data = hudsonica.2)
summary(m.2)
Anova(m.2, singular.ok = T, type = 3)

hist(resid(m.2))

par(mfrow = c(2,2))
plot(m.2)
par(mfrow = c(1,1))

m2.2 <- lm(Ctmax~treatment*tempf1_parental + sex_confirmed + length, data = hudsonica.2)
summary(m2.2)#with treatment and mean collection temperature same information is conveyed 
Anova(m2.2, singular.ok = T, type = 3)

hist(resid(m2.2))

par(mfrow = c(2,2))
plot(m2.2)
par(mfrow = c(1,1))


####Length models####

#full CTmax model with random effects
ml0 <- lmer(length~treatment*collection* sex_confirmed + (1|vial_number)+ (1|time_assay) + (1|tank_side),
           data = hudsonica)
summary(ml0)
Anova(ml0, type = 3)

# this visualizes the "random effects"
p <- plot_model(ml0, type = "re", facet.grid=FALSE) 
grid.arrange(p[[1]], p[[2]], p[[3]])# no visible strong effects....

#random effects do not improve model fit --> drop random effects

#full length model without random effects
ml0_1 <- lm(length~treatment*collection*sex_confirmed, data = hudsonica)
summary(ml0_1)
Anova(ml0_1, singular.ok = T, type = 3)

#biological CTmax model after model testing and comparing biological relevance
ml <- lm(length~collection*treatment + sex_confirmed, data = hudsonica)
summary(ml)
Anova(ml, singular.ok = T, type = 3)

#model comparison
AIC(ml, ml0, ml0_1)
BIC(ml, ml0, ml0_1)
#biological model is chosen, lowest BIC, full model has lower AIC but tradeoff with df

#diagnostics for biological Ctmax model (m)

hist(resid(ml))#normal distribution
ols_test_normality(ml)#Kolmogorov-Smirnov for observation > 50, not sig = normal distribution
ols_test_correlation(ml)#correlation between observed and expected residuals (under normality)

par(mfrow = c(2,2))
plot(m)
par(mfrow = c(1,1))

durbinWatsonTest(ml)#no autocorrelation if p > 0.05

#post-hoc pairwise comparison CTmax
ml_emmeans <- emmeans(ml, pairwise ~ collection:treatment)
summary(ml_emmeans$emmeans)
summary(ml_emmeans$contrasts)

#model substituting collection with collection temperature
ml2 <- lm(length~treatment*tempf1_parental + sex_confirmed, data = hudsonica)
summary(ml2)#with treatment and mean temperature same information is conveyed for common garden
Anova(ml2, singular.ok = T, type = 3)

hist(resid(ml2))

par(mfrow = c(2,2))
plot(ml2)
par(mfrow = c(1,1))

AIC(m, m2)#m has lower AIC
#statistical output very similar for both models

#Ctmax models excluding cold col-2
ml.2 <- lm(length~collection*treatment + sex_confirmed, data = hudsonica.2)
summary(ml.2)
Anova(ml.2, singular.ok = T, type = 3)

hist(resid(ml.2))

par(mfrow = c(2,2))
plot(ml.2)
par(mfrow = c(1,1))

ml2.2 <- lm(length~treatment*tempf1_parental + sex_confirmed + length, data = hudsonica.2)
summary(ml2.2)#with treatment and mean collection temperature same information is conveyed 
Anova(ml2.2, singular.ok = T, type = 3)

hist(resid(ml2.2))

par(mfrow = c(2,2))
plot(ml2.2)
par(mfrow = c(1,1))

#####F1 models#####

#models for Ctmax F1
m_f1 <- lm(Ctmax~treatment * tempf1_parental + sex_confirmed + length, data = f1)
summary(m_f1)#1.5 degree between cold and warm
Anova(m_f1, singular.ok = T, type = 3)#parental effect driven by cold col-2

hist(resid(m_f1))
par(mfrow = c(2,2))
plot(m_f1)
par(mfrow = c(1,1))

#build table for reaction norm usind model means
means_Ctmax <- emmeans(m_f1, ~ treatment * collection)
summary_means_Ctmax <- summary(means_Ctmax)
summary_means_Ctmax <-  summary_means_length[!(is.na(summary_means_length$emmean)),]

#remove cold col-2
m_f1.2 <- lm(Ctmax~treatment * tempf1_parental + sex_confirmed + length, data = f1.2)
summary(m_f1.2)#0.7 degree between cold and warm
Anova(m_f1.2, singular.ok = T, type = 3)# no parental effect

hist(resid(m_f1.2))#removal greatly improves distribution of residuals
par(mfrow = c(2,2))
plot(m_f1.2)
par(mfrow = c(1,1))

#models for length F1
m_lf1 <- lm(length~treatment * tempf1_parental + sex_confirmed, data = f1)
summary(m_lf1)#-57 microM between cold-warm
Anova(m_lf1, singular.ok = T, type = 3)#parental effect

hist(resid(m_lf1))
par(mfrow = c(2,2))
plot(m_lf1)
par(mfrow = c(1,1))

#build table for reaction norm usind model means
means_length <- emmeans(m_lf1, ~ treatment * collection)
summary_means_length <- summary(means_length)
summary_means_length <-  summary_means_length[!(is.na(summary_means_length$emmean)),]

#remove cold col-2
m_lf1.2 <- lm(length~treatment * tempf1_parental + sex_confirmed, data = f1.2)
summary(m_lf1.2)#-73 microM between cold-warm
Anova(m_lf1.2, singular.ok = T, type = 3)# still parental effect

hist(resid(m_lf1.2))# removal greatly improves distribution of residuals
par(mfrow = c(2,2))
plot(m_lf1.2)
par(mfrow = c(1,1))


####Figure 2 - compact letters####
#models for Ctmax and length
mF1 <- lm(Ctmax ~ collection, data = wild2)
summary(mF1); Anova(mF1, type = 3)

mF1b <- lm(length ~ collection, data = wild2)
summary(mF1b); Anova(mF1b, type = 3)

#model means and compact letters, post-hoc test
#Ctmax
cor.test(wild2$X2.week_mean, wild2$Ctmax, method = "pearson")#spearman correlation 0.766

wild2 %>% 
  group_by(collection) %>% 
  dplyr::summarise(mean = mean(Ctmax), sd = sd(Ctmax))#gives means and sd

pairwise.t.test(wild2$Ctmax, wild2$collection, p.adjust ="holm")

mmF1 <- emmeans(object = mF1,
                specs = "collection")#gets adjusted and weighted means per group

mmF1_cld <- cld(object = mmF1,
                adjust = "holm",
                Letters = letters,
                alpha = 0.05)#adds compact letters

mmF1_cld

#length
cor.test(wild2$X2.week_mean, wild2$length, method = "pearson")#spearman correlation -0.529

wild2 %>% 
  filter(is.na(wild2$length) == FALSE) %>%
  group_by(collection) %>% 
  dplyr::summarise(mean = mean(length), sd = sd(length))#gives means and sd

pairwise.t.test(wild$length, wild$collection, p.adjust ="holm")

mmF1b <- emmeans(object = mF1b,
                 specs = "collection")#gets adjusted and weighted means per group

mmF1b_cld <- cld(object = mmF1b,
                 adjust = "holm",
                 Letters = letters,
                 alpha = 0.05)#adds compact letters

mmF1b_cld


####Figure 3 - compact letters####
#for Ctmax
hudsonica %>% 
  group_by(collection:treatment) %>% 
  dplyr::summarise(mean = mean(Ctmax), sd = sd(Ctmax))#means and sd

#model
intF2 <- interaction(hudsonica$treatment, hudsonica$collection)
mF2 <- lm(Ctmax ~ intF2, data = hudsonica)
summary(mF2)
Anova(mF2, type = 3)
boxplot(Ctmax ~intF2, data = hudsonica)

#compact letters
mmF2 <- emmeans(object = mF2,
                specs = "intF2")#gets adjusted and weighted means per group

mmF2_cld <- cld(object = mmF2,
                adjust = "holm",
                Letters = letters,
                alpha = 0.05)#adds compact letters

mmF2_cld#warm ones are together, cold ones are together with wild 2

hist(resid(mF2))
ols_test_normality(mF2)
ols_test_correlation(mF2)#correlation between observed and expected residuals (under normality)

par(mfrow = c(2,2))
plot(mF2)
par(mfrow = c(1,1))

durbinWatsonTest(mF2)#no autocorrelation if p > 0.05

#post hoc test
pairwise.t.test(hudsonica$Ctmax, hudsonica$collection:hudsonica$treatment, p.adjust ="holm")

summary(glht(model= mF2, linfct= mcp(intF2 = "Tukey")))

ph_F2 <- contrast(mmF2, method = "tukey")#pairwise comparison
summary(ph_F2)$p.value#summary of p-values

#for length
hudsonica %>% 
  filter( is.na(hudsonica$length) == FALSE) %>%
  group_by(collection:treatment) %>% 
  dplyr::summarise(mean = mean(length), sd = sd(length))

mF2b <- lm(length ~ intF2, data = hudsonica)
summary(mF2b)
Anova(mF2b, type = 3)
boxplot(length ~intF2, data = hudsonica)


mmF2b <- emmeans(object = mF2b,
                 specs = "intF2")#gets adjusted and weighted means per group

mmF2b_cld <- cld(object = mmF2b,
                 adjust = "holm",
                 Letters = letters,
                 alpha = 0.05)#adds compact letters

mmF2b_cld#warm ones are together, cold ones are together with wild 2

hist(resid(mF2b))
ols_test_normality(mF2b)
ols_test_correlation(mF2b)#correlation between observed and expected residuals (under normality)

par(mfrow = c(2,2))
plot(mF2b)
par(mfrow = c(1,1))

durbinWatsonTest(mF2b)#no autocorrelation if p > 0.05

#post hoc test
pairwise.t.test(hudsonica$length, hudsonica$collection:hudsonica$treatment, p.adjust ="holm")

summary(glht(model= mF2b, linfct= mcp(intF2 = "Tukey")))

ph_F2b <- contrast(mmF2b, method = "tukey")#pairwise comparison
pval_F2b <- summary(ph_F2b)$p.value#summary of p-values







