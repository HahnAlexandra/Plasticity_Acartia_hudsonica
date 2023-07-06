##########################################################
##                        Stats                         ##
##                   Alexandra Hahn                     ##
##########################################################


setwd("~/Documents/Scripts/Thesis")

library(tidyverse)
library(olsrr)
library(car)
library(emmeans)
library(multcomp)
library(multcompView)

#import full data set
assays <- read.csv("~/Documents/Scripts/Thesis/data.csv", sep=";", header = TRUE)
assays <- assays[-c(1, 59, 136, 246), ]#remove dead or inactive
assays <-  assays[!(is.na(assays$collection)),]#remove potential NA values created when importing


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
  
  wild <- assays[which(assays$treatment == "wild"),]
  wild2 <- wild[which(wild$species != "tonsa"),]
  
} 



####Figure 1####
#models for Ctmax and length
mF1 <- lm(Ctmax ~ collection, data = wild2)
summary(mF1); Anova(mF1, type = 3)

mF1b <- lm(length ~ collection, data = wild2)
summary(mF1b); Anova(mF1b, type = 3)

#model means and compact letters, post-hoc test
#Ctmax
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


####Figure 2####
#for Ctmax
hudsonica %>% 
  group_by(collection:treatment) %>% 
  dplyr::summarise(mean = mean(Ctmax), sd = sd(Ctmax))#means and sd

#model
intF2 <- interaction(hudsonica$treatment, hudsonica$collection)
mF2 <- lm(Ctmax ~ intF2, data = data)
summary(mF2)
Anova(mF2, type = 3)
boxplot(Ctmax ~intF2, data = hudsonica)

#compact letters
mmF2 <- emmeans(object = mF2,
                specs = "intF5")#gets adjusted and weighted means per group

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
                 specs = "intF5")#gets adjusted and weighted means per group

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

#### Figure 3####
rnorm <- assays[which(assays$treatment != "wild" & assays$species != "tonsa"),]
rnorm_stat  <- rnorm[which(rnorm$collection != "5"),]#excluding col-5 bc it has no warm data

#slopes for Ctmax
models <- lapply(unique(rnorm_stat$collection), function(col) {
  lm(Ctmax ~ treatment, data = subset(rnorm_stat, collection == col))
})#fit lin regression models

slopes <- sapply(models, function(model) {
  coef(model)[2]
})#extract slopes

slopes# col-1 0.9723269, col-2 2.4471053, col-4 0.9976983

#slopes for length
models2 <- lapply(unique(rnorm_stat$collection), function(col) {
  lm(length ~ treatment, data = subset(rnorm_stat, collection == col))
})

# Extract slopes
slopes2 <- sapply(models2, function(model) {
  coef(model)[2]
})

slopes2# col-1 -82.07544, col-2 -40.3333, col-4 -70.17644

#stats for reaction norms Ctmax
m_rnorm_c <- lm(Ctmax ~ treatment*collection, rnorm)
summary(m_rnorm_c)
Anova(m_rnorm_c, type = 3, singular.ok = TRUE)

hist(resid(m_rnorm_c))
ols_test_normality(m_rnorm_c)#Kolmogorov-Smirnov for observations larger than 50, if p-value is > 0.05 --> residuals normally distributed
ols_test_correlation(m_rnorm_c)#correlation between observed and expected residuals (under normality)

#post-hoc for Ctmax
pairs(emmeans(m_rnorm_c, ~treatment|collection))#all cold-warm sig different
pairs(emmeans(m_rnorm_c, ~collection|treatment))#col 2 sig different at cold


#stats for reaction norms length
m_rnorm_l <- lm(length ~ treatment*collection, rnorm)
summary(m_rnorm_l)
Anova(m_rnorm_l, type = 3, singular.ok = TRUE)

hist(resid(m_rnorm_l))
ols_test_normality(m_rnorm_l)
ols_test_correlation(m_rnorm_l)

#post-hoc for length
pairs(emmeans(m_rnorm_l, ~treatment|collection))#all cold-warm sif different
pairs(emmeans(m_rnorm_l, ~collection|treatment))

#drop col 2 then model for Ctmax
m_rnorm_c2 <- rnorm %>%
  filter(collection != "2")%>%
  lm(Ctmax ~ treatment*collection, .)

Anova(m_rnorm_c2, type = 3, singular.ok = TRUE)# no interaction
pairs(emmeans(m_rnorm_c2, ~collection|treatment))# marginally significant

#drop col 2 then model for length
m_rnorm_l2 <- rnorm %>%
  filter(collection != "2")%>%
  lm(length ~ treatment*collection, .)

Anova(m_rnorm_l2, type = 3, singular.ok = TRUE)#no interaction
pairs(emmeans(m_rnorm_l2, ~collection|treatment))#significantly different

#### Figure 4 ####




