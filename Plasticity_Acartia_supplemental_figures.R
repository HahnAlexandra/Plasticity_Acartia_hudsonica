##########################################################
##                 Supplemental Figures                 ##
##                   Alexandra Hahn                     ##
##########################################################

#supplemental figures for Seasonal plasticity in Acartia hudsonica

setwd("~/Documents/Scripts/Thesis")

#load necessary packages
library(tidyverse)
library(wesanderson)
library(cowplot)
library(utils)


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
  data <- assays %>%
    filter(species != "tonsa")
  
  hudsonica <- assays %>%
    filter(species == "hudsonica")%>%
    filter(collection != 3)
  
  all <- assays %>%
    filter(collection != 3)
  
  wild <- assays[which(assays$treatment == "wild"),]
  wild2 <- wild[which(wild$species != "tonsa"),]
  
} 

#### Figure S.1 - heating rate #####

temp_r <- read.csv("temp_running.csv", sep = ";")#imports experimental thermometer


#only consider measurements each minute
temp_1m <- temp_r[grepl(glob2rx('*:*:00'), temp_r$time_running),]

#loop to create to data sets, before and after 22.5 °C (point of unplugging the second heater)
i = 1

lowout <- c()
highout <- c()
for(i in 1:max(temp_1m$trial, na.rm=T)){
  temp1 <-temp_1m[which(temp_1m$trial == i),]
  temp1$difference <- c(0,diff(temp1$temperature))#adapt 
  low1 <- temp1[which(temp1$temperature < 22.5),]
  high1 <- temp1[which(temp1$temperature > 22.5),]
  lowout[i]<- mean(low1$difference)
  highout[i]<- mean(high1$difference)
  
}

dev.new()#new window for plots
par(mfrow = c(1,2))

#histogram of heating rates before unplugging
hist(lowout, xlab = "mean heating rate in °C/min", breaks = 20, main = "before 22.5 °C")
abline(v = mean(lowout),                       # Add line for mean
       col = "red",
       lwd = 3)
text(x = 0.155,                   
     y = 10,
     paste("Mean =",0.207),
     col = "red",
     cex = 1)
abline(v = median(lowout),                    # Add line for median
       col = "blue",
       lwd = 3)
text(x = 0.155,                 
     y = 9,
     paste("Median =",0.210),
     col = "blue",
     cex = 1)

#histogram of heating rates after unplugging
hist(highout, xlab = "mean heating rate in °C/min", breaks = 20, main = "after 22.5 °C")
abline(v = mean(highout),                       # Add line for mean
       col = "red",
       lwd = 3)
text(x = 0.1,                   
     y = 8,
     paste("Mean =", 0.114),
     col = "red",
     cex = 1)
abline(v = median(highout),                    # Add line for median
       col = "blue",
       lwd = 3)
text(x = 0.1,                 
     y = 7,
     paste("Median =", 0.117),
     col = "blue",
     cex = 1)


#### Figure S.2 - developmental temperatures ####

# plot for daily mean temperature
group_mean <- interaction(data$daily_mean,data$sex_confirmed)#with daily mean kimocc
S2_A <- ggplot(data,aes(x=daily_mean, 
                            y=Ctmax,
                            color=sex_confirmed,
                            fill = sex_confirmed,
                            group = group_mean)) +
  theme_light(base_size = 14)+
  ylab(expression("CT"["max"]* " in °C"))+
  xlab("Developmental temperature in °C")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  scale_fill_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  geom_point(aes(shape = generation), size = 2.75)+
  geom_smooth(method = "lm", aes(group = sex_confirmed, col = sex_confirmed))


# plot for 2 week mean temperature
group_2w <- interaction(data$X2.week_mean,data$sex_confirmed)#with 2-week mean
S2_B <- ggplot(data,aes(x=X2.week_mean, 
                          y=Ctmax,
                          color=sex_confirmed,
                          fill = sex_confirmed,
                          group = group_2w)) +
  theme_light(base_size = 14)+
  ylab(expression("CT"["max"]* " in °C"))+
  xlab("Developmental temperature in °C")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  scale_fill_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  geom_point(aes(shape = generation), size = 2.75)+
  geom_smooth(method = "lm", aes(group = sex_confirmed, col = sex_confirmed))

#combine all plots
legendS2 <- get_legend(ggplot(data,aes(x=X2.week_mean, 
                        y=Ctmax,
                        color=sex_confirmed,
                        fill = sex_confirmed,
                        group = group_2w)) +
  theme_light(base_size = 14)+
  ylab(expression("CT"["max"]* " in °C"))+
  xlab("Developmental temperature in °C")+
  scale_color_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  scale_fill_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  geom_point(aes(shape = generation), size = 2.75)+
  geom_smooth(method = "lm", aes(group = sex_confirmed, col = sex_confirmed)))

#combine plots and plot in new window
dev.new()
combined_S2 <- plot_grid(S2_A, S2_B, ncol = 1, labels = c("A", "B"))
plot_grid(nrow = 1, combined_S2, legendS2, ncol = 2 , rel_widths = c(3/4, 1/4))

cor.test(data$daily_mean, data$Ctmax, method = "spearman")
cor.test(data$X2.week_mean, data$Ctmax, method = "spearman")

#### Figure S.3 - corr SST and Ctmax ####

ggplot(wild2, aes(y = Ctmax, x = X2.week_mean, col = X2.week_mean))+
  geom_point(size = 2.75)+
  geom_smooth(method = "lm", color = "grey", fill = "lightgrey")+
  scale_fill_gradientn(colors = alpha(wes_palette("Zissou1", type = "continuous")))+
  scale_color_gradientn(colors = alpha(wes_palette("Zissou1", type = "continuous"),0.5))+
  theme_light(base_size = 14)+
  theme(legend.position = "none")+
  ylab("Critical thermal maximum in °C")+
  xlab("Mean SST in °C")

cor.test(wild2$X2.week_mean, wild2$Ctmax, method = "spearman")
mS3 <- lm(Ctmax ~ X2.week_mean, data = wild2)
summary(mS3)

#### Figure S.4 - effects of developmental temperature ####

#Ctmax
S4_A <- ggplot(hudsonica, aes(y = Ctmax, x = X2.week_mean, col = sex_confirmed))+
  geom_point(aes(shape = generation), size = 2.75)+
  geom_smooth(method = "lm", aes(group = sex_confirmed, col = sex_confirmed, fill = sex_confirmed))+
  scale_color_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  scale_fill_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  scale_x_continuous()+
  theme_light(base_size = 14)+
  theme(legend.position = "none")+
  xlab("Developmental temperature in °C")+
  ylab(expression("CT"["max"]* " in °C"))

mS4_A <- lm(Ctmax ~ X2.week_mean + sex_confirmed, data = hudsonica)
summary(mS4_A)
Anova(mS4_A, type = 3)

#length
S4_B <- ggplot(hudsonica, aes(y = length, x = X2.week_mean, col = sex_confirmed))+
  geom_point(aes(shape = generation), size = 2.75)+
  geom_smooth(method = "lm", aes(group = sex_confirmed, col = sex_confirmed, fill = sex_confirmed))+
  scale_color_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  scale_fill_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  scale_x_continuous()+
  theme_light(base_size = 14)+
  theme(legend.position = "none")+
  xlab("Developmental temperature in °C")+
  ylab("Prosome length in µm")

mS4_B <- lm(length ~ X2.week_mean + sex_confirmed, data = hudsonica)
summary(mS4_B)
Anova(mS4_B, type = 3)

#combine plots and plot in new window
dev.new()
combined_S4 <- plot_grid(S4_A, S4_B, ncol = 1, labels = c("A", "B"))
plot_grid(nrow = 1, combined_S4, legendS2, ncol = 2 , rel_widths = c(3/4, 1/4))


#### Figure S.5 - Ctmax all ####
#CTmax - all collections
all2 <- distinct(all, date_sampling, treatment)%>%
  arrange(date_sampling, treatment)
all2$yloc <- max(all$Ctmax)+ .5
all2$label <- c("a", "b", "cd", "bc", "a", "cd", "de", "bc", "f", "ef", "bc", "g")#compact letter from Ctmax_stats

Ctmax_all <- ggplot(all, aes(x=date_sampling, y=Ctmax, fill=treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2)+
  theme_light(base_size = 11)+
  scale_x_discrete(labels=c("2022-04-06" = "April", "2022-05-16" = " May", "2022-06-27" = " June" , "2022-07-19" = "July"))+
  scale_y_continuous(limits = c(23,36), n.breaks = 5)+
  xlab("Sampling time")+ ylab("Critical thermal maximum in °C")+
  theme(legend.position = "right")+
  scale_fill_manual(values = c("#96B48E","#CFD4EB","#D5968F"), name = "Treatment")

Ctmax_all +
  ylim(NA, max(all$Ctmax)+ .5)+
  geom_text(data = all2, aes(y = yloc, label = label),
            position = position_dodge(width = .75))

#### Figure S.6  ####
#https://github.com/HahnAlexandra/Plasticity_Acartia_hudsonica/tree/main/genotyping
#find pipeline to create plot in directory above

#### Figure S.7 - species identity####
S7_A <- ggplot(all, aes(x = Ctmax, y = length, col = species_gen))+
  geom_point(aes(shape = generation), size = 2.75)+
  scale_color_manual(values = c("#333333","#CFD4EB", "#96B48E","#D5968F"), name = "species")+
  theme_light(base_size = 14)+
  xlab("Critical thermal maximum in °C")+
  ylab("Prosome length in µm")+
  theme(legend.position = "none")

S7_B <- ggplot(all, aes(x = Ctmax, y = length, col = species))+
  geom_point(aes(shape = generation), size = 2.75)+
  scale_color_manual(values = c("#CFD4EB","#D5968F"), name = "species")+
  theme_light(base_size = 14)+
  xlab("Critical thermal maximum in °C")+
  ylab("Prosome length in µm")+
  theme(legend.position = "none")

legendS7 <- get_legend(ggplot(all, aes(x = Ctmax, y = length, col = species_gen))+
  geom_point(aes(shape = generation), size = 2.75)+
  scale_color_manual(values = c("#333333","#CFD4EB", "#96B48E","#D5968F"), name = "species")+
  theme_light(base_size = 14))

dev.new()
combined_S7 <- plot_grid(S7_A, S7_B,  ncol = 1, labels = c("A", "B"))
plot_grid(combined_S7, legendS7, rel_widths = c(3/4, 1/4) )

#### Figure S.8 - CTD data####
p_exp <- read.csv("~/Documents/Scripts/Thesis/CTD2.csv", sep = ";")
p_exp$Date <- paste(p_exp$Year, p_exp$Month, p_exp$Day, sep = "-")

p_exp$Date <- as.Date(p_exp$Date) 

#temperature
S8_A <- ggplot(p_exp, aes(x = Date, y = Pressure..db, color = Temp..degC))+
  geom_point()+
  scale_y_reverse()+
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous"))+
  theme_light()+
  labs(title = "Water temperature at sampling dates",
       y = "Depth [m]", x = "Month", color = "Temp [°C]")

##salinity

S8_B <- ggplot(p_exp, aes(x = Date, y = Pressure..db, color = SALIN..ppt))+
  geom_point()+
  scale_y_reverse()+
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous"))+
  theme_light()+
  labs(title = "Salinity at sampling dates",
       y = "Depth [m]", x = "Month", color = "Salinity [ppt]")

dev.new()
plot_grid(S9_A, S9_B, ncol = 1, labels = c("A", "B"))
