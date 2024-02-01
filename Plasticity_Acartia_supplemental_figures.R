##########################################################
##                 Supplemental Figures                 ##
##                   Alexandra Hahn                     ##
##########################################################

#this script produces the supplemental figures for 
#"Phenotypic Plasticity Drives Seasonal Thermal Tolerance in a Baltic Copepod"

#load necessary packages
library(tidyverse)
library(wesanderson)
library(cowplot)
library(utils)
library(gridExtra)

#setwd("~/your path") set path to download folder

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
  
  f1 <- hudsonica[which(hudsonica$generation == "f1"),]
  
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


#### Figure S.2 - correlation Ctmax and SST ####

#plot for daily temperatures, all wild A. hudsonica included
FS2_daily <- ggplot(wild2, aes(y = Ctmax, x = daily_mean, col = daily_mean))+
  geom_point(size = 2.75)+
  geom_smooth(method = "lm", color = "grey", fill = "lightgrey")+
  scale_fill_gradientn(colors = alpha(wes_palette("Zissou1", type = "continuous")))+
  scale_color_gradientn(colors = alpha(wes_palette("Zissou1", type = "continuous"),0.5))+
  scale_x_continuous(limits = c(5,20), n.breaks = 4)+
  theme_light(base_size = 14)+
  theme(legend.position = "none")+
  ylab("Critical thermal maximum in °C")+
  xlab("Mean SST in °C")

#plot for 2-week mean, all wild A. hudsonica included
FS2_week <- ggplot(wild2, aes(y = Ctmax, x = X2.week_mean, col = X2.week_mean))+
  geom_point(size = 2.75)+
  geom_smooth(method = "lm", color = "grey", fill = "lightgrey")+
  scale_fill_gradientn(colors = alpha(wes_palette("Zissou1", type = "continuous")))+
  scale_color_gradientn(colors = alpha(wes_palette("Zissou1", type = "continuous"),0.5))+
  scale_x_continuous(limits = c(5,20), n.breaks = 4)+
  theme_light(base_size = 14)+
  theme(legend.position = "none")+
  ylab("Critical thermal maximum in °C")+
  xlab("Mean SST in °C")

#extract legend
legend_FS2 <- get_legend(ggplot(wild2, aes(y = Ctmax, x = X2.week_mean, col = collection))+
                           scale_color_manual(values = alpha(c("#3B9AB2", "#D5C660", "#E9C624", "#EC7B00", "#F21A00")))+
                           geom_point(size = 3))

#combine plots and plot in new window
dev.new()
combined_S2 <- plot_grid(FS2_daily,FS2_week, ncol = 1, labels = c("A", "B"))
plot_grid(nrow = 1, combined_S2,  legend_FS2, ncol = 2 , rel_widths = c(3/4, 1/4))
  

cor.test(wild2$X2.week_mean, wild2$Ctmax, method = "pearson")
cor.test(wild2$daily_mean, wild2$Ctmax, method = "pearson")
mS2 <- lm(Ctmax ~ X2.week_mean, data = wild2)
summary(mS2)

#### Figure S.3 - effects of developmental temperature ####

#Ctmax
S3_A <- ggplot(hudsonica, aes(y = Ctmax, x = X2.week_mean, col = sex_confirmed))+
  geom_point(aes(shape = generation), size = 2.75)+
  geom_smooth(method = "lm", aes(group = sex_confirmed, col = sex_confirmed, fill = sex_confirmed))+
  scale_color_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  scale_fill_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  scale_x_continuous()+
  theme_light(base_size = 14)+
  theme(legend.position = "none")+
  xlab("Developmental temperature in °C")+
  ylab(expression("CT"["max"]* " in °C"))

mS3_A <- lm(Ctmax ~ X2.week_mean + sex_confirmed, data = hudsonica)
summary(mS3_A)

#length
S3_B <- ggplot(hudsonica, aes(y = length, x = X2.week_mean, col = sex_confirmed))+
  geom_point(aes(shape = generation), size = 2.75)+
  geom_smooth(method = "lm", aes(group = sex_confirmed, col = sex_confirmed, fill = sex_confirmed))+
  scale_color_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  scale_fill_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  scale_x_continuous()+
  theme_light(base_size = 14)+
  theme(legend.position = "none")+
  xlab("Developmental temperature in °C")+
  ylab("Prosome length in µm")

mS3_B <- lm(length ~ X2.week_mean + sex_confirmed, data = hudsonica)
summary(mS3_B)

#extract legend
legendS3 <- ggplot(hudsonica, aes(y = length, x = X2.week_mean, col = sex_confirmed))+
  geom_point(aes(shape = generation), size = 2.75)+
  geom_smooth(method = "lm", aes(group = sex_confirmed, col = sex_confirmed, fill = sex_confirmed))+
  scale_color_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  scale_fill_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  theme_light(base_size = 14)
  
#combine plots and plot in new window
dev.new()
combined_S3 <- plot_grid(S3_A, S3_B, ncol = 1, labels = c("A", "B"))
plot_grid(nrow = 1, combined_S3, legendS3, ncol = 2 , rel_widths = c(3/4, 1/4))

cor.test(hudsonica$X2.week_mean, hudsonica$Ctmax, method = "pearson")
cor.test(hudsonica$X2.week_mean, hudsonica$length, method = "pearson")
cor.test(wild2$X2.week_mean, wild2$Ctmax, method = "pearson")
cor.test(wild2$X2.week_mean, wild2$length, method = "pearson")


#### Figure S.4 - Ctmax all ####
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

#### Figure S.5 - Boxplots Ctmax (sex separated) ####
#both sexes
Ctmax_both <- ggplot(hudsonica, aes(x=date_sampling, y=Ctmax, fill=treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2)+
  theme_light(base_size = 11)+
  scale_x_discrete(labels=c("2022-04-06" = "Col 1", "2022-05-16" = " Col 2", "2022-06-27" = "Col 4" , "2022-07-19" = "Col 5"))+
  scale_y_continuous(limits = c(23,31), n.breaks = 5)+
  xlab("")+ ylab("Critical thermal maximum in °C")+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#96B48E","#CFD4EB","#D5968F"), name = "Treatment")+
  ggtitle("both sexes")

#only female
Ctmax_f <- ggplot(subset(hudsonica, sex_confirmed == "f"), aes(x=date_sampling, y=Ctmax, fill=treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2)+
  theme_light(base_size = 11)+
  scale_x_discrete(labels=c("2022-04-06" = "Col 1", "2022-05-16" = " Col 2", "2022-06-27" = "Col 4" , "2022-07-19" = "Col 5"))+
  scale_y_continuous(limits = c(23,31), n.breaks = 5)+
  xlab("")+ ylab("Critical thermal maximum in °C")+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#96B48E","#CFD4EB","#D5968F"), name = "Treatment")+
  ggtitle("females")

#only male
Ctmax_m <- ggplot(subset(hudsonica, sex_confirmed == "m"), aes(x=date_sampling, y=Ctmax, fill=treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2)+
  theme_light(base_size = 11)+
  scale_x_discrete(labels=c("2022-04-06" = "Col 1", "2022-05-16" = " Col 2", "2022-06-27" = "Col 4" , "2022-07-19" = "Col 5"))+
  scale_y_continuous(limits = c(23,31), n.breaks = 5)+
  xlab("")+ ylab("Critical thermal maximum in °C")+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#96B48E","#CFD4EB","#D5968F"), name = "Treatment")+
  ggtitle("males")

legendS5 <- get_legend(ggplot(hudsonica, aes(x=date_sampling, y=Ctmax, fill=treatment))+
                       scale_fill_manual(values = c("#96B48E","#CFD4EB","#D5968F"), name = "Treatment")+
                       geom_boxplot(outlier.shape = NA)+
                       geom_point(alpha = 0.2)+
                       theme_light()+
                       theme(legend.position = "bottom"))

dev.new()
combined_S5 <- grid.arrange(Ctmax_both,Ctmax_f,Ctmax_m, ncol = 3, nrow = 1)
grid.arrange(combined_S5,legendS5, nrow = 2, heights = c(5,1))

#### Figure S.6 - Boxplots length (sex separated)####
#both sexes
Length_both <- ggplot(hudsonica, aes(x=date_sampling, y=length, fill=treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2)+
  theme_light(base_size = 11)+
  scale_y_continuous(limits = c(590,930))+
  scale_x_discrete(labels=c("2022-04-06" = "Col 1", "2022-05-16" = " Col 2", "2022-06-27" = "Col 4" , "2022-07-19" = "Col 5"))+
  xlab("")+ ylab("Prosome length in µm")+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#96B48E","#CFD4EB","#D5968F"), name = "Treatment")+
  ggtitle("both sexes")

#only female
Length_f <- ggplot(subset(hudsonica, sex_confirmed == "f"), aes(x=date_sampling, y=length, fill=treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2)+
  theme_light(base_size = 11)+
  scale_y_continuous(limits = c(590,930))+
  scale_x_discrete(labels=c("2022-04-06" = "Col 1", "2022-05-16" = " Col 2", "2022-06-27" = "Col 4" , "2022-07-19" = "Col 5"))+
  xlab("")+ ylab("Prosome length in µm")+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#96B48E","#CFD4EB","#D5968F"), name = "Treatment")+
  ggtitle("females")

#only male
Length_m <- ggplot(subset(hudsonica, sex_confirmed == "m"), aes(x=date_sampling, y=length, fill=treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2)+
  theme_light(base_size = 11)+
  scale_y_continuous(limits = c(590,930))+
  scale_x_discrete(labels=c("2022-04-06" = "Col 1", "2022-05-16" = " Col 2", "2022-06-27" = "Col 4" , "2022-07-19" = "Col 5"))+
  xlab("")+ ylab("Prosome length in µm")+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#96B48E","#CFD4EB","#D5968F"), name = "Treatment")+
  ggtitle("males")

legendS6 <- get_legend(ggplot(hudsonica, aes(x=date_sampling, y=length, fill=treatment))+
                       scale_fill_manual(values = c("#96B48E","#CFD4EB","#D5968F"), name = "Treatment")+
                       geom_boxplot(outlier.shape = NA)+
                       geom_point(alpha = 0.2)+
                       theme_light()+
                       theme(legend.position = "bottom"))

dev.new()
combined_S6 <- grid.arrange(Length_both,Length_f,Length_m, ncol = 3, nrow = 1)
grid.arrange(combined_S6,legendS6, nrow = 2, heights = c(5,1))

#### Figure S.7 - reaction norms ####
#sex separated reaction norm for length
rn_length <- ggplot(f1, aes(y = length, x = treatment, col = collection))+
  stat_summary_bin(fun=mean, geom="line", aes(group = collection), linewidth = 0.75)  + 
  facet_wrap(~sex_confirmed, labeller = labeller(sex_confirmed = c(f = "female", m = "male")))+
  theme_light(base_size = 14)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size = 10))+
  scale_color_manual(values = alpha(c("#3B9AB2", "#D5C660", "#EC7B00", "#F21A00")), name = "Collection")+
  stat_summary(fun=mean, geom="point", size = 2.75)+
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.05, size=0.75)+
  geom_point(alpha = 0.2, position = position_dodge(0.1))+
  scale_x_discrete(expand = c(0.1,0))+
  ylab("Prosome length in µm") + xlab("Treatment")+
  theme(legend.position = "none")

#sex separated reaction norm for Ctmax
rn_ctmax <- ggplot(f1, aes(y = Ctmax, x = treatment, col = collection))+
  stat_summary_bin(fun=mean, geom="line", aes(group = collection), linewidth = 0.75)  + 
  facet_wrap(~sex_confirmed, labeller = labeller(sex_confirmed = c(f = "female", m = "male")))+
  theme_light(base_size = 14)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size = 10))+
  scale_color_manual(values = alpha(c("#3B9AB2", "#D5C660", "#EC7B00", "#F21A00")), name = "Collection")+
  stat_summary(fun=mean, geom="point", size = 2.75)+
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.05, size=0.75)+
  geom_point(alpha = 0.2, position = position_dodge(0.1))+
  scale_x_discrete(expand = c(0.1,0))+
  ylab("Critical thermal maximum in °C") + xlab("Treatment")+
  theme(legend.position = "none")

#extract legend 
legendS7 <- get_legend(ggplot(f1, aes(y = Ctmax, x = treatment, col = collection))+
                        scale_color_manual(values = alpha(c("#3B9AB2", "#D5C660", "#EC7B00", "#F21A00")), name = "Collection")+
                        geom_point(size = 3)+
                        theme_light())

#combine two plots in new window 
dev.new()
plot_grid(nrow = 1, rn_ctmax, rn_length, legendS7, rel_widths = c(3/7, 3/7, 1/7), labels = c("A", "B", ""))

#### Figure S.8 - species tree ####
#https://github.com/rsbrennan/mitotype_pipeline
#pipeline above was used to create the plots

#### Figure S.9 - species identity histogram ####
#subset
known_species <- assays %>%
  filter(species_gen != "not confirmed")%>%
  filter(species_gen != "failed")

#histogram for Ctmax
ggplot() +
  geom_histogram(data = assays, aes(x = Ctmax), binwidth = 0.5, fill = "#656565", alpha = 0.7) +
  geom_histogram(data = known_species, aes(x = Ctmax, fill = species), binwidth = 0.5) +
  scale_fill_manual(values = c("#CFD4EB","#D5968F")) +
  theme_light()
#### Figure S.10 - species identity clustering####
S10_A <- ggplot(all, aes(x = Ctmax, y = length, col = species_gen))+
  geom_point(aes(shape = generation), size = 2.75)+
  scale_color_manual(values = c("#333333","#CFD4EB", "#96B48E","#D5968F"), name = "species")+
  theme_light(base_size = 14)+
  xlab("Critical thermal maximum in °C")+
  ylab("Prosome length in µm")+
  theme(legend.position = "none")

S10_B <- ggplot(all, aes(x = Ctmax, y = length, col = species))+
  geom_point(aes(shape = generation), size = 2.75)+
  scale_color_manual(values = c("#CFD4EB","#D5968F"), name = "species")+
  theme_light(base_size = 14)+
  xlab("Critical thermal maximum in °C")+
  ylab("Prosome length in µm")+
  theme(legend.position = "none")

legendS10 <- get_legend(ggplot(all, aes(x = Ctmax, y = length, col = species_gen))+
  geom_point(aes(shape = generation), size = 2.75)+
  scale_color_manual(values = c("#333333","#CFD4EB", "#96B48E","#D5968F"), name = "species")+
  theme_light(base_size = 14))

dev.new()
combined_S10 <- plot_grid(S10_A, S10_B,  ncol = 1, labels = c("A", "B"))
plot_grid(combined_S10, legendS10, rel_widths = c(3/4, 1/4) )

#### Figure S.11 - correlation Ctmax and length #####

#data set with wild individuals removed, excluding outliers in cold col-2
data_cw <- hudsonica %>%
  dplyr::filter(treatment != "wild")%>%
  dplyr::filter(collection != "2" | treatment != "cold")

wild <- ggplot(wild2, aes(y = Ctmax, x = length, col = sex_confirmed))+
  geom_point(size = 2.75, alpha = 0.5)+
  geom_smooth(method = "lm", aes(group = sex_confirmed, col = sex_confirmed, fill = sex_confirmed))+
  scale_color_manual(values = c("#D5968F","#96B3FF"), name = "sex")+
  scale_fill_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  facet_wrap(~treatment)+
  theme_light(base_size = 10)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size = 10))+
  theme(legend.position = "none")+
  xlab("Prosome length in µm")+
  ylab("Critical thermal maximum in °C")


cold_warm <- ggplot(data_cw, aes(y = Ctmax, x = length, col = sex_confirmed))+
  geom_point(shape = 17, size = 2.75, alpha = 0.4)+
  facet_wrap(~treatment)+
  theme_light(base_size = 10)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size = 10))+
  geom_smooth(method = "lm", aes(group = sex_confirmed, col = sex_confirmed, fill = sex_confirmed))+
  scale_color_manual(values = c("#D5968F","#96B3FF"), name = "sex")+
  scale_fill_manual(values = c("#D5968F","#96B3FF"), name = "sex")+
  theme(legend.position = "none")+
  ylab("Critical thermal maximum in °C")+
  xlab("Prosome length in µm")


legendS11 <- get_legend(ggplot(assays, aes(y = Ctmax, x = length, col = sex_confirmed))+
                        geom_point(aes(shape = generation), size = 2.75)+
                        geom_smooth(method = "lm", aes(group = sex_confirmed, col = sex_confirmed, fill = sex_confirmed))+
                        scale_color_manual(values = c("#D5968F","#96B3FF"), name = "sex")+
                        scale_fill_manual(values = c("#D5968F","#96B3FF"), name = "sex")+
                        theme_light()+
                        theme(legend.position = "bottom"))

dev.new()
combined_S11 <- plot_grid(wild, cold_warm, nrow = 2, labels = c("A", "B", "C"))
plot_grid(combined_S11, legendS11, nrow = 2 , rel_heights  = c(9/10, 1/10))

cor.test(wild2$length, wild2$Ctmax, method = "pearson")
#### Figure S.12 - CTD data####
p_exp <- read.csv("CTD_data.csv", sep = ";")
p_exp$Date <- paste(p_exp$Year, p_exp$Month, p_exp$Day, sep = "-")

p_exp$Date <- as.Date(p_exp$Date) 

#temperature
S12_A <- ggplot(p_exp, aes(x = Date, y = Pressure..db, color = Temp..degC))+
  geom_point()+
  scale_y_reverse()+
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous"))+
  theme_light()+
  labs(title = "Water temperature at sampling dates",
       y = "Depth [m]", x = "Month", color = "Temp [°C]")

##salinity

S12_B <- ggplot(p_exp, aes(x = Date, y = Pressure..db, color = SALIN..ppt))+
  geom_point()+
  scale_y_reverse()+
  scale_color_gradientn(colors = wes_palette("Zissou1", type = "continuous"))+
  theme_light()+
  labs(title = "Salinity at sampling dates",
       y = "Depth [m]", x = "Month", color = "Salinity [ppt]")

dev.new()
plot_grid(S12_A, S12_B, ncol = 1, labels = c("A", "B"))




