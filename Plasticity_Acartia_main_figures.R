##########################################################
##                     Main Figures                     ##
##                   Alexandra Hahn                     ##
##########################################################

#main figures for Seasonal plasticity in Acartia hudsonica

setwd("~/Documents/Scripts/Thesis")

#load necessary packages
library(tidyverse)
library(wesanderson)
library(cowplot)
library(ggplotFL)
library(Hmisc)
library(Rmisc)

#import full data set
assays <- read.csv("~/Documents/Scripts/Thesis/assays2.csv", sep=";", header = TRUE)
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

####Figure 1 - SST trend and boxplots just for wild collections####
#wild for hudsonica and tonsa, wild2 for just hudsonica
#CTmax
box <- ggplot(wild2, aes(x = mean2, y = Ctmax, fill = X2.week_mean))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(col = X2.week_mean, alpha = 0.2))+
  theme_light(base_size = 17)+
  scale_fill_gradientn(colors = alpha(wes_palette("Zissou1", type = "continuous"), 0.5))+
  scale_color_gradientn(colors = alpha(wes_palette("Zissou1", type = "continuous"), 0.1))+
  scale_x_discrete(labels =c("6.36" = "Col 1","12.86" = "Col 3","11.44" = "Col 2","16.55" = "Col 4","18.11" = "Col 5"))+
  xlab("")+ ylab(expression("CT"["max"]* " in °C"))+
  theme(legend.position = "none")

#length
box2 <- ggplot(wild2, aes(x = mean2, y = length, fill = X2.week_mean))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(col = X2.week_mean, alpha = 0.2))+
  theme_light(base_size = 17)+
  scale_fill_gradientn(colors = alpha(wes_palette("Zissou1", type = "continuous"), 0.5))+
  scale_color_gradientn(colors = alpha(wes_palette("Zissou1", type = "continuous"), 0.1))+
  scale_x_discrete(labels =c("6.36" = "Col 1","12.86" = "Col 3","11.44" = "Col 2","16.55" = "Col 4","18.11" = "Col 5"))+
  xlab("")+ ylab(expression("Prosome length in µm"))+
  theme(legend.position = "none")

#SST trend
#load temperature data
SST <- read.csv("~/Documents/Scripts/Thesis/temperature.csv", sep = ",")
SST$Date <- as.Date(SST$Date)#transfor to date format


SST_plot <- ggplot(SST,aes(Date, T..IPTS.90)) +
  geom_flquantiles(probs=c(0.025, 0.50, 0.975), fill="red", alpha=0.25) + 
  theme_light(base_size = 17)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("") +
  ylab("SST in °C")

#build boxes for 2-week period prior to each sampling day
rects <- data.frame(
  name = c('col1', 'col2', 'col3', 'col4', 'col5'),
  start = c("2022-03-23", "2022-05-02","2022-05-30", "2022-06-13", "2022-07-05"),
  end = c("2022-04-06", "2022-05-16", "2022-06-13", "2022-06-27", "2022-07-19")
)
rects$start <- as.Date(rects$start)
rects$end <- as.Date(rects$end)

#add boxes and dashed lines indicating sampling days
SST_sampling <- SST_plot + geom_rect(data = rects, inherit.aes=FALSE, mapping=aes(xmin = start, xmax = end,
                                                                    ymin = -Inf, ymax = Inf, fill = name), alpha = 0.35)+
  scale_fill_manual(values = alpha(c("#3B9AB2", "#D5C660", "#E9C624", "#EC7B00", "#F21A00")))+
  theme(legend.position = "none")+
  geom_vline(aes(xintercept = as.numeric(as.Date("2022-04-06"))), col = "black", linetype = 3, size = 1)+
  geom_vline(aes(xintercept = as.numeric(as.Date("2022-05-16"))), col = "black", linetype = 3, size = 1)+
  geom_vline(aes(xintercept = as.numeric(as.Date("2022-06-27"))), col = "black", linetype = 3, size = 1)+
  geom_vline(aes(xintercept = as.numeric(as.Date("2022-07-19"))), col = "black", linetype = 3, size = 1)+
  geom_vline(aes(xintercept = as.numeric(as.Date("2022-06-13"))), col = "darkgrey", linetype = 3, size = 1)

SST_sampling

#plot togehter with boxplots in new window
dev.new()
plot_grid(SST_sampling, box, box2, ncol = 3, labels = c("A", "B", "C"))

####Figure 2 - Boxplots for CTmax####
#boxplots Ctmax
hud2 <- distinct(hudsonica, date_sampling, treatment)%>%
  arrange(date_sampling, treatment)
hud2$yloc <- max(hudsonica$Ctmax)+ .5
hud2$label <- c("a", "b", "c", "b", "a", "c", "c", "b", "c", "c", "b")#compact letter from Ctmax_stats

Ctmax_out <- ggplot(hudsonica, aes(x=date_sampling, y=Ctmax, fill=treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2)+
  theme_light(base_size = 11)+
  scale_x_discrete(labels=c("2022-04-06" = "April 22", "2022-05-16" = " May 22", "2022-06-27" = " June 22" , "2022-07-19" = "July 22"))+
  scale_y_continuous(limits = c(23,31), n.breaks = 5)+
  xlab("")+ ylab("Critical thermal maximum in °C")+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#96B48E","#CFD4EB","#D5968F"), name = "Treatment")

#adding compact letters to display
C_out <- Ctmax_out +
  ylim(NA, max(hudsonica$Ctmax)+ .5)+
  geom_text(data = hud2, aes(y = yloc, label = label),
            position = position_dodge(width = .75))


#boxplots length - hudsonica
Length_out <- ggplot(hudsonica, aes(x=date_sampling, y=length, fill=treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2)+
  theme_light(base_size = 11)+
  scale_x_discrete(labels=c("2022-04-06" = "April 22", "2022-05-16" = " May 22", "2022-06-27" = " June 22" , "2022-07-19" = "July 22"))+
  xlab("Sampling time")+ ylab("Prosome length in µm")+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#96B48E","#CFD4EB","#D5968F"), name = "Treatment")


hud2a <- distinct(hudsonica, date_sampling, treatment)%>%
  arrange(date_sampling, treatment)
hud2a$yloc <- 950
hud2a$label <- c("f", "f", "bcd", "f", "bcde", "ab", "cde", "de", "a", "abc", "ef")#compact letter from Ctmax_stats

L_out <- Length_out +
  ylim(NA, max(hudsonica$length)+ 10)+
  geom_text(data = hud2a, aes(y = yloc, label = label),
            position = position_dodge(width = .75))

#extract legend
legend2 <- get_legend(ggplot(hudsonica, aes(x=date_sampling, y=length, fill=treatment))+
                        geom_boxplot(outlier.shape = NA)+
                        scale_fill_manual(values = c("#96B48E","#CFD4EB","#D5968F"), name = "Treatment")+
                        theme_light())

#combine plots
plotCL <- plot_grid(C_out,  L_out, ncol = 1, labels = c("A", "B"))

#plot to new window
dev.new()
plot_grid(plotCL, legend2, rel_widths = c(5/6, 1/6))


#### Figure 3 - reaction norms####
rnorm <- assays[which(assays$treatment != "wild" & assays$species != "tonsa"),]

#reaction norm for length
rn_length <- ggplot(rnorm, aes(y = length, x = treatment, col = collection))+
  stat_summary_bin(fun=mean, geom="line", aes(group = collection), size = 0.75)  + 
  scale_color_manual(values = alpha(c("#3B9AB2", "#D5C660", "#EC7B00", "#F21A00")), name = "Collection")+
  stat_summary(fun=mean, geom="point", size = 2)+
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.05, size=0.75)+
  theme_light(base_size = 14)+
  geom_point(alpha = 0.2, position = position_dodge(0.1))+
  scale_x_discrete(expand = c(0.1,0))+
  ylab("Prosome length in µm") + xlab("Treatment")+
  theme(legend.position = "none")

#reaction norm for Ctmax
rn_ctmax <- ggplot(rnorm, aes(y = Ctmax, x = treatment, col = collection))+
  stat_summary_bin(fun=mean, geom="line", aes(group = collection), size = 0.75)  + 
  scale_color_manual(values = alpha(c("#3B9AB2", "#D5C660", "#EC7B00", "#F21A00")), name = "Collection")+
  stat_summary(fun=mean, geom="point", size = 2)+
  stat_summary(fun.data = mean_cl_normal, geom="errorbar", width=0.05, size=0.75)+
  theme_light(base_size = 14)+
  geom_point(alpha = 0.2, position = position_dodge(0.1))+
  scale_x_discrete(expand = c(0.1,0))+
  ylab("Critical thermal maximum in °C") + xlab("Treatment")+
  theme(legend.position = "none") 

#extract legend 
legend3 <- get_legend(ggplot(rnorm, aes(y = Ctmax, x = treatment, col = collection))+
                       scale_color_manual(values = alpha(c("#3B9AB2", "#D5C660", "#EC7B00", "#F21A00")), name = "Collection")+
                       geom_point()+
                       theme_light())

#combine two plots in new window 
dev.new()
plot_grid(nrow = 1, rn_ctmax, rn_length, legend3, rel_widths = c(3/7, 3/7, 1/7), labels = c("A", "B", ""))


####Figure 4 - Correlation CTmax and prosome length####  

#data set with wild individuals removed, excluding outliers in cold col-2
data_cw <- hudsonica %>%
  dplyr::filter(treatment != "wild")%>%
  dplyr::filter(collection != "2" | treatment != "cold")
  

wild <- ggplot(wild2, aes(y = Ctmax, x = length, col = sex_confirmed))+
  geom_point(aes(shape = generation), size = 2.75)+
  geom_smooth(method = "lm", aes(group = sex_confirmed, col = sex_confirmed, fill = sex_confirmed))+
  scale_color_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  scale_fill_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  facet_wrap(~treatment)+
  theme_light(base_size = 9)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size = 10))+
  theme(legend.position = "none")+
  xlab("Prosome length in µm")+
  ylab("Critical thermal maximum in °C")


cold_warm <- ggplot(data_cw, aes(y = Ctmax, x = length, col = sex_confirmed))+
  geom_point(shape = 17, size = 2.75)+
  facet_wrap(~treatment)+
  theme_light(base_size = 9)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size = 10))+
  geom_smooth(method = "lm", aes(group = sex_confirmed, col = sex_confirmed, fill = sex_confirmed))+
  scale_color_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  scale_fill_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
  theme(legend.position = "none")+
  ylab("Critical thermal maximum in °C")+
  xlab("Prosome length in µm")


legend4 <- get_legend(ggplot(assays, aes(y = Ctmax, x = length, col = sex_confirmed))+
                        geom_point(aes(shape = generation), size = 2.75)+
                        geom_smooth(method = "lm", aes(group = sex_confirmed, col = sex_confirmed, fill = sex_confirmed))+
                        scale_color_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
                        scale_fill_manual(values = c("#D5968F","#CFD4EB"), name = "sex")+
                        theme_light()+
                        theme(legend.position = "bottom"))

dev.new()
plot_grid(wild, cold_warm, legend4, nrow = 3 , rel_heights  = c(5/10, 4/10, 1/10))

