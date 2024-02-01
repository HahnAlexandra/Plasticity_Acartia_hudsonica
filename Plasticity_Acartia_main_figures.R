##########################################################
##                    Main Figures                      ##
##                   Alexandra Hahn                     ##
##########################################################

#this script produces the main figures from 
#"Phenotypic Plasticity Drives Seasonal Thermal Tolerance in a Baltic Copepod"

#load necessary packages
library(tidyverse)
library(wesanderson)
library(cowplot)
library(ggplotFL)
library(emmeans)

#setwd("~/your path") set path to download folder

#import full data set
assays <- read.csv("data.csv", sep=";", header = TRUE)
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
  
  f1 <- hudsonica[which(hudsonica$generation == "f1"),]
  
} 

####Figure 1 - Conceptual Figure #####
# no code for this figure

####Figure 2 - SST trend and boxplots just for wild individuals####
#subset wild2 only includes Acartia hudsonica

# boxplot for CTmax
box <- ggplot(wild2, aes(x = mean2, y = Ctmax, fill = X2.week_mean))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(col = X2.week_mean, alpha = 0.2))+
  theme_light(base_size = 17)+
  scale_fill_gradientn(colors = alpha(wes_palette("Zissou1", type = "continuous"), 0.5))+
  scale_color_gradientn(colors = alpha(wes_palette("Zissou1", type = "continuous"), 0.1))+
  scale_x_discrete(labels =c("6.36" = "Col 1","12.86" = "Col 3","11.44" = "Col 2","16.55" = "Col 4","18.11" = "Col 5"))+
  xlab("")+ ylab(expression("CT"["max"]* " in °C"))+
  theme(legend.position = "none")

#box plot for length
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
SST <- read.csv("temperature.csv", sep = ",")
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

#plot togehter with boxplots in new window
dev.new()
plot_grid(SST_sampling, box, box2, ncol = 3, labels = c("A", "B", "C"))

####Figure 3 - Boxplots for CTmax and length####
#boxplots Ctmax
hud2 <- distinct(hudsonica, date_sampling, treatment)%>%
  arrange(date_sampling, treatment)
hud2$yloc <- max(hudsonica$Ctmax)+ .5
hud2$label <- c("a", "b", "c", "b", "a", "c", "c", "b", "c", "c", "b")#compact letter from Ctmax_stats

Ctmax_out <- ggplot(hudsonica, aes(x=date_sampling, y=Ctmax, fill=treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.2)+
  theme_light(base_size = 11)+
  scale_x_discrete(labels=c("2022-04-06" = "Col 1", "2022-05-16" = " Col 2", "2022-06-27" = "Col 4" , "2022-07-19" = "Col 5"))+
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
  scale_x_discrete(labels=c("2022-04-06" = "Col 1", "2022-05-16" = "Col 2", "2022-06-27" = "Col 4" , "2022-07-19" = "Col 5"))+
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
legend3 <- get_legend(ggplot(hudsonica, aes(x=date_sampling, y=length, fill=treatment))+
                        geom_boxplot(outlier.shape = NA)+
                        scale_fill_manual(values = c("#96B48E","#CFD4EB","#D5968F"), name = "Treatment")+
                        theme_light())

#combine plots
plotCL <- plot_grid(C_out,  L_out, ncol = 1, labels = c("A", "B"))

#plot to new window
dev.new()
plot_grid(plotCL, legend3, rel_widths = c(5/6, 1/6))


#### Figure 4 - reaction norms####

#build table from model means, see Plasticity_Acartia_stats
m_lf1 <- lm(length~treatment * collection + sex_confirmed, data = f1)
means_length <- emmeans(m_lf1, ~ treatment * collection)
summary_means_length <- summary(means_length)
summary_means_length <-  summary_means_length[!(is.na(summary_means_length$emmean)),]

#build reaction norm for length with raw data from subset f1 and model means
rn_length <- ggplot(f1, aes(y = length, x = treatment, col = collection))+
  geom_point(data = summary_means_length, 
             aes(y = emmean, x = treatment, col = collection),
             position = position_dodge(0.1), size = 3)+
  geom_errorbar(data = summary_means_length,
                aes(y = emmean, ymin = (emmean - SE), ymax = (emmean + SE), x = treatment, col = collection),
                position = position_dodge(0.1), width = 0.2) +
  geom_line(data = summary_means_length,
            aes(y = emmean, x = treatment, col = collection, group = collection),
            position = position_dodge(0.1),
            linetype = "solid") +
  scale_color_manual(values = alpha(c("#3B9AB2", "#D5C660", "#EC7B00", "#F21A00")), name = "Collection")+
  theme_light(base_size = 14)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size = 10))+
  geom_point(alpha = 0.2, position = position_dodge(0.1))+
  scale_x_discrete(expand = c(0.1,0))+
  ylab("Critical thermal maximum in °C") + xlab("Treatment")+
  theme(legend.position = "none") 

#reaction norm for Ctmax
#build table from model means, see Plasticity_Acartia_stats
m_f1 <- lm(Ctmax~treatment * collection + sex_confirmed + length, data = f1)
means_Ctmax <- emmeans(m_f1, ~ treatment * collection)
summary_means_Ctmax <- summary(means_Ctmax)
summary_means_Ctmax <-  summary_means_Ctmax[!(is.na(summary_means_Ctmax$emmean)),]

#build reaction norm for Ctmax with raw data from subset f1 and model means
rn_ctmax <- ggplot(f1, aes(y = Ctmax, x = treatment, col = collection))+
  geom_point(data = summary_means_Ctmax, 
             aes(y = emmean, x = treatment, col = collection),
             position = position_dodge(0.1), size = 3)+
  geom_errorbar(data = summary_means_Ctmax,
                aes(y = emmean, ymin = (emmean - SE), ymax = (emmean + SE), x = treatment, col = collection),
                position = position_dodge(0.1), width = 0.2) +
  geom_line(data = summary_means_Ctmax,
            aes(y = emmean, x = treatment, col = collection, group = collection),
            position = position_dodge(0.1),
            linetype = "solid") +
  scale_color_manual(values = alpha(c("#3B9AB2", "#D5C660", "#EC7B00", "#F21A00")), name = "Collection")+
  theme_light(base_size = 14)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size = 10))+
  geom_point(alpha = 0.2, position = position_dodge(0.1))+
  scale_x_discrete(expand = c(0.1,0))+
  ylab("Critical thermal maximum in °C") + xlab("Treatment")+
  theme(legend.position = "none") 

#extract legend 
legend4 <- get_legend(ggplot(f1, aes(y = Ctmax, x = treatment, col = collection))+
                       scale_color_manual(values = alpha(c("#3B9AB2", "#D5C660", "#EC7B00", "#F21A00")), name = "Collection")+
                       geom_point(size = 3)+
                       theme_light())

#combine two plots in new window 
dev.new()
plot_grid(nrow = 1, rn_ctmax, rn_length, legend4, rel_widths = c(3/7, 3/7, 1/7), labels = c("A", "B", ""))

