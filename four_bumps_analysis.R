#install package to read matrices
install.packages("R.matlab")

#load R.matlab package and load the data
library(R.matlab)
bumps_data <- readMat("C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/four_bumps_bsl.mat")
peak1<- t(bumps_data$col.peaks[,1])
peak2<- t(bumps_data$col.peaks[,2])
peak3<- t(bumps_data$col.peaks[,3])
peak4<- t(bumps_data$col.peaks[,4])

#file_names contains both the filenames and the layer
file_names_data <- read.csv("C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/filenames.csv", header = TRUE,sep = ',')


col_bumps_data <- t(cbind(peak1,peak2,peak3,peak4))

#Create bump index column
bmp_name <- c("peak1","peak2","peak3","peak4")
bumpnb <- matrix(rep(bmp_name, each =74), ncol=4)
col_bumpnb <- t(cbind(t(bumpnb[,1]),t(bumpnb[,2]),t(bumpnb[,3]),t(bumpnb[,4])))

#create channel index
channel_idx = t(t(rep(1:74, times =4)))

#data frame

org_data = data.frame(channel_idx,col_bumpnb, col_bumps_data, file_names_data)

#anova
install.packages("tidyverse")
library(tidyverse)
anovaPeaks = aov( col_bumps_data ~ col_bumpnb + Error(channel_idx), org_data)
summary(anovaPeaks)

#lmer
library(lmerTest)
library(car)
library(lme4)
linearPeaks = lmer(col_bumps_data ~ col_bumpnb +(1|channel_idx), org_data)
summary(linearPeaks)
  


#posthoc
library(lsmeans)
lsmip(linearPeaks, col_bumpnb) #got an error in this one 

   #other way for posthocs
library(multcompView)
CLD(lsmeans(linearPeaks, c("col_bumpnb")))

   #more detailed way for posthoc
CLD(emmeans(linearPeaks, ~col_bumpnb ))
library(multcomp)
posthoc = emmeans(linearPeaks, "col_bumpnb")
summary(as.glht(update(pairs(posthoc), by = NULL)), test = adjusted("free")) #got an error

    #more simple method #couldn't get it to work, need the package
posthoc(linearPeaks, alpha = 0.05,test = "pairwise", base = 1, alternative = "two.sided",
    statistic = "Wald",
    padj.method = "free")

## lmer accounting for the layers
linearPeaks = lmer(col_bumps_data ~ col_bumpnb*layer +(1|channel_idx), org_data)
summary(linearPeaks)
CLD(lsmeans(linearPeaks, c("col_bumpnb", "layer")))
CLD(emmeans(linearPeaks, ~col_bumpnb|layer ))
library(multcomp)
posthoc = emmeans(linearPeaks, "col_bumpnb", by = "layer")
summary(as.glht(update(pairs(posthoc), by = NULL)))#, test = adjusted("free")) #got an error



# Plot the data as a function of layers
org_data$layer <- factor(org_data$layer, 
                                 levels = c("P", "M", "K"),
                                 labels = c("P", "M", "K"))
org_data$col_bumpnb <- factor(org_data$col_bumpnb, 
                                 levels = c("peak1", "peak2","peak3","peak4"),
                                 labels = c("peak1", "peak2","peak3","peak4"))
library("ggpubr")
library("colorspace")
png('all_layers_plot.png')
ggboxplot(org_data, x = "col_bumpnb", y = "col_bumps_data", color = "layer")
          palette = c("#00AFBB", "#E7B800","#E495A5","#ABB065")
dev.off()

##other plots
   # organizing each leayer's data
p_data = org_data$col_bumps_data[org_data$layer=="P"]
peak_nb = org_data$col_bumpnb[org_data$layer=="P"]
p_idx = org_data$channel_idx[org_data$layer=="P"]
p_lay = data.frame(peak_nb, p_data, p_idx)

m_data = org_data$col_bumps_data[org_data$layer=="M"]
mpeak_nb = org_data$col_bumpnb[org_data$layer=="M"]
m_idx = org_data$channel_idx[org_data$layer=="M"]
m_lay = data.frame(mpeak_nb, m_data, m_idx)

k_data = org_data$col_bumps_data[org_data$layer=="K"]
kpeak_nb = org_data$col_bumpnb[org_data$layer=="K"]
k_idx = org_data$channel_idx[org_data$layer=="K"]
k_lay = data.frame(kpeak_nb, k_data, k_idx)

   #comparisons
my_comparisons = list(c("peak1", "peak2"), c("peak1", "peak3"), c("peak1","peak4"), c("peak2","peak3"), c("peak2", "peak4"), c("peak3","peak4"))


png('p_layer_plot.png')
dodge = position_dodge(0.6)
p_lay %>%
ggplot(aes(x = peak_nb, y = p_data, na.rm = TRUE),add = "dotplot")+
 geom_boxplot(width = .3, position = dodge, alpha = .5) +
  #scale_fill_viridis(discrete = TRUE) +
scale_y_continuous(limits = c(0, 150), breaks = seq(0,150, by=10)) +
stat_compare_means(comparisons = my_comparisons, data = p_lay, label = "p.signif", hide.ns = FALSE, method = "t.test")+
 theme_bw() +
  theme(
    plot.title   = element_text(color = "black", size = 15, face = "bold", hjust=0.5),
    axis.title.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    panel.border = element_blank(),
    axis.text    = element_text(size  = 9), 
    #legend.position = "none",
    #strip.background = element_rect(colour="gray94", fill="gray94")
  )+

          labs(
          title = "P layer peaks",
          x = "Peak Number",
          y = "Spike rate (spike/s)" )
          #palette = c("#00AFBB", "#E7B800","#E495A5","#ABB065")
dev.off()

 

