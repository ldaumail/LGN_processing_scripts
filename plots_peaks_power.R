library(R.matlab)
library("ggpubr")
library("colorspace")
bumps_data <- readMat("C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/individual_channels_peakadj2/mean_peak1_peak4.mat")

parts_data <- readMat("C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/power_spectrum/part1_part2_norm_power.mat")

#file_names contains both the filenames and the layer
file_names_data <- read.csv("C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/filenames.csv", header = TRUE,sep = ',', na.strings="")

peak1<- t(t(bumps_data$peak1.peak4[,1]))
peak4<- t(t(bumps_data$peak1.peak4[,2]))



#create channel index
channel_idx = t(t(rep(1:74, times =1)))

#data frame

org_data = data.frame(channel_idx,peak1, peak4, file_names_data)


library(tidyr)
#no_na_data <- drop_na(org_data)

row.has.na <- apply(org_data, 1, function(x){any(is.na(x))})
org_data.filt <- org_data[!row.has.na,]


## plot each cell class independently


library(reshape2)
library(ggplot2)

## plot the P cell class
peak1 = org_data.filt$peak1[org_data.filt$layer=="P"]
peak4 = org_data.filt$peak4[org_data.filt$layer=="P"]
p_idx =org_data.filt$channel_idx[org_data.filt$layer=="P"]
p_peakdata = data.frame(p_idx,peak1,peak4)

# convert from wide to long
plot_dat <- melt(p_peakdata, id.var='p_idx')

# plot
svg('C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/plots/p_class_plot_peak1_peak4.svg')
ggplot(plot_dat) + 
  # simple lines
  geom_line(aes(x=paste0('dots ', variable), y=value, group=p_idx, color=p_idx)) +
  scale_colour_gradient(high="#FF3366")+
  # box plots and jitter points, with modified x value
  geom_violin(aes(x=variable, y=value,fill=variable), width=0.6, trim=F)+
  scale_fill_manual(values=c("#FF3366", "#FF6699"))+
  geom_boxplot(aes(x=variable, y=value),width=0.1) +
  geom_jitter(alpha =.5, aes(x=paste0('dots ', variable), y=value, color=p_idx), position = position_jitter(width = .02)) +

  # specify x value order
 scale_x_discrete(limits=c('peak1', 'dots peak1', 'dots peak4',
                           'peak4'))+
scale_y_continuous(limits = c(0, 150), breaks = seq(0,150, by=10)) +
 theme_bw() +
  theme(
    plot.title   = element_text(color = "black", size = 15, face = "bold", hjust=0.5),
    axis.title.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    panel.border = element_blank(),
    axis.text    = element_text(size  = 9),
    axis.line = element_line(colour = "black", size = 1) 
     )+

          labs(
          title = "P cell class comparison: peak1 vs peak4 spiking activity",
          x = "Peak Number",
          y = "Spike rate (spikes/s)" )
dev.off()

## plot the M cells class
peak1 = org_data.filt$peak1[org_data.filt$layer=="M"]
peak4 = org_data.filt$peak4[org_data.filt$layer=="M"]
p_idx =org_data.filt$channel_idx[org_data.filt$layer=="M"]
m_peakdata = data.frame(p_idx,peak1,peak4)

# convert from wide to long
plot_dat <- melt(m_peakdata, id.var='p_idx')

# plot
png('C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/plots/m_class_plot_peak1_peak4.png')
ggplot(plot_dat) + 
  # simple lines
  geom_line(aes(x=paste0('dots ', variable), y=value, group=p_idx, color=p_idx)) +
  scale_colour_gradient(low="#333333")+
  # box plots and jitter points, with modified x value
  geom_violin(aes(x=variable, y=value,fill=variable), width=0.6, trim=F)+
  scale_fill_manual(values=c("#555555", "#999999"))+
  geom_boxplot(aes(x=variable, y=value),width=0.1) +
  geom_jitter(alpha =.5, aes(x=paste0('dots ', variable), y=value, color=p_idx), position = position_jitter(width = .02)) +

  # specify x value order
 scale_x_discrete(limits=c('peak1', 'dots peak1', 'dots peak4',
                           'peak4'))+
scale_y_continuous(limits = c(-20, 150), breaks = seq(-20,150, by=10)) +
 theme_bw() +
  theme(
    plot.title   = element_text(color = "black", size = 15, face = "bold", hjust=0.5),
    axis.title.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    panel.border = element_blank(),
    axis.text    = element_text(size  = 9),
    axis.line = element_line(colour = "black", size = 1) 
     )+

          labs(
          title = "M cell class comparison: peak1 vs peak4 spiking activity",
          x = "Peak Number",
          y = "Spike rate (spikes/s)" )
dev.off()


## plot the K cells class
peak1 = org_data.filt$peak1[org_data.filt$layer=="K"]
peak4 = org_data.filt$peak4[org_data.filt$layer=="K"]
p_idx =org_data.filt$channel_idx[org_data.filt$layer=="K"]
k_peakdata = data.frame(p_idx,peak1,peak4)

# convert from wide to long
plot_dat <- melt(k_peakdata, id.var='p_idx')

# plot
png('C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/plots/k_class_plot_peak1_peak4_5chan.png')
ggplot(plot_dat) + 
  # simple lines
  geom_line(aes(x=paste0('dots ', variable), y=value, group=p_idx, color=p_idx)) +
  scale_colour_gradient(low="#9CCC00")+
  # box plots and jitter points, with modified x value
  geom_violin(aes(x=variable, y=value,fill=variable), width=0.6, trim=F)+
  scale_fill_manual(values=c("#9CCC00", "#CCFF00"))+
  geom_boxplot(aes(x=variable, y=value),width=0.1) +
  geom_jitter(alpha =.5, aes(x=paste0('dots ', variable), y=value, color=p_idx), position = position_jitter(width = .02)) +

  # specify x value order
 scale_x_discrete(limits=c('peak1', 'dots peak1', 'dots peak4',
                           'peak4'))+
scale_y_continuous(limits = c(0, 120), breaks = seq(0,120, by=10)) +
 theme_bw() +
  theme(
    plot.title   = element_text(color = "black", size = 15, face = "bold", hjust=0.5),
    axis.title.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    panel.border = element_blank(),
    axis.text    = element_text(size  = 9),
    axis.line = element_line(colour = "black", size = 1) 
     )+

          labs(
          title = "K cell class comparison: peak1 vs peak4 spiking activity",
          x = "Peak Number",
          y = "Spike rate (spikes/s)" )
dev.off()


#### power plots

window1<- t(t(parts_data$parts[,1]))
window2<- t(t(parts_data$parts[,2]))
#create channel index
channel_idx = t(t(rep(1:74, times =1)))

#data frame

org_data = data.frame(channel_idx,window1, window2, file_names_data)


library(tidyr)
#no_na_data <- drop_na(org_data)

row.has.na <- apply(org_data, 1, function(x){any(is.na(x))})
org_data.filt <- org_data[!row.has.na,]

## plot the P cell class
window1 = org_data.filt$window1[org_data.filt$layer=="P"]
window2 = org_data.filt$window2[org_data.filt$layer=="P"]
p_idx =org_data.filt$channel_idx[org_data.filt$layer=="P"]
p_powerdata = data.frame(p_idx,window1,window2)

# convert from wide to long
plot_dat <- melt(p_powerdata, id.var='p_idx')

# plot
png('C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/power_spectrum/plots/p_class_plot_part1_part2.png')
ggplot(plot_dat) + 
  # simple lines
  geom_line(aes(x=paste0('dots ', variable), y=value, group=p_idx, color=p_idx)) +
  scale_colour_gradient(high="#FF3366")+
  # box plots and jitter points, with modified x value
  geom_violin(aes(x=variable, y=value,fill=variable), width=0.6, trim=F)+
  scale_fill_manual(values=c("#FF3366", "#FF6699"))+
  geom_boxplot(aes(x=variable, y=value),width=0.1) +
  geom_jitter(alpha =.5, aes(x=paste0('dots ', variable), y=value, color=p_idx), position = position_jitter(width = .02)) +

  # specify x value order
 scale_x_discrete(limits=c('window1', 'dots window1', 'dots window2',
                           'window2'))+
scale_y_continuous(limits = c(0, 1), breaks = seq(0,1, by=.1)) +
 theme_bw() +
  theme(
    plot.title   = element_text(color = "black", size = 15, face = "bold", hjust=0.5),
    axis.title.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    panel.border = element_blank(),
    axis.text    = element_text(size  = 9),
    axis.line = element_line(colour = "black", size = 1) 
     )+

          labs(
          title = "P cell class comparison: window1 vs window2 mean power",
          x = "Window Number",
          y = "Power at 4Hz (normalized)" )
dev.off()


## plot the M cell class
window1 = org_data.filt$window1[org_data.filt$layer=="M"]
window2 = org_data.filt$window2[org_data.filt$layer=="M"]
p_idx =org_data.filt$channel_idx[org_data.filt$layer=="M"]
m_powerdata = data.frame(p_idx,window1,window2)

# convert from wide to long
plot_dat <- melt(m_powerdata, id.var='p_idx')

# plot
png('C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/power_spectrum/plots/m_class_plot_part1_part2.png')
ggplot(plot_dat) + 
  # simple lines
  geom_line(aes(x=paste0('dots ', variable), y=value, group=p_idx, color=p_idx)) +
  scale_colour_gradient(low="#333333")+
  # box plots and jitter points, with modified x value
  geom_violin(aes(x=variable, y=value,fill=variable), width=0.6, trim=F)+
  scale_fill_manual(values=c("#555555", "#999999"))+
  geom_boxplot(aes(x=variable, y=value),width=0.1) +
  geom_jitter(alpha =.5, aes(x=paste0('dots ', variable), y=value, color=p_idx), position = position_jitter(width = .02)) +

  # specify x value order
 scale_x_discrete(limits=c('window1', 'dots window1', 'dots window2',
                           'window2'))+
scale_y_continuous(limits = c(0, 1), breaks = seq(0,1, by=.1)) +
 theme_bw() +
  theme(
    plot.title   = element_text(color = "black", size = 15, face = "bold", hjust=0.5),
    axis.title.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    panel.border = element_blank(),
    axis.text    = element_text(size  = 9),
    axis.line = element_line(colour = "black", size = 1) 
     )+

          labs(
          title = "M cell class comparison: window1 vs window2 mean power",
          x = "Window Number",
          y = "Power at 4Hz (normalized)" )
dev.off()


## plot the K cell class
window1 = org_data.filt$window1[org_data.filt$layer=="K"]
window2 = org_data.filt$window2[org_data.filt$layer=="K"]
p_idx =org_data.filt$channel_idx[org_data.filt$layer=="K"]
k_powerdata = data.frame(p_idx,window1,window2)

# convert from wide to long
plot_dat <- melt(k_powerdata, id.var='p_idx')

# plot
png('C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/power_spectrum/plots/k_class_plot_part1_part2_5chan.png')
ggplot(plot_dat) + 
  # simple lines
  geom_line(aes(x=paste0('dots ', variable), y=value, group=p_idx, color=p_idx)) +
  scale_colour_gradient(low="#9CCC00")+
  # box plots and jitter points, with modified x value
  geom_violin(aes(x=variable, y=value,fill=variable), width=0.6, trim=F)+
  scale_fill_manual(values=c("#9CCC00", "#CCFF00"))+
  geom_boxplot(aes(x=variable, y=value),width=0.1) +
  geom_jitter(alpha =.5, aes(x=paste0('dots ', variable), y=value, color=p_idx), position = position_jitter(width = .02)) +

  # specify x value order
 scale_x_discrete(limits=c('window1', 'dots window1', 'dots window2',
                           'window2'))+
scale_y_continuous(limits = c(0, 1), breaks = seq(0,1, by=.1)) +
 theme_bw() +
  theme(
    plot.title   = element_text(color = "black", size = 15, face = "bold", hjust=0.5),
    axis.title.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    panel.border = element_blank(),
    axis.text    = element_text(size  = 9),
    axis.line = element_line(colour = "black", size = 1) 
     )+

          labs(
          title = "K cell class comparison: window1 vs window2 mean power",
          x = "Window Number",
          y = "Power at 4Hz (normalized)" )
dev.off()



