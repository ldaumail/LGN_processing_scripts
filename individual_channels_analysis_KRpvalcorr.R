library(R.matlab)
library(lmerTest)
library(car)
library(lme4)

require(pbkrtest)


##
## create a for loop to treat all the data at once
filenames <-list.files("C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/individual_channels_peakadj2/", pattern=NULL, all.files=FALSE,  full.names=FALSE)
data_filename_list <- Sys.glob("C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/individual_channels_peakadj2/*.mat")
pvals <- array(0, c(length(data_filename_list),4))
tvals <- 0
df.KR <- 0
p.KR <- 0
for(i in 1:length(data_filename_list)){

filename <- filenames[i]
bumps_data <- readMat(data_filename_list[i])
labels <- c("P1","P2","P3","P4")
peak1<- t(bumps_data$all.pks[1,])
peak2<- t(bumps_data$all.pks[2,])
peak3<- t(bumps_data$all.pks[3,])
peak4<- t(bumps_data$all.pks[4,])
col_bumps_data <- t(cbind(peak1,peak2,peak3,peak4))

#create labels
labels <- c("P1","P2","P3","P4")
bumpnb <- matrix(rep(labels, each =length(bumps_data$all.pks[1,])), ncol=4)
p_label <- t(cbind(t(bumpnb[,1]),t(bumpnb[,2]),t(bumpnb[,3]),t(bumpnb[,4])))

#create trial index
trial_idx = t(t(rep(1:length(bumps_data$all.pks[1,]), times =4)))

#create table
org_data = data.frame(trial_idx, p_label, col_bumps_data)

##fit linear model
linearPeaks = lmer(col_bumps_data ~ p_label +(1|trial_idx), org_data)
linsummary <- summary(linearPeaks)

coefs <- data.frame(coef(summary(linearPeaks)))
df.KR <- get_ddf_Lb(linearPeaks, fixef(linearPeaks))
p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))
pvals[i,] <- p.KR

#for(j in 1:4){ 
#if you wanna consider the distribution normal (not recommended if small sample)
#Vcov <- vcov(linearPeaks, useScale = FALSE) 
#betas <- fixef(linearPeaks) 
#se <- sqrt(diag(Vcov)) 
#zval <- betas / se 
#pvals[i,] <- 2 * pnorm(abs(zval), lower.tail = FALSE)
#pvals[i,] <- 2 * (1 - pnorm(abs(coefs$t.value)))
#}

##save output in a mat file
#filename2 <- paste("C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/individual_channels/lmer_results",filename, ".mat", sep = "")
#writeMat(filename2, table = linsummary)

print(p.KR) 
#print(filename)
#print(linsummary)
}
#print(tvals)
#print(df.KR)
#print(p.KR)
print(pvals)
#save pvalues as a .mat file
filename3 <- paste("C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/individual_channels_peakadj2/lmer_results/lmer_results", ".mat", sep = "")
writeMat(filename3, pvalues = pvals)

#save pvalues as a .csv file
filename4 <- paste("C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/individual_channels_peakadj2/lmer_results/lmer_results", ".csv", sep = "")
write.csv(pvals, file = filename4)

