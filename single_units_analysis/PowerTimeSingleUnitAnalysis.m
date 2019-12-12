gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);

keepidx = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 30 ...
     31 32 33 34 35 36 37 38 39 41 42 43 44 45 46 47 48 49 50 51 53 54 55 56 57 58 59 61 62 ...
     64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81];
 contrast = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
 layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};

layer_idx = find(strcmp(layer, 'K'));

log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
 
Ses = struct();
bs_data = struct();
channum = 1: length(log_p_layer);
mean_S = nan(1147,38, length(channum));

xabs = -100:1301;

filtered_dMUA = nan(length(xabs), length(channum));
%dim 2 = channel, dim3 = trials
 Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];

clear i ;
 for i = 1:length(channum)
data = squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(400:1901,:,:));
   bsl = mean(data(1:200,:));
   norm_mean_bs = nan(length(xabs), 1,length(data(1,:)));
   norm_mean_bs(:,1,:) = data(101:end, :) - bsl;
   namelist1(1,1:length(sprintf('chan_%d',i))) = sprintf('chan_%d',i);
   bs_data(i).namelist1 = norm_mean_bs;
clear S namelist;
[S,t,f]        = mtspecgramc(norm_mean_bs(:,1,:) ,movingwin, params); 
 
namelist2(1,1:length(sprintf('S_%d',i))) = sprintf('S_%d',i);
Ses(i).namelist2 = S;
mean_S(:,:,i) = nanmean(S,3);
tvec     = t*1000 + (xabs(1));
%we can also store tvec and f in a struct, but they are all identical
 end
 
 layer_idx = find(strcmp(layer, 'K'));
 layer_mean = mean(mean_S(:,:,layer_idx),3);
 
 figure, 
 imagesc(tvec,f,layer_mean')
 ylim([3.9 60]); 
 set(gca,'ydir','normal')
 title({'DE50_NDE0_su', 'K layer power vs time', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from 71ms before stimulus starts(ms)')
    ylabel('frequency band (Hz)')
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum_plots\',contrast{2},'_power_freq_time_mean_K_layer');
saveas(gcf, strcat(filename, '.png'));  

%% plot the normalized data
layer_idx = find(strcmp(layer, 'P'));
 
 figure, 
 normspec = (nanmean(mean_S(:,:,layer_idx),3) - min(nanmean(mean_S(:,:,layer_idx),3)))./(max(nanmean(mean_S(:,:,layer_idx),3)) - min(nanmean(mean_S(:,:,layer_idx),3)));
 imagesc(tvec,f,normspec')
 ylim([3.9 150]); 
 set(gca,'ydir','normal')
  title({'DE50_NDE0_su', 'P layer power vs time normalized', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from 71ms before stimulus starts(ms)')
    ylabel('frequency band (Hz)')
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum_plots\',contrast{2},'normalized_power_freq_time_mean_P_layer');
saveas(gcf, strcat(filename, '.png')); 
 
 figure, 
 imagesc(tvec,f,nanmedian(S,3)')
 ylim([3.9 60]); 
 set(gca,'ydir','normal')
 
 %%
 %plot the mean data only in the 5Hz range 
 layer_idx = find(strcmp(layer, 'P'));
figure, 
normspec = (nanmean(mean_S(:,:,layer_idx),3) - min(nanmean(mean_S(:,:,layer_idx),3)))./(max(nanmean(mean_S(:,:,layer_idx),3)) - min(nanmean(mean_S(:,:,layer_idx),3)));
plot(tvec,squeeze(normspec(:,1))')
title({'DE50_NDE0_su', 'Mean P layer power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from 71ms before stimulus onset(ms)')
    ylabel('Normalized Power at 4Hz(no units)')
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum_plots\',contrast{2},'normalized_power_freq_time_mean_P_layer_4hz');
saveas(gcf, strcat(filename, '.png')); 
 
 
%plot the median data only in the 5Hz range 
 layer_idx = find(strcmp(layer, 'K'));
figure, 
normspec = (nanmedian(mean_S(:,:,layer_idx),3) - min(nanmedian(mean_S(:,:,layer_idx),3)))./(max(nanmedian(mean_S(:,:,layer_idx),3)) - min(nanmedian(mean_S(:,:,layer_idx),3)));
plot(tvec,squeeze(normspec(:,1))')
title({'DE50_NDE0_su', 'Median K layer power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from 71ms before stimulus onset(ms)')
    ylabel('Normalized Power at 4Hz(no units)')
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum_plots\',contrast{2},'normalized_power_freq_time_median_K_layer_4hz');
saveas(gcf, strcat(filename, '.png')); 

%% plot mean, median, error bars
layer_idx = find(strcmp(layer, 'K'));
figure, 
normspec = (nanmedian(mean_S(:,:,layer_idx),3) - min(nanmedian(mean_S(:,:,layer_idx),3)))./(max(nanmedian(mean_S(:,:,layer_idx),3)) - min(nanmedian(mean_S(:,:,layer_idx),3)));
plot(tvec,squeeze(normspec(:,1))', 'LineWidth',1)
hold on
norm_chan = nan(length(mean_S(:,1,1)), length(layer_idx));
clear i;
for i = 1:length(layer_idx)
min_chan =min(squeeze(mean_S(:,1,layer_idx(i))),[],1);
max_chan = max(squeeze(mean_S(:,1,layer_idx(i))),[],1);
norm_chan(:,i) = (squeeze(mean_S(:,1,layer_idx(i)))-min_chan)./(max_chan - min_chan);
end
%mean_norm_chan = nanmean(norm_chan,2);
normspec = (nanmean(mean_S(:,:,layer_idx),3) - min(nanmean(mean_S(:,:,layer_idx),3)))./(max(nanmean(mean_S(:,:,layer_idx),3)) - min(nanmean(mean_S(:,:,layer_idx),3)));
plot(tvec,squeeze(normspec(:,1))', 'LineWidth',1)

%plot(tvec,squeeze(mean_norm_chan)', 'LineWidth',1)
hold on

ci_low = normspec(:,1) - std(norm_chan,0,2,'omitnan');
 plot(tvec, ci_low, 'LineWidth',1,'Color', 'blue')
 hold on
 ci_high = normspec(:,1) + std(norm_chan,0,2,'omitnan');
 plot(tvec, ci_high, 'LineWidth',1,'Color', 'blue')
 
title({'DE50_NDE0_su', 'Median and mean K layer power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from 71ms before stimulus onset(ms)')
    ylabel('Normalized Power at 4Hz(no units)')
    legend('Median', 'Mean', 'Mean-std', 'Mean+std', 'Location', 'bestoutside')
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum_plots\',contrast{2},'new_normalized_power_freq_time_median_mean_std_K_layer_4hz');
saveas(gcf, strcat(filename, '.png'));

%% Kacie's code 
 load('/Users/kaciedougherty/Downloads/160623_I_p01_uclust4_cinterocdrft_stab_fft_sigmat.mat')
 
 Fs              = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];
?
 xabs         = -100:1301;
 data         = squeeze(channel_data.sdftr_chan(400:1901,:));
 bsl          = mean(data(1:200,:));
 norm_mean_bs = nan(length(xabs),length(data(1,:)));
 norm_mean_bs = data(101:end, :) - bsl;
?
 [S,t,f]        = mtspecgramc(norm_mean_bs ,movingwin, params); 
 
 tvec     = t*1000 + (xabs(1)); 
 figure, 
 imagesc(tvec,f,nanmean(S,3)')
 ylim([3.9 60]); 
 set(gca,'ydir','normal')
 
 figure, 
 normspec = (nanmean(S,3) - min(nanmean(S,3)))./(max(nanmean(S,3)) - min(nanmean(S,3)));
 imagesc(tvec,f,normspec')
 ylim([3.9 150]); 
 set(gca,'ydir','normal')
 
 
 figure, 
 imagesc(tvec,f,nanmedian(S,3)')
 ylim([3.9 60]); 
 set(gca,'ydir','normal')
 
  
figure, 
plot(tvec,squeeze(nanmedian(S(:,1,:),3))')


%% plotting mean median wit Kacie's way
normspec2 = (nanmean(mean_S(:,:,layer_idx),3) - min(nanmean(mean_S(:,:,layer_idx),3)))./(max(nanmean(mean_S(:,:,layer_idx),3)) - min(nanmean(mean_S(:,:,layer_idx),3)));
normspec = (nanmedian(mean_S(:,:,layer_idx),3) - min(nanmedian(mean_S(:,:,layer_idx),3)))./(max(nanmedian(mean_S(:,:,layer_idx),3)) - min(nanmedian(mean_S(:,:,layer_idx),3)));

