gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);

%this vector is only an indicator of the units conserved after the
%channel selection (significance at 4Hz, quality)
%keepidx = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 30 ...
 %    31 32 33 34 35 36 37 38 39 41 42 43 44 45 46 47 48 49 50 51 53 54 55 56 57 58 59 61 62 ...
  %   64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81];
 contrast = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
 %exclude 160517, (first unit, left empty, it is a K neuron)
 layer = {'','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};

%% compute power spectrum
layer_idx = find(strcmp(layer, 'K'));

log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
 
Ses = struct();
bs_data = struct();
channum = 1: length(log_p_layer);
mean_S_stim = nan(1147+128,38, length(channum));

xabs = -100:1301;

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
mean_S_stim(129:end,:,i) = nanmean(S,3);
%tvec     = t*1000 + (xabs(1));
time_adj = 1:128;
x_stim = cat(2, time_adj , t*1000) ;
%we can also store tvec and f in a struct, but they are all identical
 end
 
 %% plots
 layer_idx = find(strcmp(layer, 'K'));
 layer_mean = mean(mean_S_stim(:,:,layer_idx),3);
 
 figure, 
 imagesc(x_stim,f,layer_mean')
 ylim([3.9 60]); 
 set(gca,'ydir','normal')
 title({'DE50_NDE0_su', 'K layer power vs time', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from stimulus starts(ms)')
    ylabel('frequency band (Hz)')
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum_plots\',contrast{2},'_power_freq_time_mean_K_layer');
saveas(gcf, strcat(filename, '.png'));  

%% plot the normalized data
layer_idx = find(strcmp(layer, 'M'));
 
 figure, 
 normspec = (nanmean(mean_S_stim(:,:,layer_idx),3) - min(nanmean(mean_S_stim(:,:,layer_idx),3)))./(max(nanmean(mean_S_stim(:,:,layer_idx),3)) - min(nanmean(mean_S_stim(:,:,layer_idx),3)));
 imagesc(x_stim,f,normspec')
 ylim([3.9 150]); 
 set(gca,'ydir','normal')
  title({'DE50_NDE0_su', 'P layer power vs time normalized', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from stimulus starts(ms)')
    ylabel('frequency band (Hz)')
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum_plots\',contrast{2},'normalized_power_freq_time_mean_P_layer');
saveas(gcf, strcat(filename, '.png')); 
 
 figure, 
 imagesc(x_stim,f,nanmedian(S,3)')
 ylim([3.9 60]); 
 set(gca,'ydir','normal')
 
 %%
 %plot the mean data only in the 5Hz range 
 layer_idx = find(strcmp(layer, 'P'));
figure, 
normspec = (nanmean(mean_S_stim(:,:,layer_idx),3) - min(nanmean(mean_S_stim(:,:,layer_idx),3)))./(max(nanmean(mean_S_stim(:,:,layer_idx),3)) - min(nanmean(mean_S_stim(:,:,layer_idx),3)));
plot(x_stim,squeeze(normspec(:,1))')
title({'DE50_NDE0_su', 'Mean P layer power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from stimulus onset(ms)')
    ylabel('Normalized Power at 4Hz(no units)')
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum_plots\',contrast{2},'normalized_power_freq_time_mean_P_layer_4hz');
saveas(gcf, strcat(filename, '.png')); 
 
 
%plot the median data only in the 5Hz range 
 layer_idx = find(strcmp(layer, 'K'));
figure, 
normspec = (nanmedian(mean_S_stim(:,:,layer_idx),3) - min(nanmedian(mean_S_stim(:,:,layer_idx),3)))./(max(nanmedian(mean_S_stim(:,:,layer_idx),3)) - min(nanmedian(mean_S_stim(:,:,layer_idx),3)));
plot(x_stim,squeeze(normspec(:,1))')
title({'DE50_NDE0_su', 'Median K layer power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from stimulus onset(ms)')
    ylabel('Normalized Power at 4Hz(no units)')
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum_plots\',contrast{2},'normalized_power_freq_time_median_K_layer_4hz');
saveas(gcf, strcat(filename, '.png')); 

%% plot mean, median, error bars
layer_idx = find(strcmp(layer, 'K'));
figure, 
normspec = (nanmedian(mean_S_stim(:,:,layer_idx),3) - min(nanmedian(mean_S_stim(:,:,layer_idx),3)))./(max(nanmedian(mean_S_stim(:,:,layer_idx),3)) - min(nanmedian(mean_S_stim(:,:,layer_idx),3)));
plot(x_stim,squeeze(normspec(:,1))', 'LineWidth',1)
hold on
norm_chan = nan(length(mean_S_stim(:,1,1)), length(layer_idx));
clear i;
for i = 1:length(layer_idx)
min_chan =min(squeeze(mean_S_stim(:,1,layer_idx(i))),[],1);
max_chan = max(squeeze(mean_S_stim(:,1,layer_idx(i))),[],1);
norm_chan(:,i) = (squeeze(mean_S_stim(:,1,layer_idx(i)))-min_chan)./(max_chan - min_chan);
end
%mean_norm_chan = nanmean(norm_chan,2);
normspec = (nanmean(mean_S_stim(:,:,layer_idx),3) - min(nanmean(mean_S_stim(:,:,layer_idx),3)))./(max(nanmean(mean_S_stim(:,:,layer_idx),3)) - min(nanmean(mean_S_stim(:,:,layer_idx),3)));
plot(x_stim,squeeze(normspec(:,1))', 'LineWidth',1)

%plot(tvec,squeeze(mean_norm_chan)', 'LineWidth',1)
hold on

ci_low = normspec(:,1) - std(norm_chan,0,2,'omitnan');
 plot(x_stim, ci_low, 'LineWidth',1,'Color', 'blue')
 hold on
 ci_high = normspec(:,1) + std(norm_chan,0,2,'omitnan');
 plot(x_stim, ci_high, 'LineWidth',1,'Color', 'blue')
 
title({'DE50_NDE0_su', 'Median and mean K layer power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from 71ms before stimulus onset(ms)')
    ylabel('Normalized Power at 4Hz(no units)')
    legend('Median', 'Mean', 'Mean-std', 'Mean+std', 'Location', 'bestoutside')
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum_plots\',contrast{2},'new_normalized_power_freq_time_median_mean_std_K_layer_4hz');
saveas(gcf, strcat(filename, '.png'));

%% Plot mean with error bars before stimulation and during stimulation separately
layer_idx = find(strcmp(layer, 'M'));
log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
 
Ses = struct();
bs_data = struct();
channum = 1: length(log_p_layer);
mean_S_stim = nan(1147+128,38, length(channum));
mean_S_bl = nan(345+128,38, length(channum));
%mean_S_bl = ;


%compute the power spectrum
%dim 2 = channel, dim3 = trials
 Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];
 
xabs = -100:1301;
clear i ;
 for i = 1:length(channum)
data = squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(1:1901,:,:));
   bsl = mean(data(400:599,:));
   %stim data
   norm_mean_bs = nan(length(xabs), 1,length(data(1,:)));
   norm_mean_bs(:,1,:) = data(500:end, :) - bsl;
   %baseline data
   norm_mean_bsbl = nan(600, 1,length(data(1,:)));
   norm_mean_bsbl(:,1,:) = data(1:600, :) - bsl;

clear S namelist;
[S,t,f]        = mtspecgramc(norm_mean_bs(:,1,:) ,movingwin, params); 
[Sbl,tbl,fbl]      = mtspecgramc(norm_mean_bsbl(:,1,:) ,movingwin, params); 

mean_S_stim(129:end,:,i) = nanmean(S,3);
mean_S_bl(129:end,:,i) = nanmean(Sbl,3);

%we can also store tvec and f in a struct, but they are all identical
 end
 
 
time_adj = -99:28;
time_adj_bl = -1:-1:-128;
x_stim = cat(2, time_adj , t*1000 -100) ;
x_bl = cat(2, time_adj_bl, -tbl*1000);
%length(x_bl)=473
%127 = 600-473
rev_x_bl = -600:1:-128;
%here we compute the individual normalized units necessary for the variance
%for both the baseline data and the stimulus data
norm_chan = nan(length(mean_S_stim(:,1,1)), length(layer_idx));
norm_chan_bl = nan(length(mean_S_bl(:,1,1)), length(layer_idx));
clear i;
for i = 1:length(layer_idx)
min_chan =min(squeeze(mean_S_stim(:,1,layer_idx(i))),[],1);
max_chan = max(squeeze(mean_S_stim(:,1,layer_idx(i))),[],1);
norm_chan(:,i) = (squeeze(mean_S_stim(:,1,layer_idx(i)))-min_chan)./(max_chan - min_chan);
norm_chan_bl(:,i) = (squeeze(mean_S_bl(:,1,layer_idx(i)))-min_chan)./(max_chan - min_chan);
end

normspec = (nanmean(mean_S_stim(:,:,layer_idx),3) - min(nanmean(mean_S_stim(:,:,layer_idx),3)))./(max(nanmean(mean_S_stim(:,:,layer_idx),3)) - min(nanmean(mean_S_stim(:,:,layer_idx),3)));
normspecbl = (nanmean(mean_S_bl(:,:,layer_idx),3) - min(nanmean(mean_S_stim(:,:,layer_idx),3)))./(max(nanmean(mean_S_stim(:,:,layer_idx),3)) - min(nanmean(mean_S_stim(:,:,layer_idx),3)));

figure, 
subplot(1,3,1)
plot(rev_x_bl, squeeze(normspecbl(:,1))','LineWidth',1,'Color', [0,0.470,0.7410])
hold on
ci_low = normspecbl(:,1) - 1.96*std(norm_chan_bl,0,2,'omitnan')./sqrt(length(norm_chan_bl(1,:)));
 plot(rev_x_bl, ci_low,':', 'LineWidth',1,'Color', [.40 .40 .40])
 hold on
 ci_high = normspecbl(:,1) + 1.96*std(norm_chan_bl,0,2,'omitnan')./sqrt(length(norm_chan_bl(1,:)));
 plot(rev_x_bl, ci_high,':', 'LineWidth',1,'Color', [.40 .40 .40])
set(gca, 'linewidth',2)
 set(gca,'box','off') 
  xlim([-500 -100])
  ylim([-0.8 1.2])
  xlabel('Time from stimulus onset(ms)')
   ylabel('Normalized Power at 4Hz(no units)')
subplot(1,3,2:3)
 plot(x_stim,squeeze(normspec(:,1))', 'LineWidth',1, 'Color', [0, 0.4470, 0.7410])
 hold on
 ci_low = normspec(:,1) - 1.96*std(norm_chan,0,2,'omitnan')./sqrt(length(norm_chan(1,:)));
 plot(x_stim, ci_low,':', 'LineWidth',1,'Color', [.40 .40 .40])
 hold on
 ci_high = normspec(:,1) + 1.96*std(norm_chan,0,2,'omitnan')./sqrt(length(norm_chan(1,:)));
 plot(x_stim, ci_high,':', 'LineWidth',1,'Color', [.40 .40 .40])
 plot([0 0], ylim,'k')
 hold on
 plot([1150 1150], ylim,'k')
 xlim([-100 1250])
 ylim([-0.8 1.2])
     set(gca, 'linewidth',2)
      set(gca,'box','off') 
sgtitle({'M class cells mean power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
    legend('Mean', 'Mean-1.96*sem', 'Mean+1.96*sem', 'Location', 'bestoutside')
    
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum\plots\',contrast{2},'new_normalized_power_freq_time_mean_95ci_M_layer_4hz_stim_bl');
saveas(gcf, strcat(filename, '.svg'));
saveas(gcf, strcat(filename, '.png'));
%export_fig(gcf, '-jpg', '-transparent');

%% Plot mean with error bars before stimulation and during stimulation in the same analysis
layer_idx = find(strcmp(layer, 'K'));
log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
 
Ses = struct();
bs_data = struct();
channum = 1: length(log_p_layer);
mean_S_stim = nan(1646+128,38, length(channum));
%compute the power spectrum
%dim 2 = channel, dim3 = trials
 Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];
 

clear i ;
 for i = 1:length(channum)
data = squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(1:1901,:,:));
   bsl = mean(data(400:599,:));
   %stim and bl data
   norm_mean_bs = nan(length(data(:,1)),1,length(data(1,:)));
   norm_mean_bs(:,1,:) = data(1:end,:,:) - bsl;
  

clear S namelist;
[S,t,f]        = mtspecgramc(norm_mean_bs(:,1,:) ,movingwin, params); 

mean_S_stim(129:end,:,i) = nanmean(S,3);

%we can also store tvec and f in a struct, but they are all identical
 end
 
time_adj = -99:28;
x_stim = cat(2, time_adj-500 , t*1000 -600) ;

%here we compute the individual normalized units necessary for the variance
%for both the baseline data and the stimulus data
norm_chan = nan(length(mean_S_stim(:,1,1)), length(layer_idx));
clear i;
for i = 1:length(layer_idx)
min_chan =min(squeeze(mean_S_stim(:,1,layer_idx(i))),[],1);
max_chan = max(squeeze(mean_S_stim(:,1,layer_idx(i))),[],1);
norm_chan(:,i) = (squeeze(mean_S_stim(:,1,layer_idx(i)))-min_chan)./(max_chan - min_chan);
end

normspec = (nanmean(mean_S_stim(:,:,layer_idx),3) - min(nanmean(mean_S_stim(:,:,layer_idx),3)))./(max(nanmean(mean_S_stim(:,:,layer_idx),3)) - min(nanmean(mean_S_stim(:,:,layer_idx),3)));

figure, 
 plot(x_stim,squeeze(normspec(:,1))', 'LineWidth',1, 'Color',[167/255 185/255 54/255])
 %black = [24/255 23/255 23/255] )
 %pink = [229/255, 49/255, 90/255]) 
 hold on
 ci_low = normspec(:,1) - 1.96*std(norm_chan,0,2,'omitnan')./sqrt(length(norm_chan(1,:)));
 plot(x_stim, ci_low,':', 'LineWidth',1,'Color', [.40 .40 .40])
 hold on
 ci_high = normspec(:,1) + 1.96*std(norm_chan,0,2,'omitnan')./sqrt(length(norm_chan(1,:)));
 plot(x_stim, ci_high,':', 'LineWidth',1,'Color', [.40 .40 .40])
 plot([0 0], ylim,'k')
 hold on
 plot([1150 1150], ylim,'k')
 xlim([-600 1250])
 %ylim([-0.8 1.2])
     set(gca, 'linewidth',2)
      set(gca,'box','off') 
      xlabel('Time from stimulus onset(ms)')
   ylabel('Normalized Power at 4Hz(no units)')
title({'K class cells mean power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
    legend('Mean', 'Mean-1.96*sem', 'Mean+1.96*sem', 'Location', 'bestoutside')
    
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum\plots\',contrast{2},'new_normalized_power_freq_time_mean_95ci_K_layer_4hz_gathered_green');
saveas(gcf, strcat(filename, '.svg'));
saveas(gcf, strcat(filename, '.png'));
%export_fig(gcf, '-jpg', '-transparent');

%% Plot mean with error bars before stimulation and during stimulation in the same analysis
%% different way to normalize the data(normalize mean SUA before computing the grand cell class mean
%% plotting with spiking activity significant changes

pvaluesdir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj2\lmer_results\';
 pvalfilename = [pvaluesdir 'lmer_results.csv'];
 pvalues = dlmread(pvalfilename, ',', 1,1);
channeldir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj2\';
 peakvals = load([channeldir 'all_data_peaks']);

layer_idx = find(strcmp(layer, 'P'));
log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
 
Ses = struct();
bs_data = struct();
channum = 1: length(log_p_layer);
mean_S_stim = nan(1646+128,38, length(channum));
%compute the power spectrum
%dim 2 = channel, dim3 = trials
 Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];
 

clear i ;
 for i = 1:length(channum)
data = squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(1:1901,:,:));
   bsl = mean(data(400:599,:));
   %stim and bl data
   norm_mean_bs = nan(length(data(:,1)),1,length(data(1,:)));
   norm_mean_bs(:,1,:) = data(1:end,:,:) - bsl;
  

clear S namelist;
[S,t,f]        = mtspecgramc(norm_mean_bs(:,1,:) ,movingwin, params); 

mean_S_stim(129:end,:,i) = nanmean(S,3);

%we can also store tvec and f in a struct, but they are all identical
 end
 
time_adj = -99:28;
x_stim = cat(2, time_adj-500 , t*1000 -600) ;

%here we compute the individual normalized units necessary for the variance
%for both the baseline data and the stimulus data
norm_chan = nan(length(mean_S_stim(:,1,1)), length(layer_idx));
clear i;
for i = 1:length(layer_idx)
min_chan =min(squeeze(mean_S_stim(:,1,layer_idx(i))),[],1);
max_chan = max(squeeze(mean_S_stim(:,1,layer_idx(i))),[],1);
norm_chan(:,i) = (squeeze(mean_S_stim(:,1,layer_idx(i)))-min_chan)./(max_chan - min_chan);
end

normspec = nanmean(norm_chan,2);

figure, 
 plot(x_stim,normspec', 'LineWidth',1, 'Color',[229/255, 49/255, 90/255])
 xlim([-600 1250])
 ylim([-0.1 1])
 %green[167/255 185/255 54/255])
 %black = [24/255 23/255 23/255] )
 %pink = [229/255, 49/255, 90/255]) 
 hold on
 ci_low = normspec(:,1) - 1.96*std(norm_chan,0,2,'omitnan')./sqrt(length(norm_chan(1,:)));
 plot(x_stim, ci_low,':', 'LineWidth',1,'Color', [.40 .40 .40])
 hold on
 ci_high = normspec(:,1) + 1.96*std(norm_chan,0,2,'omitnan')./sqrt(length(norm_chan(1,:)));
 plot(x_stim, ci_high,':', 'LineWidth',1,'Color', [.40 .40 .40])
 plot([0 0], ylim,'k')
 hold on
 plot([1150 1150], ylim,'k')
 hold on 
 
 cnt = 0;
 all_mean_data = nan(4, length(layer_idx));
  for nunit = 1:length(layer_idx)
 mean_data = nanmean(peakvals.data_peaks(layer_idx(nunit)).namelist,2);
   all_mean_data(:,nunit) = mean_data;
  if all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
      cnt= cnt+1;
      sig_su(:,cnt) = norm_chan(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
  
  end
  mean_sig_su = mean(sig_su,2);
  plot(x_stim, mean_sig_su,  'LineWidth',1)
  
 %ylim([-0.8 1.2])
     set(gca, 'linewidth',2)
      set(gca,'box','off') 
      xlabel('Time from stimulus onset(ms)')
   ylabel('Normalized Power at 4Hz(no units)')
title({'P class cells mean power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
    legend('Mean', 'Mean-1.96*sem', 'Mean+1.96*sem', 'Mean significant decrease', 'Location', 'bestoutside')
    
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum\plots\',contrast{2},'indiv_normalized_power_freq_time_mean_95ci_P_layer_4hz_gathered_pink_sig_suamean');
saveas(gcf, strcat(filename, '.svg'));
saveas(gcf, strcat(filename, '.png'));
%export_fig(gcf, '-jpg', '-transparent');

%% Plot mean with error bars before stimulation and during stimulation in the same analysis
%% different way to normalize the data(normalize mean SUA before computing the grand cell class mean
%% plotting with spiking activity significant changes


channeldir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj2\';
 peakvals = load([channeldir 'all_data_peaks']);
gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
 sig95_idx = load( strcat(gooddatadir,'roc_results95.mat'));
 

layer_idx = find(strcmp(layer, 'P'));
log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
 
Ses = struct();
bs_data = struct();
channum = 1: length(log_p_layer);
mean_S_stim = nan(1646+128,38, length(channum));
%compute the power spectrum
%dim 2 = channel, dim3 = trials
 Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];
 

clear i ;
 for i = 1:length(channum)
data = squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(1:1901,:,:));
   bsl = mean(data(400:599,:));
   %stim and bl data
   norm_mean_bs = nan(length(data(:,1)),1,length(data(1,:)));
   norm_mean_bs(:,1,:) = data(1:end,:,:) - bsl;
  

clear S namelist;
[S,t,f]        = mtspecgramc(norm_mean_bs(:,1,:) ,movingwin, params); 

mean_S_stim(129:end,:,i) = nanmean(S,3);

%we can also store tvec and f in a struct, but they are all identical
 end
 
time_adj = -99:28;
x_stim = cat(2, time_adj-500 , t*1000 -600) ;

%here we compute the individual normalized units necessary for the variance
%for both the baseline data and the stimulus data
norm_chan = nan(length(mean_S_stim(:,1,1)), length(layer_idx));
clear i;
for i = 1:length(layer_idx)
min_chan =min(squeeze(mean_S_stim(:,1,layer_idx(i))),[],1);
max_chan = max(squeeze(mean_S_stim(:,1,layer_idx(i))),[],1);
norm_chan(:,i) = (squeeze(mean_S_stim(:,1,layer_idx(i)))-min_chan)./(max_chan - min_chan);
end

normspec = nanmean(norm_chan,2);

figure, 
 plot(x_stim,normspec', 'LineWidth',1, 'Color',[229/255, 49/255, 90/255])
 xlim([-600 1250])
 ylim([-0.1 1])
 %green[167/255 185/255 54/255])
 %black = [24/255 23/255 23/255] )
 %pink = [229/255, 49/255, 90/255]) 
 hold on
 ci_low = normspec(:,1) - 1.96*std(norm_chan,0,2,'omitnan')./sqrt(length(norm_chan(1,:)));
 plot(x_stim, ci_low,':', 'LineWidth',1,'Color', [.40 .40 .40])
 hold on
 ci_high = normspec(:,1) + 1.96*std(norm_chan,0,2,'omitnan')./sqrt(length(norm_chan(1,:)));
 plot(x_stim, ci_high,':', 'LineWidth',1,'Color', [.40 .40 .40])
 plot([0 0], ylim,'k')
 hold on
 plot([1150 1150], ylim,'k')
 hold on 
 
 cnt = 0;
  for nunit = 1:length(layer_idx)
 
   
part1 = nanmean(norm_chan(601:1147,nunit),1);
part2 = nanmean(norm_chan(1148:1722,nunit), 1);


if part1 > part2 && sig95_idx.all_sigs95(layer_idx(nunit)) ==1
    cnt = cnt +1;
     sig_su(:,cnt) = norm_chan(:,nunit); 
       % plot(x_stim,norm_chan(:, nunit)')
     %hold on
end
     
   
  end
  
  mean_sig_su = mean(sig_su,2);
  plot(x_stim, mean_sig_su,  'LineWidth',1, 'Color',[141/255 140/255 140/255] )
  
 %ylim([-0.8 1.2])
     set(gca, 'linewidth',2)
      set(gca,'box','off') 
      xlabel('Time from stimulus onset(ms)')
   ylabel('Normalized Power at 4Hz(no units)')
title({'P class cells mean power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
    legend('Mean', 'Mean-1.96*sem', 'Mean+1.96*sem', 'Mean significant decrease', 'Location', 'bestoutside')
    
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum\plots\',contrast{2},'indiv_normalized_power_freq_time_mean_95ci_P_layer_4hz_gathered_pink_sig_suamean_pow');
saveas(gcf, strcat(filename, '.svg'));
saveas(gcf, strcat(filename, '.png'));
%export_fig(gcf, '-jpg', '-transparent');


%% Plot mean with error bars before stimulation and during stimulation in the same analysis for one single unit only

Ses = struct();
bs_data = struct();

mean_S_stim = nan(1646+128,38,1);
%compute the power spectrum
%dim 2 = channel, dim3 = trials
 Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];
 


data = squeeze(data_file.good_data(23).channel_data.hypo{1,2}.cont_su(1:1901,:,:));
   bsl = mean(data(400:599,:));
   %stim and bl data
   norm_mean_bs = nan(length(data(:,1)),1,length(data(1,:)));
   norm_mean_bs(:,1,:) = data(1:end,:,:) - bsl;
  

clear S namelist;
[S,t,f]        = mtspecgramc(norm_mean_bs(:,1,:) ,movingwin, params); 

mean_S_stim(129:end,:,1) = nanmean(S,3);

%we can also store tvec and f in a struct, but they are all identical

time_adj = -99:28;
x_stim = cat(2, time_adj-500 , t*1000 -600) ;

%here we compute the individual normalized units necessary for the variance
%for both the baseline data and the stimulus data
norm_chan = nan(length(mean_S_stim(:,1,1)), 1);

min_chan =min(squeeze(mean_S_stim(:,1,1)),[],1);
max_chan = max(squeeze(mean_S_stim(:,1,1)),[],1);
norm_chan(:,1) = (squeeze(mean_S_stim(:,1,1))-min_chan)./(max_chan - min_chan);


normspec = (mean_S_stim(:,:,1) - min(mean_S_stim(:,:,1)))./(max(mean_S_stim(:,:,1)) - min(mean_S_stim(:,:,1)));

figure, 
 plot(x_stim,squeeze(normspec(:,1))', 'LineWidth',1)
 %green = 'Color',[167/255 185/255 54/255])
 %black = [24/255 23/255 23/255] )
 %pink = [229/255, 49/255, 90/255]) 
 hold on
 plot([0 0], ylim,'k')
 hold on
 plot([1150 1150], ylim,'k')
 xlim([-600 1250])
 %ylim([-0.8 1.2])
     set(gca, 'linewidth',2)
      set(gca,'box','off') 
      xlabel('Time from stimulus onset(ms)')
   ylabel('Normalized Power at 4Hz(no units)')
title({'SU Mean power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
   
    
filename = strcat('C:\Users\maier\Documents\adaptation_LGN_abstract_poster\poster\',contrast{2},'new_normalized_power_freq_time_mean_160623_I_p03_uclust3_unfiltered_4hz_gathered');
saveas(gcf, strcat(filename, '.svg'));
saveas(gcf, strcat(filename, '.png'));
%export_fig(gcf, '-jpg', '-transparent');




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
 
 x_stim     = t*1000 + (xabs(1)); 
 figure, 
 imagesc(x_stim,f,nanmean(S,3)')
 ylim([3.9 60]); 
 set(gca,'ydir','normal')
 
 figure, 
 normspec = (nanmean(S,3) - min(nanmean(S,3)))./(max(nanmean(S,3)) - min(nanmean(S,3)));
 imagesc(x_stim,f,normspec')
 ylim([3.9 150]); 
 set(gca,'ydir','normal')
 
 
 figure, 
 imagesc(x_stim,f,nanmedian(S,3)')
 ylim([3.9 60]); 
 set(gca,'ydir','normal')
 
  
figure, 
plot(x_stim,squeeze(nanmedian(S(:,1,:),3))')


%% plotting mean median wit Kacie's way
normspec2 = (nanmean(mean_S_stim(:,:,layer_idx),3) - min(nanmean(mean_S_stim(:,:,layer_idx),3)))./(max(nanmean(mean_S_stim(:,:,layer_idx),3)) - min(nanmean(mean_S_stim(:,:,layer_idx),3)));
normspec = (nanmedian(mean_S_stim(:,:,layer_idx),3) - min(nanmedian(mean_S_stim(:,:,layer_idx),3)))./(max(nanmedian(mean_S_stim(:,:,layer_idx),3)) - min(nanmedian(mean_S_stim(:,:,layer_idx),3)));

