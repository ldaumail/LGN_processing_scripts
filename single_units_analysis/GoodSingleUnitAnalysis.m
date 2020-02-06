%this script was developped after "significant_single_units_analysis to
%plot the data

gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);

%keepidx = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 30 ...
 %    31 32 33 34 35 36 37 38 39 41 42 43 44 45 46 47 48 49 50 51 53 54 55 56 57 58 59 61 62 ...
 %    64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81];
 f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
 %% save channel mean data time locked on stimulus onset one channel per plot
    %/ normalized between -1 and 1 for each channel (within channel)
 
  channum = 1: length(data_file.good_data);
  f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 %for n = 1:3 
for chan = 1:24:length(channum)
h = figure;
xabs = -50:1301;
idx = [1 3 5 7 9 11 13 15 17 19 21 23 2 4 6 8 10 12 14 16 18 20 22 24];
nyq = 15000;
norm_mean_data = nan(length(xabs), length(1:24));
clear i ;
 for i = 1:24
 
   mean_data = mean(squeeze(data_file.good_data(keepidx(chan+i-1)).channel_data.hypo{1,2}.cont_su(550:1901,:,:)),2);
   
   norm_mean_data(:,i) = (mean_data);
  
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpdSUA      = filtfilt(bwb,bwa, norm_mean_data(:,i));
   
    sp = subplot(length(1:12), 2, idx(i));
    plot(xabs, lpdSUA)
    hold on
    plot([0 0], ylim,'k')
    hold on
    plot([1150 1150], ylim,'k')

    if i == length(12)/2
        ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}Spike Rate (spikes/s)'});
    end
   if i < 12 || (i >= 13)&&(i < 24)
        set(sp, 'XTick', [])
   end
   pvalue = data_file.good_data(chan+i-1).channel_data.hypo{1,2}.cont_stats.pvalue;
      %ylim([-2 2]); 
      ylabelh = text(max(xabs), mean(lpdSUA,1), strcat(num2str(keepidx(chan+i-1)),' | ', num2str(pvalue)),'HorizontalAlignment','left','FontName', 'Arial','FontSize', 10);
      
      
      set(gca, 'linewidth',2)
      set(gca,'box','off')
     % h = subplot(1,1,1); 
     %set(h,'position',get(h,'position').*[1 1 1 1.2]);
  end

    sgtitle({f{2}, 'all good responses, p<0.05'}, 'Interpreter', 'none')
    xlabel('Time from -50ms from stimulus onset (ms)')
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 8.5 11])
    filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\power_channels\good_single_units_data_4bumps_more\plots\',strcat(f{2}, sprintf('_%d_better_raw_data_pvalue05', keepidx(chan))));
   saveas(gcf, strcat(filename, '.png'));

end

 %% Plot the data time locked to the first peak 
 
xabs = 1:1302;
%x = 1:1702;
nyq = 15000;
channum = 1: length(data_file.good_data);
 
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dSUA = nan(length(xabs), length(channum));
all_locsdSUA_filtered = nan(length(channum),1);
up_dist = nan(1, length(channum));

pvalues = nan(1,length(channum));
mean_linreg = nan(2,1);
clear i ;
 for i = 1:length(channum)
     
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(600:1901,:,:)),2);
   
   raw_mean(:,i) = mean_data;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpdSUA      = filtfilt(bwb,bwa, raw_mean(:,i));
  
   %{
   [pkssu, locssu] = findpeaks(lpdMUA(50:1201));
x = xabs - locssu(1);
plot(x, lpdMUA)
hold on
%}
    
 filtered_dSUA(:,i) = lpdSUA;
 %find peaks for every channel that is not accounting for any maxima in the
 %beginning of the trace
 clear len
  for len = 30:550
            if filtered_dSUA(len,i) < filtered_dSUA(len+1,i)
   locsdSUA_filtered = findpeaks(filtered_dSUA(len:1201,i));
       break
            end
      
  end
         

 %store first peak location
 if exist('locsdSUA_filtered', 'var') == 1
 all_locsdSUA_filtered(i) = locsdSUA_filtered.loc(1)+len;
 
 %compute the distance between the first peak and stimulus onset and store  
 %in a matrix
 up_dist(:,i)= length(xabs)- all_locsdSUA_filtered(i);
 end
 end
  h = figure;
 %get the max distance 
 max_low_dist = max(all_locsdSUA_filtered);
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist);
 fp_locked_data = nan(new_dist,length(channum));
 
 for n = 1:length(channum)
 lower_bound =max_low_dist-all_locsdSUA_filtered(n)+1;
 upper_bound =max_low_dist-all_locsdSUA_filtered(n)+length(xabs);
 
 fp_locked_data(lower_bound:upper_bound,n) = filtered_dSUA(:,n);
 x =1:length(fp_locked_data(:,1)); 
 plot(x, fp_locked_data(:,n))
 hold on
 end
  
  mean_filtered = mean(fp_locked_data, 2);
 locsdSUA_mean = findpeaks(mean_filtered(1:length(mean_filtered)));
 x1 = 1:length(mean_filtered);
 plot(x1,mean_filtered,'LineWidth',1, 'Color', 'black')
txt1 = 'mean';
text(x1(400), mean_filtered(400), txt1)
 hold on
 
 median = nanmedian(fp_locked_data, 2);
 median_filtered      = filtfilt(bwb,bwa, median);
 locsdMUA_median = findpeaks(median_filtered(1:length(mean_filtered)));
 x2 = 1:length(mean_filtered); 
 plot(x2,median_filtered,'LineWidth',1, 'Color', 'red')

txt2 = 'median';
 textColor = 'red';
    text(x2(680), median_filtered(680), txt2, 'Color', textColor)
   title({'DE50_NDE0_su', 'all single units of good_data'}, 'Interpreter', 'none')
    xlabel('Time from first peak (ms)')
    ylabel('Spike rate (spike/sec)')
    
    
 %% looking at each bump
   
xabs = 0:1302;
xabs2 = 1:4;
nyq = 15000;
channum = 1: length(data_file.good_data);
 
norm_mean_bs = nan(length(xabs), length(channum));
filtered_dSUA = nan(length(xabs), length(channum));
all_pksdMUA = nan(4, length(channum));
all_locsdMUA = nan(4, length(channum));

clear i ;
 for i = 1:length(channum)
    
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(550:1901,:,:)),2);
   bsl = mean(mean_data(1:50));
   norm_mean_bs(:,i) = mean_data(50:end) - bsl;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpdSUA      = filtfilt(bwb,bwa, norm_mean_bs(:,i));
   filtered_dSUA(:,i) = lpdSUA; 
 for len = 30:550
            if lpdSUA(len) < lpdSUA(len+1)
   locs = findpeaks(lpdSUA(len:1200));
        break
            end
        
  end
         
         if length(locs.loc) >= 4
             %adjust location to the first data point of lpsu (+len),
    lpsulocs = locs.loc(1:4) + len;
   all_pksdSUA(:,i) = lpdSUA(lpsulocs(1:4));
   all_locsdMUA(:,i) = lpsulocs(1:4);
         end 
%locsdSUA = findpeaks(lpdSUA(100:1251));
%if i ~= 10 && i ~= 26 && i ~= 64 %&& i ~= 39


end   
 
  h = figure;
x = [1 2 3 4];% {'peak1', 'peak2','peak3', 'peak4'};
 boxplot(all_pksdSUA', x, 'Labels',  {'peak1', 'peak2','peak3', 'peak4'} )
 ylim([-20 100])
 title({'DE50_NDE0_su', 'peaks distributions'}, 'Interpreter', 'none')
 ylabel('Spike rate (spike/sec)') 
 xlabel('peak number')

filename = [gooddatadir 'four_bumps_bsl'];
col_peaks = all_pksdMUA';
save(strcat(filename,'.mat'), 'col_peaks');

%% save the filenames
filenames = struct();
for i = 1:length(channum)
    filenames(i).filename = data_file.good_data(i).channel_data.filename;
end
 save(strcat(filename,'.mat'), 'filenames');
   
 
 %% Layers analysis
 

 
%% Plot the data time locked to the first peak according to the layer
  
gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);
pvaluesdir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj2\lmer_results\';
 pvalfilename = [pvaluesdir 'lmer_results.csv'];
 pvalues = dlmread(pvalfilename, ',', 1,1);
channeldir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj2\';
 peakvals = load([channeldir 'all_data_peaks']);
 %exclude 160517, (first unit, left empty, it is a K neuron)
 layer = {'','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};

 
layer_idx = find(strcmp(layer, 'M'));

log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);

contrast = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
xabs = -199:1300;
nyq = 15000;
channum = 1: length(data_file.good_data);
 
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dSUA = nan(length(xabs), length(channum));
all_locsdSUA_filtered = nan(length(channum),1);
up_dist = nan(1, length(channum));

clear i ;
 for i = 1:length(channum)
     
     if log_p_layer(i) == 1
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(401:1900,:,:)),2);
   bsl = mean(mean_data(1:200));
   raw_mean_bs(:,i) = mean_data(1:end)- bsl;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpdSUA      = filtfilt(bwb,bwa, raw_mean_bs(:,i));
  
   %{
   [pkssu, locssu] = findpeaks(lpdMUA(50:1201));
x = xabs - locssu(1);
plot(x, lpdMUA)
hold on
%}
    
 filtered_dSUA(:,i) = lpdSUA;
 %find peaks for every channel
  for len = 30:550
            if filtered_dSUA(200+len,i) < filtered_dSUA(200+len+1,i)
   locsdSUA_filtered = findpeaks(filtered_dSUA(200+len:1350,i));
        break
            end
      
  end
         

 %store first peak location
 all_locsdSUA_filtered(i) = locsdSUA_filtered.loc(1)+200+len;
 
  %compute the distance between the first peak and the last datapoint and store  
 %in a matrix
 up_dist(:,i)= length(xabs)- all_locsdSUA_filtered(i);
     end
 end
 
 %get the max distance between the first peak and the stimulus onset
 max_low_dist = max(all_locsdSUA_filtered);
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist);
 clear fp_locked_data norm_fp_locked
 fp_locked_data = nan(new_dist,length(layer_idx));
 norm_fp_locked = nan(new_dist,length(layer_idx));
 for n = 1:length(layer_idx)
    
 lower_bound =max_low_dist-all_locsdSUA_filtered(layer_idx(n))+1;
 upper_bound =max_low_dist-all_locsdSUA_filtered(layer_idx(n))+length(xabs);
 
 fp_locked_data(lower_bound:upper_bound,n) = filtered_dSUA(:,layer_idx(n));
 %x = 1:length(fp_locked_data(:,1));
 %plot(x, fp_locked_data(:,n))
 %hold on
  norm_fp_locked(lower_bound:upper_bound,n) = (fp_locked_data(lower_bound:upper_bound,n)-min(fp_locked_data(lower_bound:upper_bound,n)))/(max(fp_locked_data(lower_bound:upper_bound,n))-min(fp_locked_data(lower_bound:upper_bound,n)));
 end
 
  %get the significant adapting single units 
 clear sig_su mean_sig_su
  cnt = 0;
 all_mean_data = nan(4, length(layer_idx));
 %sig_su = nan(length(fp_locked_data(:,1)),length(layer_idx));
  for nunit = 1:length(layer_idx)
 mean_data = nanmean(peakvals.data_peaks(layer_idx(nunit)).namelist,2);
   all_mean_data(:,nunit) = mean_data;
  if all_mean_data(4,nunit) < all_mean_data(1,nunit) && pvalues(layer_idx(nunit),4) < .05
      cnt= cnt+1;
      sig_su(:,cnt) = norm_fp_locked(:,nunit); 
     % plot(x_stim,norm_chan(:, nunit)')
     %hold on
  end
  
  end
  mean_sig_su = mean(sig_su,2);
 
 h = figure();
  mean_filtered = mean(norm_fp_locked, 2);
 locsdSUA_mean = findpeaks(mean_filtered(1:end));
 x1 = -locsdSUA_mean.loc(1)+1:length(mean_filtered)-locsdSUA_mean.loc(1);
 plot(x1,mean_filtered,'LineWidth',1, 'Color',[24/255 23/255 23/255] )
 ylim([-.1 1.1])
 %green= [167/255 185/255 54/255])
 %black = [24/255 23/255 23/255] )
 %pink = [229/255, 49/255, 90/255])
 hold on
 ci_low = mean_filtered - 1.96*std(norm_fp_locked,0,2)./sqrt(length(norm_fp_locked(1,:)));
 plot(x1, ci_low,':', 'LineWidth',.7,'Color', [24/255 23/255 23/255] )
 hold on
 ci_high = mean_filtered + 1.96*std(norm_fp_locked,0,2)./sqrt(length(norm_fp_locked(1,:)));
 plot(x1, ci_high,':', 'LineWidth',.7,'Color',[24/255 23/255 23/255] )
  hold on
  maxindex = find(~isnan(mean_filtered),1, 'last');
  minindex = find(~isnan(mean_filtered),1,'first');
  sig_su_plot = mean_sig_su(minindex:maxindex);
  plot(x1(minindex:maxindex), sig_su_plot,  'LineWidth',1, 'Color',[141/255 140/255 140/255] )
  hold on
  plot([0 0], ylim,'k')

      set(gca, 'linewidth',2)
      set(gca,'box','off')
   title({'M class cells mean spiking activity'}, 'Interpreter', 'none')
    xlabel('Time from first{\bf peak} (ms)')
    ylabel('Spike rate (spike/sec)')
    legend( 'Mean','Mean-1.96*sem','Mean+1.96*sem','Mean significant decrease','Location', 'bestoutside')
     xlim([-300 1100])
     
    filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\plots\', strcat(contrast{2},'mean_M_class_fpaligned_black_mean_sigdecsua_norm'));
   saveas(gcf, strcat(filename, '.svg'));
   saveas(gcf, strcat(filename, '.png'));
   
  
%% linear regression analysis

layer_idx = find(strcmp(layer, 'K'));

log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
 
xabs = 0:1301;
nyq = 15000;
channum = 1: length(data_file.good_data);
 
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dSUA = nan(length(xabs), length(channum));
all_locsdSUA_filtered = nan(length(channum),1);
up_dist = nan(1, length(channum));

pvalues = nan(1,length(channum));
mean_linreg = nan(2,1);
clear i ;
 for i = 1:length(channum)
     
     if log_p_layer(i) == 1
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(550:1901,:,:)),2);
   bsl = mean(mean_data(1:50));
   raw_mean_bs(:,i) = mean_data(51:end)- bsl;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpdSUA      = filtfilt(bwb,bwa, raw_mean_bs(:,i));
    
 filtered_dSUA(:,i) = lpdSUA;
 %find peaks for every channel
 %find peaks for every channel
  for len = 30:550
            if filtered_dSUA(len,i) < filtered_dSUA(len+1,i)
   locsdSUA_filtered = findpeaks(filtered_dSUA(len:1251,i));
        break
            end
 
  end
         

 %store first peak location
 all_locsdSUA_filtered(i) = locsdSUA_filtered.loc(1)+len;
 
 %compute the distance between the first peak and the last datapoint and store  
 %in a matrix
 up_dist(:,i)= length(xabs)- all_locsdSUA_filtered(i);
     end
 end
 
 %get the max distance between the first peak and the stimulus onset
 max_low_dist = max(all_locsdSUA_filtered);
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist);
 fp_locked_data = nan(new_dist,length(layer_idx));
 for n = 1:length(layer_idx)
    
 lower_bound =max_low_dist-all_locsdSUA_filtered(layer_idx(n))+1;
 upper_bound =max_low_dist-all_locsdSUA_filtered(layer_idx(n))+length(xabs);
 fp_locked_data(lower_bound:upper_bound,n) = filtered_dSUA(:,layer_idx(n));
 
 %{
 x = 1:length(fp_locked_data(:,1));
 plot(x, fp_locked_data(:,n))
 hold on
 %}
 end
 h = figure();
  mean_filtered = mean(fp_locked_data, 2);
 locsdSUA_mean = findpeaks(mean_filtered(1:end));
 x1 = -locsdSUA_mean.loc(1)+1:length(mean_filtered)-locsdSUA_mean.loc(1);
 plot(x1,mean_filtered,'LineWidth',1, 'Color', 'black')
 hold on
  linreg1 = fitlm(locsdSUA_mean.loc(1:3), mean_filtered(locsdSUA_mean.loc(1:3)), 'y ~ x1');
   linreg_coeff1 = table2array(linreg1.Coefficients(1:2,1));
   pvalues1 = table2array(linreg1.Coefficients(2,4));
   ylinreg1 = linreg_coeff1(2,1) .* locsdSUA_mean.loc + linreg_coeff1(1,1);
 
   hold on
  xreg1 =locsdSUA_mean.loc - locsdSUA_mean.loc(1);
  plot(xreg1, ylinreg1, 'Color', 'black');
  txtreg1 = ['y'];%, pvalue =' num2str(pvalues1)];
  text(xreg1(1), ylinreg1(1), txtreg1)
%{  
 ci_low = mean_filtered - std(fp_locked_data,0,2);
 plot(x1, ci_low, 'LineWidth',1,'Color', 'blue')
 hold on
 
 ci_high = mean_filtered + std(fp_locked_data,0,2);
 plot(x1, ci_high, 'LineWidth',1,'Color', 'blue')
 %}

   title({'DE50_NDE0_su', 'K layers single units'}, 'Interpreter', 'none')
    xlabel('Time from first peak (ms)')
    ylabel('Spike rate (spike/sec)')
    legend( 'k mean', ['regression on the mean', newline, ...
        'y = (' num2str(linreg_coeff1(2,1)) ')x + (' num2str(linreg_coeff1(1,1)) ')'], ...
       'Location', 'bestoutside') % 'mean-std','mean+std'
    %% Data analysis with the data aligned on the fourth peak
    
  layer_idx = find(strcmp(layer, 'K'));

log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);


xabs = 0:1301;
%x = 1:1702;
nyq = 15000;
channum = 1: length(data_file.good_data);
 
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dSUA = nan(length(xabs), length(channum));
all_locsdSUA_filtered = nan(length(channum),1);
up_dist = nan(1, length(channum));

pvalues = nan(1,length(channum));
mean_linreg = nan(2,1);
clear i ;
 for i = 1:length(channum)
     
     if log_p_layer(i) == 1
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(550:1901,:,:)),2);
   bsl = mean(mean_data(1:50));
   raw_mean_bs(:,i) = mean_data(51:end)- bsl;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpdSUA      = filtfilt(bwb,bwa, raw_mean_bs(:,i));
  
   %{
   [pkssu, locssu] = findpeaks(lpdMUA(50:1201));
x = xabs - locssu(1);
plot(x, lpdMUA)
hold on
%}
    
 filtered_dSUA(:,i) = lpdSUA;
 %find peaks for every channel
  %find peaks for every channel
  for len = 30:550
            if filtered_dSUA(len,i) < filtered_dSUA(len+1,i)
   locsdSUA_filtered = findpeaks(filtered_dSUA(len:1251,i));
        break
            end
   
  end
         

 %store first peak location
 all_locsdSUA_filtered(i) = locsdSUA_filtered.loc(4)+len;
 
 %compute the distance between the first peak and the last datapoint and store  
 %in a matrix
 up_dist(:,i)= length(xabs)- all_locsdSUA_filtered(i);
    end
end

 
 %get the max distance between the first peak and the stimulus onset
 max_low_dist = max(all_locsdSUA_filtered);
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist);
 fp_locked_data = nan(new_dist,length(layer_idx));
 for n = 1:length(layer_idx)
    if ~isnan(all_locsdSUA_filtered(layer_idx(n)))
 lower_bound =max_low_dist-all_locsdSUA_filtered(layer_idx(n))+1;
 upper_bound =max_low_dist-all_locsdSUA_filtered(layer_idx(n))+length(xabs);
 
 fp_locked_data(lower_bound:upper_bound,n) = filtered_dSUA(:,layer_idx(n));
 x = 1:length(fp_locked_data(:,1));
 plot(x, fp_locked_data(:,n))
 hold on

    end
 end
 clear n
 for n = 1:length(fp_locked_data(1,:))
  if isnan(fp_locked_data(:,n))
     fp_locked_data(:,n) = [];
  end
 end

  mean_filtered = mean(fp_locked_data, 2);
 locsdSUA_mean = findpeaks(mean_filtered);
 x1 = -locsdSUA_mean.loc(4)+1:length(mean_filtered)-locsdSUA_mean.loc(4);

 h = figure();
 
 plot(x1,mean_filtered,'LineWidth',1, 'Color', 'black')
txt1 = 'mean';
text(x1(400), mean_filtered(400), txt1)
 hold on
 
 ci_low = mean_filtered - std(fp_locked_data,0,2);
 plot(x1, ci_low, 'LineWidth',1,'Color', 'blue')
 hold on
 
 ci_high = mean_filtered + std(fp_locked_data,0,2);
 plot(x1, ci_high, 'LineWidth',1,'Color', 'blue')
 %}

   title({'DE50_NDE0_su', 'K layers single units'}, 'Interpreter', 'none')
    xlabel('Time from fourth peak (ms)')
    ylabel('Spike rate (spike/sec)')
    legend( 'mean','Location', 'bestoutside')  
    
    %% Align the data to every peak
    layer_idx = find(strcmp(layer, 'K'));

log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
    
  xabs = -199:1300;
%x = 1:1702;
nyq = 15000;
channum = 1: length(log_p_layer);
 
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dSUA = nan(length(xabs), length(channum),4);
all_locsdSUA_filtered = nan(length(channum),4);
up_dist = nan(length(channum), 4);

 clear peakalign;
for peakalign = 1:4
    clear i ;
 for i = 1:length(log_p_layer)
     
     if log_p_layer(i) == 1
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(401:1900,:,:)),2);
   bsl = mean(mean_data(1:200));
   raw_mean_bs(:,i) = mean_data(1:end)- bsl;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpdSUA      = filtfilt(bwb,bwa, raw_mean_bs(:,i));
  
    
 filtered_dSUA(:,i, peakalign) = lpdSUA;
 
  %find peaks for every channel
  for len = 30:550
            if filtered_dSUA(200+len,i) < filtered_dSUA(200+len+1,i)
   locsdSUA_filtered = findpeaks(filtered_dSUA(200+len:1350,i));
        break
            end
   
  end
         

 %store first peak location
 if length(locsdSUA_filtered.loc) >= 4
 all_locsdSUA_filtered(i, peakalign) = locsdSUA_filtered.loc(peakalign)+200+len;
 
 %compute the distance between the peakalign and the last datapoint and store  
 %in a matrix
 up_dist(i,peakalign)= length(xabs)- all_locsdSUA_filtered(i, peakalign);
 
 end
      end
 end
 end

 %get the max distance between the peakalign and the stimulus onset
 max_low_dist = max(all_locsdSUA_filtered, [], 'all');
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist, [], 'all');
 fp_locked_data = nan(new_dist,length(layer_idx),4);
  %h1 = figure();
 clear peakalign;
for peakalign = 1:4 
   clear n
 for n = 1:length(fp_locked_data(1,:,1))
    if ~isnan(all_locsdSUA_filtered(layer_idx(n),peakalign))
 lower_bound =max_low_dist-all_locsdSUA_filtered(layer_idx(n),peakalign)+1;
 upper_bound =max_low_dist-all_locsdSUA_filtered(layer_idx(n),peakalign)+length(xabs);
 
 fp_locked_data(lower_bound:upper_bound,n, peakalign) = filtered_dSUA(:,layer_idx(n), peakalign);

 %{
 x = 1:length(fp_locked_data(:,1));
 plot(x, fp_locked_data(:,n))
 hold on
%}
    end
 end


clear m;
 for m = length(fp_locked_data(1,:,peakalign)):-1:1
  if isnan(fp_locked_data(:,m, peakalign))
     fp_locked_data(:,m, :) = [];
  end
 end



  mean_filtered = mean(fp_locked_data(:,:,peakalign), 2);
 locsdSUA_mean = findpeaks(mean_filtered);
 x1 = -locsdSUA_mean.loc(peakalign)+1:length(mean_filtered)-locsdSUA_mean.loc(peakalign);
%{
 h = figure();
 
 plot(x1,mean_filtered,'LineWidth',1, 'Color', 'black')
%txt1 = 'mean';
%text(x1(400), mean_filtered(400), txt1)
 hold on
 %{
 ci_low = mean_filtered - std(fp_locked_data,0,2);
 plot(x1, ci_low, 'LineWidth',1,'Color', 'blue')
 hold on
 
 ci_high = mean_filtered + std(fp_locked_data,0,2);
 plot(x1, ci_high, 'LineWidth',1,'Color', 'blue')
 %}

   title({'DE50_NDE0_su', 'P layers single units', sprintf('aligned on peak %d', peakalign)}, 'Interpreter', 'none')
    xlabel(sprintf('Time from %d peak (ms)', peakalign))
    ylabel('Spike rate (spike/sec)')
    legend( 'mean','Location', 'bestoutside')  
   % filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\plots\',strcat(f{2}, sprintf('_%d_aligned_bsl_mean_M_layer', peakalign)));
   %saveas(gcf, strcat(filename, '.png'));
%}
end

 
%% compare peak 1 to peak 4 with the aligned data
f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
   %compare peak 1 with peak 4
   h = figure();
   univx = -150:150;  
  mean_pk1 = mean(fp_locked_data(:,:,1), 2);
     %catch peak 1
  locs1 = findpeaks(mean_pk1);
   x = locs1.loc(1)-150:locs1.loc(1)+150;
   plot(univx, mean_pk1(x),'LineWidth',1, 'Color', '#e6550d')
   hold on
   ci_low = mean_pk1 - 1.96*std(fp_locked_data(:,:,1),0,2)./sqrt(length(fp_locked_data(1,:,1)));
    plot(univx, ci_low(x),':', 'LineWidth',.7,'Color', '#e6550d')
   hold on
    
   ci_high = mean_pk1 + 1.96*std(fp_locked_data(:,:,1),0,2)./sqrt(length(fp_locked_data(1,:,1)));
 plot(univx, ci_high(x),':', 'LineWidth',.7,'Color', '#e6550d')
   hold on
 
   mean_pk4 = mean(fp_locked_data(:,:,4), 2);
   %catch peak 4
   locs2 = findpeaks(mean_pk4);
   x2 =locs2.loc(4)-150:locs2.loc(4)+150;
   plot(univx, mean_pk4(x2),'LineWidth',1, 'Color', '#3182bd')
   
   ci_low4 = mean_pk4 - 1.96*std(fp_locked_data(:,:,4),0,2)./sqrt(length(fp_locked_data(1,:,4)));
    plot(univx, ci_low4(x),':', 'LineWidth',.7,'Color', '#3182bd')
   hold on
    
   ci_high4 = mean_pk4 + 1.96*std(fp_locked_data(:,:,4),0,2)./sqrt(length(fp_locked_data(1,:,4)));
 plot(univx, ci_high4(x),':', 'LineWidth',.7,'Color', '#3182bd')
   hold on
    plot([0 0], ylim,'k')
    
   set(gca, 'linewidth',2)
      set(gca,'box','off') 
   title({'K class cells mean spiking activity', 'peak aligned'}, 'Interpreter', 'none')
   xlabel('Time from peak (ms)')
   ylabel('Spike rate (spikes/sec)')
   legend( 'meanp1','meanp1-1.96*sem1','meanp1+1.96*sem1','meanp4','meanp4-1.96*sem4','meanp4+1.96*sem4','Location', 'bestoutside')  
   filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\plots\',strcat(f{2},'_aligned_bsl_mean_K_layer_peak1_peak4_95ci'));
   saveas(gcf, strcat(filename, '.png'));
   saveas(gcf, strcat(filename, '.svg'));

%% compare peak 1 to every other peak with the aligned data

   %compare peak 1 with peak n
   
  mean_pk1 = mean(fp_locked_data(:,:,1), 2);
     %catch peak 1
   locs1 = findpeaks(mean_pk1);
 for peaknb = 2:4
  mean_pk2 = mean(fp_locked_data(:,:,peaknb), 2);
   h = figure();
   univx = -150:150;
   x = locs1.loc(1)-150:locs1.loc(1)+150;
   plot(univx, mean_pk1(x),'LineWidth',1, 'Color', 'black')
   hold on
   %catch peak 2
   locs2 = findpeaks(mean_pk2);
   x2 =locs2.loc(peaknb)-150:locs2.loc(peaknb)+150;
   plot(univx, mean_pk2(x2),'LineWidth',1, 'Color', 'red')
   title({'DE50_NDE0_su', 'M layers single units', 'peak aligned'}, 'Interpreter', 'none')
    xlabel('Time from peak (ms)')
    ylabel('Spike rate (spike/sec)')
    legend( 'peak1',sprintf('peak%d', peaknb),'Location', 'bestoutside')  
   %filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\plots\',strcat(f{2},sprintf('_aligned_bsl_mean_P_layer_peak%d_peak%d', 1, peaknb )));
  % saveas(gcf, strcat(filename, '.png'));
 end

 %% Align data to every trough
 
    layer_idx = find(strcmp(layer, 'M'));

log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
    
  xabs = -100:1301;
%x = 1:1702;
nyq = 15000;
channum = 1: length(log_p_layer);
 
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dSUA = nan(length(xabs), length(channum),3);
all_locsdSUA_filtered = nan(length(channum),3);
up_dist = nan(length(channum), 3);

pvalues = nan(1,length(channum));
mean_linreg = nan(2,1);
 clear troughalign;
for troughalign = 1:3
    clear i ;
 for i = 1:length(log_p_layer)
     
     if log_p_layer(i) == 1
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(400:1901,:,:)),2);
   bsl = mean(mean_data(1:200));
   raw_mean_bs(:,i) = mean_data(101:end)- bsl;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpdSUA      = filtfilt(bwb,bwa, raw_mean_bs(:,i));
  
    
 filtered_dSUA(:,i, troughalign) = lpdSUA;
 %find peaks for every channel
 locsdSUA_filtered = findpeaks(-filtered_dSUA(100:1251,i, troughalign));
 %store first trough location

   for len = 1:400
   
            if lpdSUA(len) < lpdSUA(len+1)
   locsp = findpeaks(lpdSUA(len:1200));
    
   locst = findpeaks(-lpdSUA(len+locsp.loc(1):1200));
       break
            end
    end
        
         if length(locst.loc) >= 3
             %adjust location to the first data point of lpsu (+len),
    lpsulocs = locst.loc(1:3) +locsp.loc(1) + len;
    all_locsdSUA_filtered(i, troughalign) = lpsulocs(troughalign);
         end 
 %compute the distance between the peakalign and the last datapoint and store  
 %in a matrix
 up_dist(i,troughalign)= length(xabs)- all_locsdSUA_filtered(i, troughalign);
     
     end
 end
end
 %get the max distance between the peakalign and the stimulus onset
 max_low_dist = max(all_locsdSUA_filtered, [], 'all');
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist, [], 'all');
 fp_locked_data = nan(new_dist,length(layer_idx),4);
  %h1 = figure();

 clear troughalign;
for troughalign = 1:3 
   clear n
 for n = 1:length(fp_locked_data(1,:,1))
    if ~isnan(all_locsdSUA_filtered(layer_idx(n),troughalign))
 lower_bound =max_low_dist-all_locsdSUA_filtered(layer_idx(n),troughalign)+1;
 upper_bound =max_low_dist-all_locsdSUA_filtered(layer_idx(n),troughalign)+length(xabs);
 
 fp_locked_data(lower_bound:upper_bound,n, troughalign) = filtered_dSUA(:,layer_idx(n), troughalign);

 %{
 x = 1:length(fp_locked_data(:,1));
 plot(x, fp_locked_data(:,n))
 hold on
%}
    end
 end


clear m;
 for m = length(fp_locked_data(1,:,troughalign)):-1:1
  if isnan(fp_locked_data(:,m, troughalign))
     fp_locked_data(:,m, :) = [];
  end
 end



  mean_filtered = mean(fp_locked_data(:,:,troughalign), 2);
 locsdSUA_mean = findpeaks(-mean_filtered);
 x1 = -locsdSUA_mean.loc(troughalign)+1:length(mean_filtered)-locsdSUA_mean.loc(troughalign);

   h = figure();
 
 plot(x1,mean_filtered,'LineWidth',1, 'Color', 'black')
txt1 = 'mean aligned to trough';
text(x1(400), mean_filtered(400), txt1)
 hold on
 %{
 ci_low = mean_filtered - std(fp_locked_data,0,2);
 plot(x1, ci_low, 'LineWidth',1,'Color', 'blue')
 hold on
 
 ci_high = mean_filtered + std(fp_locked_data,0,2);
 plot(x1, ci_high, 'LineWidth',1,'Color', 'blue')
 %}

   title({'DE50_NDE0_su', 'M layers single units', sprintf('aligned on troughs')}, 'Interpreter', 'none')
    xlabel(sprintf('Time from %d trough (ms)', troughalign))
    ylabel('Spike rate (spike/sec)')
    legend( 'mean','Location', 'bestoutside')  
   filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\plots\',strcat(f{2}, sprintf('_trough_%d_aligned_bsl_mean_M_layer', troughalign)));
  % saveas(gcf, strcat(filename, '.png'));

end

%% plotting all the curves together
layer_idx = find(strcmp(layer, 'K'));

log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
    
  xabs = -100:1301;
%x = 1:1702;
nyq = 15000;
channum = 1: length(log_p_layer);
 
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dSUA = nan(length(xabs), length(channum),3);
all_locsdSUA_filtered = nan(length(channum),3);
up_dist = nan(length(channum), 3);

pvalues = nan(1,length(channum));
mean_linreg = nan(2,1);
 clear troughalign;
for troughalign = 1:3
    clear i ;
 for i = 1:length(log_p_layer)
     
     if log_p_layer(i) == 1
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(400:1901,:,:)),2);
   bsl = mean(mean_data(1:200));
   raw_mean_bs(:,i) = mean_data(101:end)- bsl;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpdSUA      = filtfilt(bwb,bwa, raw_mean_bs(:,i));
  
    
 filtered_dSUA(:,i, troughalign) = lpdSUA;
 %find troughs for every channel

 
    for len = 1:400
   
            if lpdSUA(len) < lpdSUA(len+1)
   locsp = findpeaks(lpdSUA(len:1200));
    
   locst = findpeaks(-lpdSUA(len+locsp.loc(1):1200));
       break
            end
    end
        
         if length(locst.loc) >= 3
       %adjust location to the first data point of lpsu (+len),
    lpsulocs = locst.loc(1:3) +locsp.loc(1) + len;
    all_locsdSUA_filtered(i, troughalign) = lpsulocs(troughalign);
         end
 %compute the distance between the peakalign and the last datapoint and store
 %in a matrix
 up_dist(i,troughalign)= length(xabs)- all_locsdSUA_filtered(i, troughalign);
      end
     end
 end

 %get the max distance between the peakalign and the stimulus onset
 max_low_dist = max(all_locsdSUA_filtered, [], 'all');
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist, [], 'all');
 fp_locked_data = nan(new_dist,length(layer_idx),4);
  %h1 = figure();

 clear troughalign;
for troughalign = 1:3 
   clear n
 for n = 1:length(fp_locked_data(1,:,1))
    if ~isnan(all_locsdSUA_filtered(layer_idx(n),troughalign))
 lower_bound =max_low_dist-all_locsdSUA_filtered(layer_idx(n),troughalign)+1;
 upper_bound =max_low_dist-all_locsdSUA_filtered(layer_idx(n),troughalign)+length(xabs);
 
 fp_locked_data(lower_bound:upper_bound,n, troughalign) = filtered_dSUA(:,layer_idx(n), troughalign);

 %{
 x = 1:length(fp_locked_data(:,1));
 plot(x, fp_locked_data(:,n))
 hold on
%}
    end
 end


clear m;
 for m = length(fp_locked_data(1,:,troughalign)):-1:1
  if isnan(fp_locked_data(:,m, troughalign))
     fp_locked_data(:,m, :) = [];
  end
 end

end


   h = figure();
for troughalign = 1:3
  mean_filtered = mean(fp_locked_data(:,:,troughalign), 2);
 locsdSUA_mean = findpeaks(mean_filtered);
 x1 = -locsdSUA_mean.loc(1)+1:length(mean_filtered)-locsdSUA_mean.loc(1);

 
 plot(x1,mean_filtered,'LineWidth',1) % 'Color', 'black')
%txt1 = sprintf('trough %d', troughalign);
%text(x1(400), mean_filtered(400), txt1)
 hold on
 %{
 ci_low = mean_filtered - std(fp_locked_data,0,2);
 plot(x1, ci_low, 'LineWidth',1,'Color', 'blue')
 hold on
 
 ci_high = mean_filtered + std(fp_locked_data,0,2);
 plot(x1, ci_high, 'LineWidth',1,'Color', 'blue')
 %}

   title({'DE50_NDE0_su', 'M layers single units', sprintf('aligned on troughs')}, 'Interpreter', 'none')
    xlabel('Time from first peak (ms)')
    ylabel('Spike rate (spike/sec)')
    legend( 'mean aligned on trough 1','mean aligned on trough 2', 'mean aligned on trough 3','Location', 'bestoutside')  
   

end
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\plots\',strcat(f{2}, '_trough_aligned_bsl_mean_M_layer'));
saveas(gcf, strcat(filename, '.png'));


%% compare trough 1 to every other trough with the aligned data

   %compare trough 1 with trough n
   
  mean_pk1 = mean(fp_locked_data(:,:,1), 2);
     %catch trough 1
   locs1 = findpeaks(-mean_pk1);
 for troughnb = 2:3
  mean_pk2 = mean(fp_locked_data(:,:,troughnb), 2);
   h = figure();
   univx = -150:150;
   x = locs1.loc(1)-150:locs1.loc(1)+150;
   plot(univx, mean_pk1(x),'LineWidth',1, 'Color', 'black')
   hold on
   %catch trough 2
   locs2 = findpeaks(-mean_pk2);
   x2 =locs2.loc(troughnb)-150:locs2.loc(troughnb)+150;
   plot(univx, mean_pk2(x2),'LineWidth',1, 'Color', 'red')
   title({'DE50_NDE0_su', 'K layers single units', 'trough aligned'}, 'Interpreter', 'none')
    xlabel('Time from trough (ms)')
    ylabel('Spike rate (spike/sec)')
    legend( 'trough1',sprintf('trough%d', troughnb),'Location', 'bestoutside')  
   filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\plots\',strcat(f{2},sprintf('_aligned_bsl_mean_K_layer_trough%d_trough%d', 1, troughnb )));
   %saveas(gcf, strcat(filename, '.png'));
 end