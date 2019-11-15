gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);

keepidx = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 30 ...
     31 32 33 34 35 36 37 38 39 41 42 43 44 45 46 47 48 49 50 51 53 54 55 56 57 58 59 61 62 ...
     64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81];
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
   lpdMUA      = filtfilt(bwb,bwa, norm_mean_data(:,i));
   
    sp = subplot(length(1:12), 2, idx(i));
    plot(xabs, lpdMUA)
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
      ylabelh = text(max(xabs), mean(lpdMUA,1), strcat(num2str(keepidx(chan+i-1)),' | ', num2str(pvalue)),'HorizontalAlignment','left','FontName', 'Arial','FontSize', 10);
      
      
      set(gca, 'linewidth',2)
      set(gca,'box','off')
     % h = subplot(1,1,1); 
     %set(h,'position',get(h,'position').*[1 1 1 1.2]);
  end

    sgtitle({f{2}, 'all good responses, p<0.05'}, 'Interpreter', 'none')
    xlabel('Time from -50ms from stimulus onset (ms)')
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 8.5 11])
    filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\plots\',strcat(f{2}, sprintf('_%d_better_raw_data_pvalue05', keepidx(chan))));
   saveas(gcf, strcat(filename, '.png'));

end

 %% Plot the data time locked to the first peak 
   h = figure;
xabs = -50:1301;
x = 1:1702;
nyq = 15000;
channum = 1: length(data_file.good_data);
 
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dMUA = nan(length(xabs), length(channum));
all_locsdMUA_filtered = nan(length(channum),1);
up_dist = nan(1, length(channum));

pvalues = nan(1,length(channum));
mean_linreg = nan(2,1);
clear i ;
 for i = 1:length(channum)
     
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(550:1901,:,:)),2);
   
   raw_mean_bs(:,i) = mean_data;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpdMUA      = filtfilt(bwb,bwa, raw_mean_bs(:,i));
  
   %{
   [pkssu, locssu] = findpeaks(lpdMUA(50:1201));
x = xabs - locssu(1);
plot(x, lpdMUA)
hold on
%}
    
 filtered_dMUA(:,i) = lpdMUA;
 %find peaks for every channel
 [pksdMUA_filtered, locsdMUA_filtered] = findpeaks(filtered_dMUA(50:1201,i));
 %store first peak location
 all_locsdMUA_filtered(i) = locsdMUA_filtered(1);
 %compute the distance between the first peak and stimulus onset and store  
 %in a matrix
 up_dist(:,i)= length(xabs)- all_locsdMUA_filtered(i);
 end
 
 %get the max distance 
 max_low_dist = max(all_locsdMUA_filtered);
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist);
 fp_locked_data = nan(new_dist,length(channum));
 for n = 1:length(channum)
 lower_bound =max_low_dist-all_locsdMUA_filtered(n)+1;
 upper_bound =max_low_dist-all_locsdMUA_filtered(n)+length(xabs);
 
 fp_locked_data(lower_bound:upper_bound,n) = filtered_dMUA(:,n);
  
 plot(x, fp_locked_data(:,n))
 hold on
 end
 
  mean_filtered = mean(fp_locked_data, 2);
 [pksdMUA_mean, locsdMUA_mean] = findpeaks(mean_filtered(1:1586));
 x1 = 1:1702;
 plot(x1,mean_filtered,'LineWidth',1, 'Color', 'black')
txt1 = 'mean';
text(x1(400), mean_filtered(400), txt1)
 hold on
 
 median = nanmedian(fp_locked_data, 2);
 median_filtered      = filtfilt(bwb,bwa, median);
 [pksdMUA_median, locsdMUA_median] = findpeaks(median_filtered(1:1586));
 x2 = 1:1702; 
 plot(x2,median_filtered,'LineWidth',1, 'Color', 'red')

txt2 = 'median';
 textColor = 'red';
    text(x2(680), median_filtered(680), txt2, 'Color', textColor)
   title({'DE50_NDE0_su', 'all single units of good_data'}, 'Interpreter', 'none')
    xlabel('Time from first peak (ms)')
    ylabel('Spike rate (spike/sec)')
    
    
 %% looking at each bump
   
xabs = -100:1301;
xabs2 = 1:4;
nyq = 15000;
channum = 1: length(data_file.good_data);
 
norm_mean_bs = nan(length(xabs), length(channum));
filtered_dMUA = nan(length(xabs), length(channum));
all_pksdMUA = nan(4, length(channum));
all_locsdMUA = nan(4, length(channum));

clear i ;
 for i = 1:length(channum)
    
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(400:1901,:,:)),2);
   bsl = mean(mean_data(1:200));
   norm_mean_bs(:,i) = mean_data(101:end) - bsl;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpdMUA      = filtfilt(bwb,bwa, norm_mean_bs(:,i));
   filtered_dMUA(:,i) = lpdMUA; 
 
[pksdMUA, locsdMUA] = findpeaks(lpdMUA(100:1251));
if i ~= 10 && i ~= 26 && i ~= 39 && i ~= 64
all_pksdMUA(:,i) = pksdMUA(1:4);
all_locsdMUA(:,i) = locsdMUA(1:4);
end   
 end
  h = figure;
x = [1 2 3 4];% {'peak1', 'peak2','peak3', 'peak4'};
 boxplot(all_pksdMUA', x, 'Labels',  {'peak1', 'peak2','peak3', 'peak4'} )
 ylim([0 100])
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
 layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};

layer_idx = find(strcmp(layer, 'M'));

log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);

 
%% Plot the data time locked to the first peak according to the layer
  
xabs = -100:1301;
%x = 1:1702;
nyq = 15000;
channum = 1: length(data_file.good_data);
 
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dMUA = nan(length(xabs), length(channum));
all_locsdMUA_filtered = nan(length(channum),1);
up_dist = nan(1, length(channum));

pvalues = nan(1,length(channum));
mean_linreg = nan(2,1);
clear i ;
 for i = 1:length(channum)
     
     if log_p_layer(i) == 1
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(400:1901,:,:)),2);
   bsl = mean(mean_data(1:200));
   raw_mean_bs(:,i) = mean_data(101:end)- bsl;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpdMUA      = filtfilt(bwb,bwa, raw_mean_bs(:,i));
  
   %{
   [pkssu, locssu] = findpeaks(lpdMUA(50:1201));
x = xabs - locssu(1);
plot(x, lpdMUA)
hold on
%}
    
 filtered_dMUA(:,i) = lpdMUA;
 %find peaks for every channel
 [pksdMUA_filtered, locsdMUA_filtered] = findpeaks(filtered_dMUA(100:1251,i));
 %store first peak location
 all_locsdMUA_filtered(i) = locsdMUA_filtered(1);
 %compute the distance between the first peak and the last datapoint and store  
 %in a matrix
 up_dist(:,i)= length(xabs)- all_locsdMUA_filtered(i);
     end
 end
 
 %get the max distance between the first peak and the stimulus onset
 max_low_dist = max(all_locsdMUA_filtered);
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist);
 fp_locked_data = nan(new_dist,length(layer_idx));
 for n = 1:length(layer_idx)
    
 lower_bound =max_low_dist-all_locsdMUA_filtered(layer_idx(n))+1;
 upper_bound =max_low_dist-all_locsdMUA_filtered(layer_idx(n))+length(xabs);
 
 fp_locked_data(lower_bound:upper_bound,n) = filtered_dMUA(:,layer_idx(n));
 x = 1:length(fp_locked_data(:,1));
 plot(x, fp_locked_data(:,n))
 hold on
 
 end
 h = figure();
  mean_filtered = mean(fp_locked_data, 2);
 [pksdMUA_mean, locsdMUA_mean] = findpeaks(mean_filtered(1:end));
 x1 = -locsdMUA_mean(1)+1:length(mean_filtered)-locsdMUA_mean(1);
 plot(x1,mean_filtered,'LineWidth',1, 'Color', 'black')
txt1 = 'mean';
text(x1(400), mean_filtered(400), txt1)
 hold on
 
 ci_low = mean_filtered - std(fp_locked_data,0,2);
 plot(x1, ci_low, 'LineWidth',1,'Color', 'blue')
 hold on
 
 ci_high = mean_filtered + std(fp_locked_data,0,2);
 plot(x1, ci_high, 'LineWidth',1,'Color', 'blue')
 
 %{
 median = nanmedian(fp_locked_data, 2);
 median_filtered      = filtfilt(bwb,bwa, median);
 [pksdMUA_median, locsdMUA_median] = findpeaks(median_filtered(1:end));
 x2 = 1:length(median_filtered); 
 plot(x2,median_filtered,'LineWidth',1, 'Color', 'red')

txt2 = 'median';
 textColor = 'red';
    text(x2(680), median_filtered(680), txt2, 'Color', textColor)
 %}
   title({'DE50_NDE0_su', 'M layers single units'}, 'Interpreter', 'none')
    xlabel('Time from first peak (ms)')
    ylabel('Spike rate (spike/sec)')
    legend( 'M','Location', 'bestoutside')
    
%% linear regreassion analysis

layer_idx = find(strcmp(layer, 'K'));

log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
 
xabs = -100:1301;
nyq = 15000;
channum = 1: length(data_file.good_data);
 
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dMUA = nan(length(xabs), length(channum));
all_locsdMUA_filtered = nan(length(channum),1);
up_dist = nan(1, length(channum));

pvalues = nan(1,length(channum));
mean_linreg = nan(2,1);
clear i ;
 for i = 1:length(channum)
     
     if log_p_layer(i) == 1
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(400:1901,:,:)),2);
   bsl = mean(mean_data(1:200));
   raw_mean_bs(:,i) = mean_data(101:end)- bsl;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpdMUA      = filtfilt(bwb,bwa, raw_mean_bs(:,i));
    
 filtered_dMUA(:,i) = lpdMUA;
 %find peaks for every channel
 [pksdMUA_filtered, locsdMUA_filtered] = findpeaks(filtered_dMUA(100:1251,i));
 %store first peak location
 all_locsdMUA_filtered(i) = locsdMUA_filtered(1);
 %compute the distance between the first peak and the last datapoint and store  
 %in a matrix
 up_dist(:,i)= length(xabs)- all_locsdMUA_filtered(i);
     end
 end
 
 %get the max distance between the first peak and the stimulus onset
 max_low_dist = max(all_locsdMUA_filtered);
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist);
 fp_locked_data = nan(new_dist,length(layer_idx));
 for n = 1:length(layer_idx)
    
 lower_bound =max_low_dist-all_locsdMUA_filtered(layer_idx(n))+1;
 upper_bound =max_low_dist-all_locsdMUA_filtered(layer_idx(n))+length(xabs);
 fp_locked_data(lower_bound:upper_bound,n) = filtered_dMUA(:,layer_idx(n));
 
 %{
 x = 1:length(fp_locked_data(:,1));
 plot(x, fp_locked_data(:,n))
 hold on
 %}
 end
 h = figure();
  mean_filtered = mean(fp_locked_data, 2);
 [pksdMUA_mean, locsdMUA_mean] = findpeaks(mean_filtered(1:end));
 x1 = -locsdMUA_mean(1)+1:length(mean_filtered)-locsdMUA_mean(1);
 plot(x1,mean_filtered,'LineWidth',1, 'Color', 'black')
 hold on
  linreg1 = fitlm(locsdMUA_mean(1:3), pksdMUA_mean(1:3), 'y ~ x1');
   linreg_coeff1 = table2array(linreg1.Coefficients(1:2,1));
   pvalues1 = table2array(linreg1.Coefficients(2,4));
   ylinreg1 = linreg_coeff1(2,1) .* locsdMUA_mean + linreg_coeff1(1,1);
 
   hold on
  xreg1 =locsdMUA_mean - locsdMUA_mean(1);
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


xabs = -100:1301;
%x = 1:1702;
nyq = 15000;
channum = 1: length(data_file.good_data);
 
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dMUA = nan(length(xabs), length(channum));
all_locsdMUA_filtered = nan(length(channum),1);
up_dist = nan(1, length(channum));

pvalues = nan(1,length(channum));
mean_linreg = nan(2,1);
clear i ;
 for i = 1:length(channum)
     
     if log_p_layer(i) == 1
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(400:1901,:,:)),2);
   bsl = mean(mean_data(1:200));
   raw_mean_bs(:,i) = mean_data(101:end)- bsl;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpdMUA      = filtfilt(bwb,bwa, raw_mean_bs(:,i));
  
   %{
   [pkssu, locssu] = findpeaks(lpdMUA(50:1201));
x = xabs - locssu(1);
plot(x, lpdMUA)
hold on
%}
    
 filtered_dMUA(:,i) = lpdMUA;
 %find peaks for every channel
 [pksdMUA_filtered, locsdMUA_filtered] = findpeaks(filtered_dMUA(100:1251,i));
 %store first peak location
 if i~=26 
 all_locsdMUA_filtered(i) = locsdMUA_filtered(4);
 %compute the distance between the first peak and the last datapoint and store  
 %in a matrix
 up_dist(:,i)= length(xabs)- all_locsdMUA_filtered(i);
 end
     end
 end
 
 %get the max distance between the first peak and the stimulus onset
 max_low_dist = max(all_locsdMUA_filtered);
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist);
 fp_locked_data = nan(new_dist,length(layer_idx));
 for n = 1:length(layer_idx)
    if ~isnan(all_locsdMUA_filtered(layer_idx(n)))
 lower_bound =max_low_dist-all_locsdMUA_filtered(layer_idx(n))+1;
 upper_bound =max_low_dist-all_locsdMUA_filtered(layer_idx(n))+length(xabs);
 
 fp_locked_data(lower_bound:upper_bound,n) = filtered_dMUA(:,layer_idx(n));
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
 [pksdMUA_mean, locsdMUA_mean] = findpeaks(mean_filtered);
 x1 = -locsdMUA_mean(4)+1:length(mean_filtered)-locsdMUA_mean(4);

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
    layer_idx = find(strcmp(layer, 'M'));

log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
    
  xabs = -100:1301;
%x = 1:1702;
nyq = 15000;
channum = 1: length(log_p_layer);
 
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dMUA = nan(length(xabs), length(channum),4);
all_locsdMUA_filtered = nan(length(channum),4);
up_dist = nan(length(channum), 4);

pvalues = nan(1,length(channum));
mean_linreg = nan(2,1);
 clear peakalign;
for peakalign = 1:4
    clear i ;
 for i = 1:length(log_p_layer)
     
     if log_p_layer(i) == 1
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(400:1901,:,:)),2);
   bsl = mean(mean_data(1:200));
   raw_mean_bs(:,i) = mean_data(101:end)- bsl;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpdMUA      = filtfilt(bwb,bwa, raw_mean_bs(:,i));
  
    
 filtered_dMUA(:,i, peakalign) = lpdMUA;
 %find peaks for every channel
 [pksdMUA_filtered, locsdMUA_filtered] = findpeaks(filtered_dMUA(100:1251,i, peakalign));
 %store first peak location
      if i~=26 
 all_locsdMUA_filtered(i, peakalign) = locsdMUA_filtered(peakalign);
 %compute the distance between the peakalign and the last datapoint and store  
 %in a matrix
 up_dist(i,peakalign)= length(xabs)- all_locsdMUA_filtered(i, peakalign);
      end
     end
 end
end
 %get the max distance between the peakalign and the stimulus onset
 max_low_dist = max(all_locsdMUA_filtered, [], 'all');
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist, [], 'all');
 fp_locked_data = nan(new_dist,length(layer_idx),4);
  %h1 = figure();
 clear peakalign;
for peakalign = 1:4 
   clear n
 for n = 1:length(fp_locked_data(1,:,1))
    if ~isnan(all_locsdMUA_filtered(layer_idx(n),peakalign))
 lower_bound =max_low_dist-all_locsdMUA_filtered(layer_idx(n),peakalign)+1;
 upper_bound =max_low_dist-all_locsdMUA_filtered(layer_idx(n),peakalign)+length(xabs);
 
 fp_locked_data(lower_bound:upper_bound,n, peakalign) = filtered_dMUA(:,layer_idx(n), peakalign);

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
 [pksdMUA_mean, locsdMUA_mean] = findpeaks(mean_filtered);
 x1 = -locsdMUA_mean(peakalign)+1:length(mean_filtered)-locsdMUA_mean(peakalign);

 h = figure();
 
 plot(x1,mean_filtered,'LineWidth',1, 'Color', 'black')
txt1 = 'mean';
text(x1(400), mean_filtered(400), txt1)
 hold on
 %{
 ci_low = mean_filtered - std(fp_locked_data,0,2);
 plot(x1, ci_low, 'LineWidth',1,'Color', 'blue')
 hold on
 
 ci_high = mean_filtered + std(fp_locked_data,0,2);
 plot(x1, ci_high, 'LineWidth',1,'Color', 'blue')
 %}

   title({'DE50_NDE0_su', 'M layers single units', sprintf('aligned on peak %d', peakalign)}, 'Interpreter', 'none')
    xlabel(sprintf('Time from %d peak (ms)', peakalign))
    ylabel('Spike rate (spike/sec)')
    legend( 'mean','Location', 'bestoutside')  
   % filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\plots\',strcat(f{2}, sprintf('_%d_aligned_bsl_mean_M_layer', peakalign)));
   %saveas(gcf, strcat(filename, '.png'));

end

%% compare each peak to every other peak with the aligned data

   %compare peak 1 with peak n
   
  mean_pk1 = mean(fp_locked_data(:,:,1), 2);
     %catch peak 1
   [pks1, locs1] = findpeaks(mean_pk1);
 for peaknb = 2:4
  mean_pk2 = mean(fp_locked_data(:,:,peaknb), 2);
   h = figure();
   univx = -150:150;
   x = locs1(1)-150:locs1(1)+150;
   plot(univx, mean_pk1(x),'LineWidth',1, 'Color', 'black')
   hold on
   %catch peak 2
   [pks2, locs2] = findpeaks(mean_pk2);
   x2 =locs2(peaknb)-150:locs2(peaknb)+150;
   plot(univx, mean_pk2(x2),'LineWidth',1, 'Color', 'red')
   title({'DE50_NDE0_su', 'M layers single units', 'peak aligned'}, 'Interpreter', 'none')
    xlabel('Time from peak (ms)')
    ylabel('Spike rate (spike/sec)')
    legend( 'peak1',sprintf('peak%d', peaknb),'Location', 'bestoutside')  
   filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\plots\',strcat(f{2},sprintf('_aligned_bsl_mean_M_layer_peak%d_peak%d', 1, peaknb )));
   saveas(gcf, strcat(filename, '.png'));
 end
