%{
channeldir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj2\';
filename = [channeldir 'all_norm_data_peaks.mat'];
peak_values = load(filename);

mean_peak_vals = nan(4,length(peak_values.norm_data_peaks));
for i = 1:length(peak_values.norm_data_peaks)
    mean_peak_vals(:,i) = nanmean(peak_values.norm_data_peaks(i).namelist,2);
end

 allfilename = [channeldir 'all_mean_norm_data_peaks'];
 save(strcat(allfilename, '.mat'), 'mean_peak_vals');
 %}
channeldir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj2\';

 gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);


contrast = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
xabs = -199:1300;
nyq = 15000;
channum = 1: length(data_file.good_data);
 
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dSUA = nan(length(xabs), length(channum));
norm_filtered_dSUA= nan(length(xabs), length(channum));
zscore= nan(length(xabs), length(channum));
up_dist = nan(1, length(channum));

clear i ;
 for i = 1:length(channum)
     
    
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(401:1900,:,:)),2);
   bsl = mean(mean_data(1:200));
   raw_mean_bs(:,i) = mean_data(1:end)- bsl;
   variance = (var(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(401:1900,:,:)),0,2));
   zscore(:,i) = raw_mean_bs(:,i) ./variance(:);
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
 
 end
 
 

all_locsdSUA_filtered = nan(4,length(channum));
all_norm_pks = nan(4,length(channum));
all_zscore_pks = nan(4,length(channum));
all_bs_pks = nan(4,length(channum));
 clear i
 for i = 1:length(channum)
     
 %norm_filtered_dSUA(:,i) = (filtered_dSUA(:,i) - min(filtered_dSUA,[],'all'))/(max(filtered_dSUA,[],'all')-min(filtered_dSUA,[],'all'));
 %find peaks for every channel
  for len = 30:750
            if filtered_dSUA(200+len,i) < filtered_dSUA(200+len+1,i)
   locsdSUA_filtered = findpeaks(filtered_dSUA(200+len:1350,i));
        break
            end
  end
         

 %store peaks location
    if length(locsdSUA_filtered.loc)>=4 
    all_locsdSUA_filtered(:,i) = locsdSUA_filtered.loc(1:4)+200+len;
    %all_norm_pks(:,i) = norm_filtered_dSUA(all_locsdSUA_filtered(:,i),i);
    all_bs_pks(:,i) = filtered_dSUA(all_locsdSUA_filtered(:,i),i);
    end
    
 end
 
   allfilename = [channeldir 'all_mean_bs_data_peaks'];
 save(strcat(allfilename, '.mat'), 'all_bs_pks');

 