%this script is intended to plot mean single units activity with both
%unaligned data and aligned data (to the first peak)

%loading the clean data
newdatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
channelfilename = [newdatadir 'clean_norm_SUA_sup_50']; 
data_file = load(channelfilename);


%exclude 160517, (first unit, left empty, it is a K neuron)
%Reject 180806 p1 uclust17, M cell, as doesn't seem well triggered (46)
%Reject 181207 (B) uclust22, M cell, as doesn't seem well triggered (55)
 layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];
 f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
 
 pvaluesdir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
 pvalfilename = [pvaluesdir 'lmer_results_norm.csv'];
 pvalues = dlmread(pvalfilename, ',', 1,1);
 
%% Rough plots of peaks with pvalues (not very representative, as mean unit activity)
  channum = 1: length(data_file.clean_high_SUA);

 %for n = 1:3 
for chan = 1:12:length(channum)
h = figure;
xabs = -199:1300;
idx = [1 3 5 7 9 11 2 4 6 8 10 12];
nyq = 500;
all_mean_data = nan(length(xabs), length(1:12));
clear i ;
 for i = 1:12
 
   mean_data = mean(data_file.clean_high_SUA(chan+i-1).namelist(1:1500,:),2);
   
  
   lpc       = 4.5; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   
  if ~all(isnan(mean_data))
   lpsu      = filtfilt(bwb,bwa, mean_data);
   
    sp = subplot(length(1:6), 2, idx(i));
    plot(xabs, lpsu)
    hold on
    plot([0 0], ylim,'k')
    hold on
    plot([1150 1150], ylim,'k')

    if i == length(6)/2
        ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}Spike Rate (spikes/s)'});
    end
   if i < 6 || (i >= 7)&&(i < 12)
        set(sp, 'XTick', [])
   end
      ylabelh = text(max(xabs), mean(lpsu,1), strcat(num2str(chan+i-1),' | ', layer(chan+i-1)),'HorizontalAlignment','left','FontName', 'Arial','FontSize', 10);
   for npeak = 1:4
         for len = 231:480 %from 200 + 30 (lgn response onset) 
            if lpsu(len) < lpsu(len+1)
   locs = findpeaks(lpsu(len:1450));
        break
            end
         end
         
         if length(locs.loc) >= 4
             %adjust location to the first data point of lpsu (+len), then adjust
             %to xabs (-200)
   xlocation = locs.loc(npeak)+len-200;
         end 
            
   text(xlocation, mean(lpsu,1), strcat(num2str(sprintf('%.2f', pvalues(chan+i-1,npeak)))),'HorizontalAlignment','center','FontName', 'Arial','FontSize', 7);
   end 

 end
 end
      set(gca, 'linewidth',2)
      set(gca,'box','off')

    sgtitle({f{2}, 'all good responses, p<0.05, associated to adaptation pvalues'}, 'Interpreter', 'none')
    xlabel('Time from -50ms from stimulus onset (ms)')
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 8.5 11])
    %filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\plots\',strcat(f{2}, sprintf('x_%d_better_raw_data_peakspvalues_2dec', chan+i-1)));
    %saveas(gcf, strcat(filename, '.png'));

end