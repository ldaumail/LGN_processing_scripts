gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);

keepidx = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 30 ...
     31 32 33 34 35 36 37 38 39 41 42 43 44 45 46 47 48 49 50 51 53 54 55 56 57 58 59 61 62 ...
     64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81];
 f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
 

 xabs = -200:1301;
idx = [1 3 2 4];
channum = 1: length(data_file.good_data);
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dSUA = nan(length(xabs), length(channum));


nyq = 15000;
norm_mean_data = nan(length(xabs), length(1:4));
clear i ;
for chan = 1:4:length(channum)
   h = figure;
 for i = 1:4
 pvalue = data_file.good_data(chan+i-1).channel_data.hypo{1,2}.cont_stats.pvalue;

    if pvalue <= 0.05 && ~isnan(pvalue)
   mean_data = mean(squeeze(data_file.good_data(chan+i-1).channel_data.hypo{1,2}.cont_su(400:1901,:,:)),2);
   
   norm_mean_data(:,i) = (mean_data);
  %{
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpdSUA      = filtfilt(bwb,bwa, norm_mean_data(:,i));
   %}
    sp = subplot(length(1:2), 2, idx(i) );
    plot(xabs, mean_data)
    hold on
    plot([0 0], ylim,'k')
    hold on
    plot([1150 1150], ylim,'k')

    if i == length(2)/2
        ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}Spike Rate (spikes/s)'});
    end
   if i < 2 || (i >= 3)&&(i < 4)
        set(sp, 'XTick', [])
   end
   channame = data_file.good_data(chan+i-1).channel_data.filename;
      %ylim([-2 2]); 
      ylabelh = text(max(xabs), mean(mean_data,1), strcat(num2str(chan+i-1),' | ', num2str(channame)),'HorizontalAlignment','left','FontName', 'Arial', 'Interpreter','none','FontSize', 10);
      
      
      set(gca, 'linewidth',2)
      set(gca,'box','off')
     % h = subplot(1,1,1); 
     %set(h,'position',get(h,'position').*[1 1 1 1.2]);
    end
 end
    sgtitle({f{2}, 'all good responses, p<0.05'}, 'Interpreter', 'none')
    xlabel('Time from -50ms from stimulus onset (ms)')
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 8.5 11])
   filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\plots\',strcat(f{2}, sprintf('trigger_raw_data_%d', chan)));
  % saveas(gcf, strcat(filename, '.png'));
end  

%% plot a nice unit


xabs = -200:1301;

filtered_dSUA = nan(length(xabs), length(channum));

nyq = 15000;
   mean_data = mean(squeeze(data_file.good_data(23).channel_data.hypo{1,2}.cont_su(400:1901,:,:)),2);
 %
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpdSUA      = filtfilt(bwb,bwa, mean_data);
%}
   h = figure();
    plot(xabs, lpdSUA)
    hold on
    plot([0 0], ylim,'k')
    hold on
    plot([1150 1150], ylim,'k')
       set(gca, 'linewidth',2)
      set(gca,'box','off')
    sgtitle('Low pass filtered single unit time series', 'Interpreter', 'none')
    xlabel('Time from stimulus onset (ms)')
    ylabel({'\fontsize{9}Contacts','\fontsize{9}Spike Rate (spikes/s)'})
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 8.5 11])
   
   
  %%
   mean_data = mean(squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(300:800,:,:)),2);
   
   raw_mean_bs(:,i) = mean_data;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpdSUA      = filtfilt(bwb,bwa, raw_mean_bs(:,i));
  
   %{
   locssu = findpeaks(lpdMUA(50:1201));
x = xabs - locssuloc(1);
plot(x, lpdMUA)
hold on
%}
    
 filtered_dSUA(:,i) = lpdSUA;
 %find peaks for every channel that is not accounting for any maxima in the
 %beginning of the trace
 
 plot(xabs, raw_mean_bs(:,i))
 hold on
 end
 
%{ 
  mean_filtered = mean(fp_locked_data, 2);
 locsdMUA_mean = findpeaks(mean_filtered(1:1586));
 x1 = 1:length(mean_filtered);
 plot(x1,mean_filtered,'LineWidth',1, 'Color', 'black')
txt1 = 'mean';
text(x1(400), mean_filtered(400), txt1)
 hold on
 
 median = nanmedian(fp_locked_data, 2);
 median_filtered      = filtfilt(bwb,bwa, median);
 locsdMUA_median = findpeaks(median_filtered(1:1586));
 x2 = 1:length(mean_filtered); 
 plot(x2,median_filtered,'LineWidth',1, 'Color', 'red')

txt2 = 'median';
 textColor = 'red';
    text(x2(680), median_filtered(680), txt2, 'Color', textColor)
    %}
   title({'DE50_NDE0_su', 'all SUA baselines of good_data'}, 'Interpreter', 'none')
    xlabel('Time from stimulus onset (ms)')
    ylabel('Spike rate (spike/sec)')