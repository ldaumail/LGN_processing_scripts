npmkdir    = 'C:\Users\maier\Documents\MATLAB\NPMK-master\'; 
nbanadir   = 'C:\Users\maier\Documents\bootcamp-selected\nbanalysis\'; 

directory  = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\';
addpath(genpath(directory))
addpath(genpath(npmkdir))
addpath(genpath(nbanadir))

files = dir(directory);
data = struct();

%BRdatafile = cell(1,365);

for file = 3:length(files)
%BRdatafile{file-2} = {files(file).name};
STIMfilename   = [directory files(file).name]; 
data.datafile(file-2) = load(strcat(STIMfilename)); 
end

%%Plot channels 
 channum = 1: file(1) -2;
 
  f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};

 %for n = 1:3
 for chan = 1:24:length(channum)
h = figure;
xabs = -50:1301;
idx = [1 3 5 7 9 11 13 15 17 19 21 23 2 4 6 8 10 12 14 16 18 20 22 24];
nyq = 15000;
norm_mean_percentch = nan(length(xabs), length(1:24));
clear i ;
 for i = 1:24
     pvalue = data.datafile(chan+i-1).channel_data.hypo{1,2}.cont_stats.pvalue;
  
    if pvalue <= 0.05 && ~isnan(pvalue)

   mean_data = mean(squeeze(data.datafile(chan+i-1).channel_data.hypo{1,2}.cont_su(550:1901,:,:)),2);
   %{
   basedata = mean_data(25:75);
    mean_bp = mean(basedata,1);
   norm_mean_percentch(:,i) = (mean_data - mean_bp)/(mean_bp);
  %}
    sp = subplot(length(1:12), 2, idx(i));
    plot(xabs, mean_data)
    hold on
    plot([0 0], ylim,'k')
    hold on
    plot([1150 1150], ylim,'k')

    if i == length(12)/2
        ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}(no unit)'});
    end
   if i < 12 || (i >= 13)&&(i < 24)
        set(sp, 'XTick', [])
   end
   pvalue2 = data.datafile(chan+i-1).channel_data.hypo{1,2}.cont_stats.pvalue;
      %ylim([-2 2]); 
      ylabelh = text(max(xabs), max(ylim)/2, strcat(num2str(chan+i-1),' | ', num2str(pvalue2)),'HorizontalAlignment','left','FontName', 'Arial','FontSize', 10);
      
      
      set(gca, 'linewidth',2)
      set(gca,'box','off')
     % h = subplot(1,1,1); 
     %set(h,'position',get(h,'position').*[1 1 1 1.2]);
     end
   end
    sgtitle({f{2}, 'all responses'}, 'Interpreter', 'none')
    xlabel('Time from -50ms from stimulus onset (ms)')
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 8.5 11])
    filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\plots\',strcat(f{2}, sprintf('raw_data_%d_pvalue05', chan)));
   saveas(gcf, strcat(filename, '.png'));
 end
 
 %% plot power as a function of the frequency to spot the channels that don't show a peak at 4hz
 
 mean_data = mean(squeeze(data.datafile(28).channel_data.hypo{1,2}.cont_su(550:1901,:,:)),2);
 
 [power, freq] = calcFFT(squeeze(mean_data(50:1150)));
h = figure();
 plot(freq, power)
 xlim([0 10])
 
 %% save the good data in a different file
 
 keepidx = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 30 ...
     31 32 33 34 35 36 37 38 39 41 42 43 44 45 46 47 48 49 50 51 53 54 55 56 57 58 59 61 62 ...
     64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81];
 
 good_data = data.datafile(keepidx);
   for i = 1:length(keepidx)
    filename = files(keepidx(i)+2).name;
   good_data(i).channel_data.filename = filename;
   end
   invertedchanneldir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
  
   channelfilename = [invertedchanneldir 'good_single_units_data_4bmpmore']; 
   save(strcat(channelfilename, '.mat'), 'good_data');
 save(strcat(channelfilename, '.mat'), 'good_data', '-v7.3');
 
 