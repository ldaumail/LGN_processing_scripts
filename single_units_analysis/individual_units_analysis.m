gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);
channeldir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj2\';

keepidx = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 30 ...
     31 32 33 34 35 36 37 38 39 41 42 43 44 45 46 47 48 49 50 51 53 54 55 56 57 58 59 61 62 ...
     64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81];
 
 layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
 f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
 
 %% find peaks in order to analyze them on R and fit a LMER
 data_peaks = struct();
 
 for i = 1:length(keepidx)
 data = squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(600:1901,:,:));
 
 filename = [data_file.good_data(i).channel_data.filename, f{2}];
 filename(strfind(filename, 'mat')) = [];
 filename(strfind(filename, '.')) = []; 
 
 all_locs = nan(4,length(data(1,:)));
 all_pks = nan(4,length(data(1,:)));
 all_lpsu= nan(length(data(:,1)),length(data(1,:)));
 nyq = 15000;
 %all_locs =nan(4,length(channel_data(1,1,:)));
 for n =1:length(data(1,:))
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpsu      = filtfilt(bwb,bwa, data(:,n));
   %all_lpsu(:,n) = lpsu;
  % plot(lpsu)
 if ~all(lpsu==0) 
  for len = 30:550
            if lpsu(len) < lpsu(len+1)
   locs = findpeaks(lpsu(len:1200));
        break
            end
  end
         
         if length(locs.loc) >= 4
             %adjust location to the first data point of lpsu (+len),
    lpsulocs = locs.loc(1:4) + len;
   all_pks(:,n) = lpsu(lpsulocs(1:4));
  
         end 
 end
% plot(locs.loc +len, lpsu(locs.loc +len))
 end
 
 namelist(1,1:length(sprintf('chan_%d',i))) = sprintf('chan_%d',i);
 data_peaks(i).namelist = all_pks;
channelfilename = [channeldir filename];
save(strcat(channelfilename, '.mat'), 'all_pks');
 end
 allfilename = [channeldir 'all_data_peaks'];
 save(strcat(allfilename, '.mat'), 'data_peaks');
 
 
 %% plotting channels with pvalues on 3 last peaks computed on R with LMER with Kenward-Roger approximation

 
 pvaluesdir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj2\lmer_results\';
 pvalfilename = [pvaluesdir 'lmer_results.csv'];
 pvalues = dlmread(pvalfilename, ',', 1,1);
 
  channum = 1: length(data_file.good_data);
  f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 %for n = 1:3 
for chan = 1:12:length(channum)
h = figure;
xabs = -50:1301;
idx = [1 3 5 7 9 11 2 4 6 8 10 12];
nyq = 15000;
all_mean_data = nan(length(xabs), length(1:12));
clear i ;
 for i = 1:12
 
   mean_data = mean(squeeze(data_file.good_data(keepidx(chan+i-1)).channel_data.hypo{1,2}.cont_su(550:1901,:,:)),2);
   
   all_mean_data(:,i) = (mean_data);
  
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpsu      = filtfilt(bwb,bwa, all_mean_data(1:1352,i));
   
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
      ylabelh = text(max(xabs), mean(lpsu,1), strcat(num2str(keepidx(chan+i-1)),' | ', layer(chan+i-1)),'HorizontalAlignment','left','FontName', 'Arial','FontSize', 10);
   for npeak = 1:4
         for len = 80:250
            if lpsu(len) < lpsu(len+1)
   locs = findpeaks(lpsu(len:1200));
        break
            end
         end
         
         if length(locs.loc) >= 4
             %adjust location to the first data point of lpsu (+len), then adjust
             %to xabs (-50)
   xlocation = locs.loc(npeak)+len-50;
         end 
            
   text(xlocation, mean(lpsu,1), strcat(num2str(sprintf('%.2f', pvalues(i,npeak)))),'HorizontalAlignment','center','FontName', 'Arial','FontSize', 7);
   end 

 end
 
      set(gca, 'linewidth',2)
      set(gca,'box','off')

    sgtitle({f{2}, 'all good responses, p<0.05, associated to adaptation pvalues'}, 'Interpreter', 'none')
    xlabel('Time from -50ms from stimulus onset (ms)')
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 8.5 11])
    filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj2\plots\',strcat(f{2}, sprintf('x_%d_better_raw_data_peakspvalues_2dec', keepidx(chan))));
    saveas(gcf, strcat(filename, '.png'));

end


%% Plot 3 weird single units 26 38 64
gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
channeldir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj\';
data_file = load(channelfilename);

keepidx = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 30 ...
     31 32 33 34 35 36 37 38 39 41 42 43 44 45 46 47 48 49 50 51 53 54 55 56 57 58 59 61 62 ...
     64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81];
 
 layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
 f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
h = figure;
chanidx =[26 38 64];
for chan = 1:3
    
xabs = -50:1301;
%idx = [1 2 3];
nyq = 15000;
all_mean_data = nan(length(xabs), length(1:3));

    mean_data = mean(squeeze(data_file.good_data(chanidx(chan)).channel_data.hypo{1,2}.cont_su(550:1901,:,:)),2);
   
   all_mean_data(:,chan) = (mean_data);
  
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpsu      = filtfilt(bwb,bwa, all_mean_data(1:1352,chan));
   
    sp = subplot(length(1:3), 1, chan);
    plot(xabs, lpsu)
    hold on
    plot([0 0], ylim,'k')
    hold on
    plot([1150 1150], ylim,'k')

    if i == length(4)/2
        ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}Spike Rate (spikes/s)'});
    end
   if i < 3
        set(sp, 'XTick', [])
   end
      ylabelh = text(max(xabs), mean(lpsu,1), strcat(num2str(chanidx(chan)),' | ', layer(chanidx(chan))),'HorizontalAlignment','left','FontName', 'Arial','FontSize', 10);
   for npeak = 1:4
         for len = 80:250
            if lpsu(len) < lpsu(len+1)
   locs = findpeaks(lpsu(len:1200));
        break
            end
         end
         
         if length(locs.loc) >= 4
             %adjust location to the first data point of lpsu (+len), then adjust
             %to xabs (-50)
   xlocation = locs.loc(npeak)+len-50;
         end 
            
   text(xlocation, mean(lpsu,1), strcat(num2str(sprintf('%.2f', pvalues(chanidx(chan),npeak)))),'HorizontalAlignment','center','FontName', 'Arial','FontSize', 7);
   end 

 
 
      set(gca, 'linewidth',2)
      set(gca,'box','off')

    sgtitle({f{2}, 'good responses, p<0.05, associated to adaptation pvalues'}, 'Interpreter', 'none')
    xlabel('Time from -50ms from stimulus onset (ms)')
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 8.5 11])
    %filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj\plots\',strcat(f{2}, sprintf('x_%d_better_raw_data_peakspvalues_2dec', keepidx(chan))));
    %saveas(gcf, strcat(filename, '.png'));

end

%% compute proportion of significant adaptation per peak and proportion of neurons adapting for a certain amount of 
%peak from peak 2 to 4
%%here we exclude units 26 38 and 64

gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
channeldir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj2\';
data_file = load(channelfilename);
pvaluesdir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_peakadj2\lmer_results\';
 pvalfilename = [pvaluesdir 'lmer_results.csv'];
 pvalues = dlmread(pvalfilename, ',', 1,1);
   % only apply this line once (exclude 26, 38,64 with indices 26, 38 and 57
 % respectively)
 pvalues([25,35,57],:) = [];
 %                        %
  %                     %
%keepidx2 doesn't include 26, 38, 64
keepidx2 = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 30 ...
     31 32 33 34 35 36 37 39 41 42 43 44 45 46 47 48 49 50 51 53 54 55 56 57 58 59 61 62 ...
     65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81];
 %layer2 doesn't include 26, 38, 64
 layer2 = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','M','','P', ...
'P','','','K','P','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};

 layer_idx = find(strcmp(layer2, 'K'));

log_p_layer = zeros(length(layer2),1);
log_p_layer(layer_idx) = logical(layer_idx);

 f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
 

 nyq = 15000;
 all_locs = nan(4,length(layer_idx));
 all_pks = nan(4,length(layer_idx));
 
 channum = 1: length(layer_idx);
 
 for i = 1:length(channum)
 
 mean_data = mean(squeeze(data_file.good_data(layer_idx(i)).channel_data.hypo{1,2}.cont_su(550:1901,:,:)),2);
   
   all_mean_data(:,i) = (mean_data);
  
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpsu      = filtfilt(bwb,bwa, all_mean_data(1:1352,i));
   
     for len = 80:250
            if lpsu(len) < lpsu(len+1)
   locs = findpeaks(lpsu(len:1200));
        break
            end
     end
         
         if length(locs.loc) >= 4
             %adjust location to the first data point of lpsu (+len),
    lpsulocs = locs.loc(1:4) + len;
   all_pks(:,i) = lpsu(lpsulocs(1:4));
  
         end 
 end
 
 cntpk2 = 0;
 cntpk3 = 0;
 cntpk4 = 0;
 cntpk2pk3 = 0;
 cntpk2pk3pk4 = 0;
 
 for nunit = 1:length(layer_idx)
     if all_pks(2,nunit) < all_pks(1,nunit) && pvalues(layer_idx(nunit),2) < .05
         cntpk2 = cntpk2 +1;
     end
     if all_pks(3,nunit) < all_pks(1,nunit) && pvalues(layer_idx(nunit),3) < .05
         cntpk3 = cntpk3 +1;
     end
     if all_pks(4,nunit) < all_pks(1,nunit) && pvalues(layer_idx(nunit),4) < .05
         cntpk4 = cntpk4 +1;
     end
     if all_pks(2,nunit) < all_pks(1,nunit) && pvalues(layer_idx(nunit),2) < .05 && ...
             all_pks(3,nunit) < all_pks(1,nunit) && pvalues(layer_idx(nunit),3) < .05
         cntpk2pk3 = cntpk2pk3 +1;
     end
     if all_pks(2,nunit) < all_pks(1,nunit) && pvalues(layer_idx(nunit),2) < .05 && ...
             all_pks(3,nunit) < all_pks(1,nunit) && pvalues(layer_idx(nunit),3) < .05 ...
             && all_pks(4,nunit) < all_pks(1,nunit) && pvalues(layer_idx(nunit),4) < .05
         cntpk2pk3pk4 = cntpk2pk3pk4 +1;
     end
 end
 
 percentpk2 = cntpk2*100/length(all_pks(1,:));
 percentpk3 = cntpk3*100/length(all_pks(1,:)); 
 percentpk4 = cntpk4*100/length(all_pks(1,:));
 percentpk2pk3 = cntpk2pk3*100/length(all_pks(1,:));
 percentpk2pk3pk4 = cntpk2pk3pk4*100/length(all_pks(1,:));
 
 %% Replicate the analysis on the troughs instead of the peaks
 
 gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);
channeldir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_troughadj\';


keepidx = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 30 ...
     31 32 33 34 35 36 37 38 39 41 42 43 44 45 46 47 48 49 50 51 53 54 55 56 57 58 59 61 62 ...
     64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81];
 
 layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
 f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
 
  %% find troughs in order to analyze them on R and fit a LMER
 data_troughs = struct();
 
 for i = 1:length(keepidx)
 data = squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(600:1901,:,:));
 
 filename = [data_file.good_data(i).channel_data.filename, f{2}];
 filename(strfind(filename, 'mat')) = [];
 filename(strfind(filename, '.')) = []; 
 
 all_locs = nan(3,length(data(1,:)));
 all_trghs = nan(3,length(data(1,:)));
 all_lpsu= nan(length(data(:,1)),length(data(1,:)));
 nyq = 15000;
 %all_locs =nan(4,length(channel_data(1,1,:)));
 for n =1:length(data(1,:))
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpsu      = filtfilt(bwb,bwa, data(:,n));
   all_lpsu(:,n) = lpsu;
   %plot(lpsu)
   %hold on
   
  if ~all(lpsu==0)
    for len = 1:400
   
            if lpsu(len) < lpsu(len+1)
   locsp = findpeaks(lpsu(len:1200));
    
   locst = findpeaks(-lpsu(len+locsp.loc(1):1200));
       break
            end
    end
        
         if length(locst.loc) >= 3
             %adjust location to the first data point of lpsu (+len),
    lpsulocs = locst.loc(1:3) +locsp.loc(1) + len;
   all_trghs(:,n) = lpsu(lpsulocs(1:3));
  
         end 
  end
     %plot(locsp.loc, lpsu(locsp.loc))
     %hold on
     %plot(locst.loc + locsp.loc(1) +len, lpsu(locst.loc +locsp.loc(1)+len))

end
 
 namelist(1,1:length(sprintf('chan_%d',i))) = sprintf('chan_%d',i);
 data_troughs(i).namelist = all_trghs;
channelfilename = [channeldir filename];
save(strcat(channelfilename, '.mat'), 'all_trghs');
 end
 allfilename = [channeldir 'all_data_troughs'];
 save(strcat(allfilename, '.mat'), 'data_troughs');
 
  %% plotting channels with pvalues on 2 last troughs computed on R with LMER with Kenward-Roger approximation

 
 pvaluesdir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_troughadj\lmer_results\';
 pvalfilename = [pvaluesdir 'lmer_results.csv'];
 pvalues = dlmread(pvalfilename, ',', 1,1);
 
  channum = 1: length(data_file.good_data);
  f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 %for n = 1:3 
for chan = 1:12:length(channum)
h = figure;
xabs = -50:1301;
idx = [1 3 5 7 9 11 2 4 6 8 10 12];
nyq = 15000;
all_mean_data = nan(length(xabs), length(1:12));
clear i ;
 for i = 1:12
 
   mean_data = mean(squeeze(data_file.good_data(keepidx(chan+i-1)).channel_data.hypo{1,2}.cont_su(550:1901,:,:)),2);
   
   all_mean_data(:,i) = (mean_data);
  
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpsu      = filtfilt(bwb,bwa, all_mean_data(1:1352,i));
   
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
      ylabelh = text(max(xabs), mean(lpsu,1), strcat(num2str(keepidx(chan+i-1)),' | ', layer(chan+i-1)),'HorizontalAlignment','left','FontName', 'Arial','FontSize', 10);
   for ntrough = 1:3
         for len = 80:250
            if lpsu(len) < lpsu(len+1)
   locsp = findpeaks(lpsu(len:1200));
   locst = findpeaks(-lpsu(len+locsp.loc(1):1200));
        break
            end
         end
         
         if length(locst.loc) >= 3
             %adjust location to the first data point of lpsu (+len), then adjust
             %to xabs (-50)
   xlocation = locst.loc(ntrough)+len +locsp.loc(1)-50;
         end 
            
   text(xlocation, mean(lpsu,1), strcat(num2str(sprintf('%.2f', pvalues(i,ntrough)))),'HorizontalAlignment','center','FontName', 'Arial','FontSize', 7);
   end 

 end
 
      set(gca, 'linewidth',2)
      set(gca,'box','off')

    sgtitle({f{2}, 'all good responses, p<0.05, associated to adaptation pvalues on troughs'}, 'Interpreter', 'none')
    xlabel('Time from -50ms from stimulus onset (ms)')
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 8.5 11])
    filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_troughadj\plots\',strcat(f{2}, sprintf('x_%d_better_raw_data_troughspvalues_2dec', keepidx(chan))));
    saveas(gcf, strcat(filename, '.png'));

end
 
%% compute proportion of significant adaptation per trough and proportion of neurons adapting for a certain amount of 
%troughs from trough 2 to 4
%%here we exclude units 26 38 and 64

gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
channeldir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_troughadj\';
data_file = load(channelfilename);
pvaluesdir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\individual_channels_troughadj\lmer_results\';
 pvalfilename = [pvaluesdir 'lmer_results.csv'];
 pvalues = dlmread(pvalfilename, ',', 1,1);
   % only apply this line once (exclude 26, 38,64 with indices 26, 38 and 57
 % respectively)
 pvalues([25,35,57],:) = [];
 %                        %
  %                     %
%keepidx2 doesn't include 26, 38, 64
keepidx2 = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 30 ...
     31 32 33 34 35 36 37 39 41 42 43 44 45 46 47 48 49 50 51 53 54 55 56 57 58 59 61 62 ...
     65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81];
 %layer2 doesn't include 26, 38, 64
 layer2 = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','M','','P', ...
'P','','','K','P','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};

 layer_idx = find(strcmp(layer2, 'M'));



 f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
 

 nyq = 15000;
 all_locs = nan(3,length(layer_idx));
 all_trghs = nan(3,length(layer_idx));
 
 channum = 1: length(layer_idx);
 
 for i = 1:length(channum)
 
 mean_data = mean(squeeze(data_file.good_data(layer_idx(i)).channel_data.hypo{1,2}.cont_su(550:1901,:,:)),2);
   
   all_mean_data(:,i) = (mean_data);
  
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpsu      = filtfilt(bwb,bwa, all_mean_data(1:1352,i));
   
     for len = 80:250
            if lpsu(len) < lpsu(len+1)
   locsp = findpeaks(lpsu(len:1200));
   locst = findpeaks(-lpsu(len+locsp.loc(1):1200));
        break
            end
     end
         
         if length(locst.loc) >= 3
             %adjust location to the first data point of lpsu (+len),
    lpsulocs = locst.loc(1:3) + len +locsp.loc(1);
 
   all_trghs(:,i) = lpsu(lpsulocs(1:3));
  
         end 
 end
 
 cnttrgh2 = 0;
 cnttrgh3 = 0;
 cnttrgh2trgh3 = 0;
 
 
 for nunit = 1:length(layer_idx)
     if all_trghs(2,nunit) < all_trghs(1,nunit) && pvalues(layer_idx(nunit),2) < .05
         cnttrgh2 = cnttrgh2 +1;
     end
     if all_trghs(3,nunit) < all_trghs(1,nunit) && pvalues(layer_idx(nunit),3) < .05
         cnttrgh3 = cnttrgh3 +1;
     end
    
     if all_trghs(2,nunit) < all_trghs(1,nunit) && pvalues(layer_idx(nunit),2) < .05 && ...
             all_trghs(3,nunit) < all_trghs(1,nunit) && pvalues(layer_idx(nunit),3) < .05
         cnttrgh2trgh3 = cnttrgh2trgh3 +1;
     end
     
 end
 
 percenttrgh2 = cnttrgh2*100/length(all_trghs(1,:));
 percenttrgh3 = cnttrgh3*100/length(all_trghs(1,:)); 
 percenttrgh2trgh3 = cnttrgh2trgh3*100/length(all_trghs(1,:));

