gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);

%keepidx = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 30 ...
  %   31 32 33 34 35 36 37 38 39 41 42 43 44 45 46 47 48 49 50 51 53 54 55 56 57 58 59 61 62 ...
  %   64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81];
 contrast = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
 %exclude 160517, (first unit, left empty, it is a K neuron)
 layer = {'','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};

layer_idx = find(strcmp(layer, 'K'));

log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);


%% compute the power spectrum using mtspecgramc function

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
%tvec     = t*1000 + (xabs(1));
time_adj = 1:128;
tvec = cat(2, time_adj , t*1000) ;
%we can also store tvec and f in a struct, but they are all identical
 end
 
%plot the mean data only in the 5Hz range 
 layer_idx = find(strcmp(layer, 'P'));
figure, 
normspec = (nanmean(mean_S(:,:,layer_idx),3) - min(nanmean(mean_S(:,:,layer_idx),3)))./(max(nanmean(mean_S(:,:,layer_idx),3)) - min(nanmean(mean_S(:,:,layer_idx),3)));
time_adj_data = nan(length(tvec),1);
time_adj_data(129:end,1) = normspec(:,1);
plot(tvec,squeeze(time_adj_data(:,1))')
title({'DE50_NDE0_su', 'Mean P layer power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from 71ms before stimulus onset(ms)')
    ylabel('Normalized Power at 4Hz(no units)')
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum_plots\',contrast{2},'normalized_power_freq_time_mean_P_layer_4hz');
%saveas(gcf, strcat(filename, '.png')); 
 
%% perform the Receiver Operating Characteristics analysis

reps   = 10000;
all_sigs95 = nan(length(Ses),1);
all_sigs90 = nan(length(Ses),1);
for i = 1:length(Ses)
part1 = nanmean(squeeze(Ses(i).namelist2(1:547,1,:)), 1);
part2 = nanmean(squeeze(Ses(i).namelist2(548:1122,1,:)), 1);

    if nanmean(part1) > nanmean(part2)
    cond1                  = part1;
    cond2               = part2;
    else
    cond1               = part2;
    cond2               = part1;
    end

    [X,Y,T,AUC]           = perfcurve([ones(length(cond1),1); repmat(2,length(cond2),1)],[cond1 cond2],1);
    NP                    = length(cond2);
    PR                    = length(cond1);
    catdat                = [cond1 cond2];


    for r       = 1:reps
       clear shufNP shufPR
       shufPR         = catdat(randperm(length(catdat),PR));
       shufNP         = catdat(randperm(length(catdat),NP));
       [~,~,~,...
       shufAUC(r)]    = perfcurve([ones(PR,1); repmat(2,NP,1)],[shufPR shufNP],1);
    end

    critT95         = quantile(shufAUC,.95);
    critT90         = quantile(shufAUC,.90);

    if AUC > critT95
       sig95          = 1;
    else
       sig95          = 0;
    end

    all_sigs95(i) = sig95;
    
    if AUC > critT90
       sig90        = 1;
    else
       sig90         = 0;
    end
    all_sigs90(i) = sig90;
end

save( strcat(gooddatadir,'roc_results90.mat'), 'all_sigs90');

sig95_idx = find(all_sigs95);
sig90_idx = find(all_sigs90);

%% Analyze ROC analysis results per layer

gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);


sig95_idx = load( strcat(gooddatadir,'roc_results95.mat'));
 
layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer_idx = find(strcmp(layer, 'M'));
log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);
Ses = struct();
bs_data = struct();
%channum = 1: length(log_p_layer);
%add +128 ms to adjust the time to -100ms before stim onset
mean_S = nan(1147+128,38, length(layer_idx));
Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];
xabs = -100:1301;
xpow = -99:1175;

idx = [1 3 2 4];

clear i
clear chan
for chan =1:4:length(layer_idx)
 figure, 
    for i = 1:4

data = squeeze(data_file.good_data(layer_idx(chan+i-1)).channel_data.hypo{1,2}.cont_su(400:1901,:,:));
   bsl = mean(data(1:200,:));
   norm_mean_bs = nan(length(xabs), 1,length(data(1,:)));
   norm_mean_bs(:,1,:) = data(101:end, :) - bsl;
 
clear S namelist;
[S,t,f]        = mtspecgramc(norm_mean_bs(:,1,:) ,movingwin, params); 
 

mean_S(129:end,:,chan+i-1) = nanmean(S,3);
%tvec     = t*1000 + (xabs(1));
normchan = (mean_S(:,:,chan+i-1) - min(mean_S(:,:,chan+i-1)))./(max(mean_S(:,:,chan+i-1)) - min(mean_S(:,:,chan+i-1)));
xpow = -99:1175; 
sp = subplot(length(1:2), 2, idx(i) );
plot(xpow,normchan(:,1))
    hold on
    plot([0 0], ylim,'k')
    hold on
    plot([1150 1150], ylim,'k')
    xlim([-100 1174])
     if i == length(2)/2
        ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}Normalized Power at 4Hz (no units)'});
     end
    

      ylabelh = text(mean(xpow)/2, max(normchan(:,1))+0.05, strcat(num2str(layer_idx(chan+i-1)),' | ROC sign ', num2str(sig95_idx.all_sigs95(layer_idx(chan+i-1)))),'HorizontalAlignment','left','FontName', 'Arial', 'Interpreter','none','FontSize', 10);
      
      set(gca, 'linewidth',2)
      set(gca,'box','off')
    end
    sgtitle({'DE50_NDE0_su', 'Mean M layer SUA power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
   xlabel('Time from stimulus onset(ms)')
   % ylabel('Normalized Power at 4Hz(no units)')
   filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum_plots\',strcat('DE50_NDE0_', sprintf('mean_sua_power_4hz_M_layer_%d', chan)));
   saveas(gcf, strcat(filename, '.png')); 
end 


%% Proportions of significant units per layer

sig95_idx = load( strcat(gooddatadir,'roc_results95.mat'));
 
%lets remove the first unit, a 'K' unit as it was not well triggered
%(recording started after stimulus onset
layer = {'','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer_idx = find(strcmp(layer, 'K'));
log_p_layer = zeros(length(layer),1);
log_p_layer(layer_idx) = logical(layer_idx);



%dim 2 = channel, dim3 = trials
 Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];

cntdecrease = 0;
cntincrease = 0;
clear i ;
 for i = 1:length(layer_idx)
data = squeeze(data_file.good_data(layer_idx(i)).channel_data.hypo{1,2}.cont_su(400:1901,:,:));
   bsl = mean(data(1:200,:));
   norm_mean_bs = nan(length(xabs), 1,length(data(1,:)));
   norm_mean_bs(:,1,:) = data(101:end, :) - bsl;
   
clear S namelist;
[S,t,f]        = mtspecgramc(norm_mean_bs(:,1,:) ,movingwin, params); 
 

part1 = nanmean(nanmean(squeeze(S(1:547,1,:)), 1),2);

part2 = nanmean(nanmean(squeeze(S(548:1122,1,:)), 1),2);


if part1 > part2 && sig95_idx.all_sigs95(layer_idx(i)) ==1
    cntdecrease = cntdecrease +1;
end

if part1 < part2 && sig95_idx.all_sigs95(layer_idx(i)) ==1
    cntincrease = cntincrease +1;
end

end
 
 percentdecrease = cntdecrease*100./length(layer_idx);
 percentincrease = cntincrease*100./length(layer_idx);
 
 %% plot mean part1 and mean part2 for each cell type
 
 xabs = -100:1301;
  Fs = 1000;
 movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
 params.tapers   = [2 3];
 params.Fs       = Fs;
 params.fpass    = [1 150];


 part1 = nan(length(data_file.good_data),1);
 part2 = nan(length(data_file.good_data),1);
for i = 1:length(part1)
   data = squeeze(data_file.good_data(i).channel_data.hypo{1,2}.cont_su(400:1901,:,:));
   bsl = mean(data(1:200,:));
   norm_mean_bs = nan(length(xabs), 1,length(data(1,:)));
   norm_mean_bs(:,1,:) = data(101:end, :) - bsl;
   
clear S namelist;
[S,t,f]        = mtspecgramc(norm_mean_bs(:,1,:) ,movingwin, params); 
  
normS = (nanmean(S(:,1,:),3) - min(nanmean(S(:,1,:),3)))./(max(nanmean(S(:,1,:),3)) - min(nanmean(S(:,1,:),3)));

    part1(i,1) =nanmean(normS(1:547,1),1);
    part2(i,1) = nanmean(normS(548:1122,1), 1);
end

%parts = cat(2, part1, part2);
%filename = [ gooddatadir, 'part1_part2_norm_power'];
%save(strcat(filename, '.mat'), 'parts')

layer_idx = find(strcmp(layer, 'P'));
d{1} = part1(layer_idx)';
d{2} = part2(layer_idx)';

 %green[167/255 185/255 54/255] [225/255 225/255 129/255]
 %black = [24/255 23/255 23/255] )
 %pink = [229/255, 49/255, 90/255] [1, 119/255, 160/255]
 fig_position = [200 200 600 400]; % coordinates for figures
f4 = figure('Position', fig_position);
subplot(1, 2, 1)
h1 = raincloud_plot(d{1}, 'box_on', 1,'box_dodge', 1, 'box_dodge_amount',...
0, 'dot_dodge_amount', .3, 'color',[229/255, 49/255, 90/255], 'cloud_edge_col', [229/255, 49/255, 90/255]);
title('Part 1 values')
set(gca, 'linewidth',2)
set(gca,'box','off') 
set(gca, 'XLim', [-75 225]);
%K cells [-250 600]);
box off
view([90 -90]);
axis ij
xlabel('Mean power at 4 Hz')
ylabel('Density')
subplot(1, 2, 2)
h2 = raincloud_plot(d{2}, 'box_on', 1, 'box_dodge', 1, 'box_dodge_amount',...
0, 'dot_dodge_amount', .3, 'color', [1, 119/255, 160/255], 'cloud_edge_col', [1, 119/255, 160/255]);
title('Part 2 values')
set(gca, 'linewidth',2)
set(gca,'box','off') 
set(gca, 'XLim', [-75 225]);
box off
view([90 -90]);
axis ij
xlabel('Mean power at 4 Hz')
ylabel('Density')
sgtitle({'P cell class receiver operating characteristics', 'analysis: part 1 vs part 2 mean power'});
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum\plots\',contrast{2},'mean_part1part2_95ci_P_layer_4hz_gathered_pink');
saveas(gcf, strcat(filename, '.svg'));
saveas(gcf, strcat(filename, '.png'));