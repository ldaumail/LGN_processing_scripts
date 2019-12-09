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

layer_idx = find(strcmp(layer, 'P'));

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
 
%plot the mean data only in the 5Hz range 
 layer_idx = find(strcmp(layer, 'P'));
figure, 
normspec = (nanmean(mean_S(:,:,layer_idx),3) - min(nanmean(mean_S(:,:,layer_idx),3)))./(max(nanmean(mean_S(:,:,layer_idx),3)) - min(nanmean(mean_S(:,:,layer_idx),3)));
plot(tvec,squeeze(normspec(:,1))')
title({'DE50_NDE0_su', 'Mean P layer power at 4Hz vs time normalized', sprintf('')}, 'Interpreter', 'none')
    xlabel('Time from 71ms before stimulus onset(ms)')
    ylabel('Normalized Power at 4Hz(no units)')
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\power_spectrum_plots\',contrast{2},'normalized_power_freq_time_mean_P_layer_4hz');
%saveas(gcf, strcat(filename, '.png')); 
 


reps   = 10000;
all_sigs95 = nan(length(Ses),1);
all_sigs90 = nan(length(Ses),1);
for i = 1:length(Ses)
deresp = nanmean(squeeze(Ses(i).namelist2(1:548,1,:)), 1);
biresp = nanmean(squeeze(Ses(i).namelist2(549:end,1,:)), 1);

    if nanmean(deresp) > nanmean(biresp)
    pref                  = deresp;
    nonpref               = biresp;
    else
    pref                  = biresp;
    nonpref               = deresp;
    end

    [X,Y,T,AUC]           = perfcurve([ones(length(pref),1); repmat(2,length(nonpref),1)],[pref nonpref],1);
    NP                    = length(nonpref);
    PR                    = length(pref);
    catdat                = [pref nonpref];


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

