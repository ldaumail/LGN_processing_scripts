
npmkdir    = 'C:\Users\maier\Documents\MATLAB\NPMK-master\'; 
nbanadir   = 'C:\Users\maier\Documents\bootcamp-selected\nbanalysis\'; 

directory  = 'C:\Users\maier\Documents\LGN_data\single_units\';
addpath(genpath(directory))
addpath(genpath(npmkdir))
addpath(genpath(nbanadir))

powerdir = 'C:\Users\maier\Documents\LGN_data\single_units\power\';
channeldir = 'C:\Users\maier\Documents\LGN_data\single_units\power_channels\';
invertedchanneldir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\';


files = dir(strcat(directory, '*cinterocdrft*'));
for file = 1:length(files)
    
STIMBRdatafile = files(file).name;

STIMfilename   = [directory STIMBRdatafile]; 
powerfilename = [powerdir strcat('power',STIMBRdatafile)];

STIMdMUA = load(strcat(STIMfilename)); 

%preallocate power and frequency, power at 4 hz
power = nan(length(STIMdMUA.STIM.sdftr(1,:)),1025);
freq = nan(length(STIMdMUA.STIM.sdftr(1,:)),1025);
fourhzpower =nan(length(STIMdMUA.STIM.sdftr(1,:)));

%experiment conditions
data_header.contrast1 = STIMdMUA.STIM.contrast ==0 & STIMdMUA.STIM.fixedc >= 0.5;
data_header.contrast2 = STIMdMUA.STIM.contrast >=0.5 & STIMdMUA.STIM.fixedc == 0;
data_header.contrast3 = STIMdMUA.STIM.contrast >=0.5 & STIMdMUA.STIM.fixedc >= 0.5;
data_header.contrast4 = STIMdMUA.STIM.contrast ==0 & STIMdMUA.STIM.fixedc == 0;
f = {'_DE0_NDE50','_DE50_NDE0','_DE50_NDE50'};

%preallocate stats results
signi = nan(1,3);
pvalue = nan(1, 3);
ci = nan(1,3, 2);
stats = struct();

hypo = {'_DE0_NDE50_vsB','_DE50_NDE0_vsB','_DE50_NDE50_vsB'};

%for chan = 1:length(STIMdMUA.STIM.sdftr(1,:,1))
    for tr = 1:length(STIMdMUA.STIM.sdftr(1,:))
%compute fft for every trial of a channel, for every channel   
[power(tr,:), freq(tr,:)] = calcFFT(squeeze(STIMdMUA.STIM.sdftr(600:1700,tr)));

%find the index of the frequency vector closest to 4hz and point to the
%power value of this index for every trial, and store the value in
%fourhzpower
[val,index] = min(abs(4-freq(tr,:)));
fourhzpower(tr) = power(tr,index);

    end
    %perform t-test accross :
         %- DE50-NDE0 vs DE0-NDE0
         %- NDE50-DE0 vs DE0-NDE0
         %- NDE50-DE50 vs DE0-NDE0
    cont_power = struct();
    cont_power.cont1_power = fourhzpower(data_header.contrast1(1:length(fourhzpower)));
    cont_power.cont2_power = fourhzpower(data_header.contrast2(1:length(fourhzpower)));
    cont_power.cont3_power = fourhzpower(data_header.contrast3(1:length(fourhzpower)));
    cont_power.cont4_power = fourhzpower(data_header.contrast4(1:length(fourhzpower)));
    
    cont_power_cell = struct2cell(cont_power);
    for testnb = 1:3
    X = cont_power_cell{testnb};
    Y = cont_power.cont4_power;
    [signi(testnb), pvalue(testnb)] = ttest2(X,Y);
    statistics = struct();
    statistics.significance = signi;
    statistics.pvalues = pvalue;
   
    %{
    [signi(chan, testnb), pvalue(chan,testnb), ci(chan,testnb,:), stats] = ttest2(X,Y);
  
    all_stats.strcat('stats', sprintf('channel %d',chan), hypo{testnb}) = stats;
    statistics = struct();
    statistics.significance = signi;
    statistics.pvalues = pvalue;
    statistics.ConfidenceInterval = ci;
    statistics.teststats = stats;
    %}
    %if test result is significant, save the channel dMUA data and power
    %data
    if signi(testnb) ~= 0
        
        channel_data = struct();
        channel_data.hypo{1}.cont_su = STIMdMUA.STIM.sdftr(:,data_header.contrast1(1:length(STIMdMUA.STIM.sdftr(1,:))));
        channel_data.hypo{2}.cont_su = STIMdMUA.STIM.sdftr(:,data_header.contrast2(1:length(STIMdMUA.STIM.sdftr(1,:))));
        channel_data.hypo{3}.cont_su = STIMdMUA.STIM.sdftr(:,data_header.contrast3(1:length(STIMdMUA.STIM.sdftr(1,:))));
        channel_data.hypo{1}.cont_power = power(data_header.contrast1(1:length(STIMdMUA.STIM.sdftr(1,:))));
        channel_data.hypo{2}.cont_power = power(data_header.contrast2(1:length(STIMdMUA.STIM.sdftr(1,:))));
        channel_data.hypo{3}.cont_power = power(data_header.contrast3(1:length(STIMdMUA.STIM.sdftr(1,:))));
        
        
        channel_data.hypo{1}.cont_stats.significance = statistics.significance(:);
        channel_data.hypo{2}.cont_stats.significance = statistics.significance(:);
        channel_data.hypo{3}.cont_stats.significance = statistics.significance(:);
        channel_data.hypo{1}.cont_stats.pvalue = statistics.pvalues(1);
        channel_data.hypo{2}.cont_stats.pvalue = statistics.pvalues(2);
        channel_data.hypo{3}.cont_stats.pvalue = statistics.pvalues(3);
        channel_data.onsets = STIMdMUA.STIM.onsets;
        channel_data.offsets = STIMdMUA.STIM.offsets;
        channel_data.trstart = STIMdMUA.STIM.trstart;
        channel_data.trend = STIMdMUA.STIM.trend;
        channel_data.presnum = STIMdMUA.STIM.presnum;
        channel_data.trial = STIMdMUA.STIM.trial;
        channel_data.tilt = STIMdMUA.STIM.tilt;
        channel_data.sf = STIMdMUA.STIM.sf;
        channel_data.contrast = STIMdMUA.STIM.contrast;
        channel_data.fixedc = STIMdMUA.STIM.fixedc;
        channel_data.diameter = STIMdMUA.STIM.diameter;
        channel_data.eye = STIMdMUA.STIM.eye;
        channel_data.oridist = STIMdMUA.STIM.oridist;
        channel_data.phase = STIMdMUA.STIM.phase;
        channel_data.temporal_freq = STIMdMUA.STIM.temporal_freq;
        channel_data.xpos = STIMdMUA.STIM.xpos;
        channel_data.ypos = STIMdMUA.STIM.ypos;
        if exist('STIMdMUA.STIM.fix_x')
       channel_data.fix_x = STIMdMUA.STIM.fix_x;
        end
        if exist('STIMdMUA.STIM.fix_y')
        channel_data.fix_y = STIMdMUA.STIM.fix_y;
        end
        channel_data.photo_on = STIMdMUA.STIM.photo_on;
        channel_data.trg_photo = STIMdMUA.STIM.trg_photo;
        channel_data.refresh = STIMdMUA.STIM.refresh;
        channel_data.measured_refresh = STIMdMUA.STIM.measured_refresh;
        channel_data.paradigm = STIMdMUA.STIM.paradigm;
        channel_data.sdftr_chan = STIMdMUA.STIM.sdftr(:,:);
        channel_data.spk_bin_chan = STIMdMUA.STIM.spk_bin(:,:);
        channel_data.wf = STIMdMUA.STIM.wf;
        if exist('STIMdMUA.STIM.wfWin')
        channel_data.wfWin = STIMdMUA.STIM.wfWin;
        end
        channel_data.clust = STIMdMUA.STIM.clust;
        if exist('STIMdMUA.STIM.BRdatafile')
        channel_data.BRdatafile = STIMdMUA.STIM.BRdatafile;
        end
        channel_data.norm_autocorr = STIMdMUA.STIM.norm_autocorr;
        channel_data.lags = STIMdMUA.STIM.lags;
        channel_data.autocorrFrac = STIMdMUA.STIM.autocorrFrac;
        channel_data.chan = STIMdMUA.STIM.chan;
        
        %{
        channel_data.hypo{1}.cont1_stats_chan.ConfidenceInterval = statistics.ConfidenceInterval.ci(chan,:);
        channel_data.hypo{2}.cont2_stats_chan.ConfidenceInterval = statistics.ConfidenceInterval.ci(chan,:);
        channel_data.hypo{3}.cont3_stats_chan.ConfidenceInterval = statistics.ConfidenceInterval.ci(chan,:);       
        %}
        
     %change contrast name and data if the significance is NDE50-
     %0DE vs Blank condition as it means the DE and NDE are inverted
     
     if testnb == 1
        channel_data.hypo{1}.cont_su = channel_data.hypo{2}.cont_su;
        channel_data.hypo{2}.cont_su = channel_data.hypo{1}.cont_su;
        channel_data.hypo{1}.cont_power = channel_data.hypo{2}.cont_power;
        channel_data.hypo{2}.cont_power = channel_data.hypo{1}.cont_power;
        channel_data.hypo{1}.cont_stats.significance = channel_data.hypo{2}.cont_stats.significance;
        channel_data.hypo{2}.cont_stats.significance = channel_data.hypo{1}.cont_stats.significance;
        channel_data.hypo{1}.cont_stats.pvalue = channel_data.hypo{2}.cont_stats.pvalue;
        channel_data.hypo{2}.cont_stats.pvalue = channel_data.hypo{1}.cont_stats.pvalue;
        channel_data.fixedc = channel_data.contrast;
        channel_data.contrast = channel_data.fixedc;
       
     end
    
     STIMBRdatafile(strfind(STIMBRdatafile, '.mat')) = [];
     channelfilename = [invertedchanneldir STIMBRdatafile];
     save(strcat(channelfilename, '.mat'), 'channel_data');

    end
    end 
%end
    
%save all the power values of each session
save(powerfilename, 'power');
end
