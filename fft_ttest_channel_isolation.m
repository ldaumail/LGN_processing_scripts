
%load data (loop through the folder)

npmkdir    = 'C:\Users\maier\Documents\MATLAB\NPMK-master\'; 
nbanadir   = 'C:\Users\maier\Documents\bootcamp-selected\nbanalysis\'; 

directory  = 'C:\Users\maier\Documents\LGN_data\kacie_preproc\';
addpath(genpath(directory))
addpath(genpath(npmkdir))
addpath(genpath(nbanadir))

powerdir = 'C:\Users\maier\Documents\LGN_data\kacie_preproc\power\';
channeldir = 'C:\Users\maier\Documents\LGN_data\kacie_preproc\power_channels\';


files = dir(strcat(directory, '*cinterocdrft*'));
for file = 1:length(files)
    
STIMBRdatafile = files(file).name;
STIMfilename   = [directory STIMBRdatafile]; 
powerfilename = [powerdir strcat('power',STIMBRdatafile)];

STIMdMUA = load(strcat(STIMfilename)); 

%preallocate power and frequency, power at 4 hz
power = nan(length(STIMdMUA.STIM.sdftr(1,:,1)),length(STIMdMUA.STIM.sdftr(1,1,:)),1025);
freq = nan(length(STIMdMUA.STIM.sdftr(1,:,1)),length(STIMdMUA.STIM.sdftr(1,1,:)),1025);
fourhzpower =nan(length(STIMdMUA.STIM.sdftr(1,:,1)),length(STIMdMUA.STIM.sdftr(1,1,:)));

%experiment conditions
data_header.contrast1 = STIMdMUA.STIM.contrast ==0 & STIMdMUA.STIM.fixedc >= 0.5;
data_header.contrast2 = STIMdMUA.STIM.contrast >=0.5 & STIMdMUA.STIM.fixedc == 0;
data_header.contrast3 = STIMdMUA.STIM.contrast >=0.5 & STIMdMUA.STIM.fixedc >= 0.5;
data_header.contrast4 = STIMdMUA.STIM.contrast ==0 & STIMdMUA.STIM.fixedc == 0;
f = {'_DE0_NDE50','_DE50_NDE0','_DE50_NDE50'};

%preallocate stats results
signi = nan(length(STIMdMUA.STIM.sdftr(1,:,1)), 3);
pvalue = nan(length(STIMdMUA.STIM.sdftr(1,:,1)), 3);
ci = nan(length(STIMdMUA.STIM.sdftr(1,:,1)),3, 2);
stats = struct();

hypo = {'_DE0_NDE50_vsB','_DE50_NDE0_vsB','_DE50_NDE50_vsB'};

for chan = 1:length(STIMdMUA.STIM.sdftr(1,:,1))
    for tr = 1:length(STIMdMUA.STIM.sdftr(1,1,:))
%compute fft for every trial of a channel, for every channel   
[power(chan,tr,:), freq(chan,tr,:)] = calcFFT(squeeze(STIMdMUA.STIM.sdftr(600:1700,chan,tr)));

%find the index of the frequency vector closest to 4hz and point to the
%power value of this index for every trial, and store the value in
%fourhzpower
[val,index] = min(abs(4-freq(chan,tr,:)));
fourhzpower(chan,tr) = power(chan,tr,index);

    end
    %perform t-test accross :
         %- DE50-NDE0 vs DE0-NDE0
         %- NDE50-DE0 vs DE0-NDE0
         %- NDE50-DE50 vs DE0-NDE0
    cont_power = struct();
    cont_power.cont1_power = fourhzpower(chan,data_header.contrast1(1:length(fourhzpower(chan,:))));
    cont_power.cont2_power = fourhzpower(chan,data_header.contrast2(1:length(fourhzpower(chan,:))));
    cont_power.cont3_power = fourhzpower(chan,data_header.contrast3(1:length(fourhzpower(chan,:))));
    cont_power.cont4_power = fourhzpower(chan,data_header.contrast4(1:length(fourhzpower(chan,:))));
    
    cont_power_cell = struct2cell(cont_power);
    for testnb = 1:3
    X = cont_power_cell{testnb};
    Y = cont_power.cont4_power;
    [signi(chan, testnb), pvalue(chan,testnb)] = ttest2(X,Y);
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
    if signi(chan,testnb) ~= 0
        
        channel_data = struct();
        channel_data.hypo{1}.cont1_dMUA_chan = STIMdMUA.STIM.sdftr(1:length(STIMdMUA.STIM.sdftr(:,1,1)),:,data_header.contrast1(1:length(STIMdMUA.STIM.sdftr(1,chan,:))));
        channel_data.hypo{2}.cont2_dMUA_chan = STIMdMUA.STIM.sdftr(1:length(STIMdMUA.STIM.sdftr(:,1,1)),:,data_header.contrast2(1:length(STIMdMUA.STIM.sdftr(1,chan,:))));
        channel_data.hypo{3}.cont3_dMUA_chan = STIMdMUA.STIM.sdftr(1:length(STIMdMUA.STIM.sdftr(:,1,1)),:,data_header.contrast3(1:length(STIMdMUA.STIM.sdftr(1,chan,:))));
        channel_data.hypo{1}.cont1_power_chan = power(data_header.contrast1(1:length(STIMdMUA.STIM.sdftr(1,chan,:))));
        channel_data.hypo{2}.cont2_power_chan = power(data_header.contrast2(1:length(STIMdMUA.STIM.sdftr(1,chan,:))));
        channel_data.hypo{3}.cont3_power_chan = power(data_header.contrast3(1:length(STIMdMUA.STIM.sdftr(1,chan,:))));
        
        
        channel_data.hypo{1}.cont1_stats_chan.significance = statistics.significance(chan,:);
        channel_data.hypo{2}.cont2_stats_chan.significance = statistics.significance(chan,:);
        channel_data.hypo{3}.cont3_stats_chan.significance = statistics.significance(chan,:);
        channel_data.hypo{1}.cont1_stats_chan.pvalue = statistics.pvalues(chan,:);
        channel_data.hypo{2}.cont2_stats_chan.pvalue = statistics.pvalues(chan,:);
        channel_data.hypo{3}.cont3_stats_chan.pvalue = statistics.pvalues(chan,:);
        %{
        channel_data.hypo{1}.cont1_stats_chan.ConfidenceInterval = statistics.ConfidenceInterval.ci(chan,:);
        channel_data.hypo{2}.cont2_stats_chan.ConfidenceInterval = statistics.ConfidenceInterval.ci(chan,:);
        channel_data.hypo{3}.cont3_stats_chan.ConfidenceInterval = statistics.ConfidenceInterval.ci(chan,:);       
        %}
        
     %change contrast name and data if the significance is NDE50-
     %0DE vs Blank condition
     if testnb == 1
        channel_data.hypo{1}.cont1_dMUA_chan = channel_data.hypo{2}.cont2_dMUA_chan;
        channel_data.hypo{2}.cont2_dMUA_chan = channel_data.hypo{1}.cont1_dMUA_chan;
        channel_data.hypo{1}.cont1_power_chan = channel_data.hypo{2}.cont2_power_chan;
        channel_data.hypo{2}.cont2_power_chan = channel_data.hypo{1}.cont1_power_chan;
        channel_data.hypo{1}.cont1_stats_chan.significance = channel_data.hypo{2}.cont2_stats_chan.significance;
        channel_data.hypo{2}.cont2_stats_chan.significance = channel_data.hypo{1}.cont1_stats_chan.significance;
        channel_data.hypo{1}.cont1_stats_chan.pvalue = channel_data.hypo{2}.cont2_stats_chan.pvalue;
        channel_data.hypo{2}.cont2_stats_chan.pvalue = channel_data.hypo{1}.cont1_stats_chan.pvalue;
        
     end
     STIMBRdatafile(strfind(STIMBRdatafile, '.mat')) = [];
     channelfilename = [channeldir strcat(STIMBRdatafile, sprintf('channel_%d',chan))];
     save(strcat(channelfilename, '.mat'), 'channel_data');

    end
    end
    
    
   
end
    
%save all the power values of each session
save(powerfilename, 'power');


end




% only save all the trials of channels fo which it was significant

         
 
