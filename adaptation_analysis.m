npmkdir    = 'C:\Users\maier\Documents\MATLAB\NPMK-master\'; 
nbanadir   = 'C:\Users\maier\Documents\bootcamp-selected\nbanalysis\'; 

directory  = 'C:\Users\maier\Documents\LGN_data\kacie_preproc\power_channels\';
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

%% Plotting the data

%%Linear regressions
h = figure;
%unitnames = {'190326','190326','190213','190213','190210','190208','190120', ...
    %'190119','190119','190124', 'mean'};
xabs = -50:1301;
xabs2 = 0:1150;
nyq = 15000;

channum = 1: file(1) -2;
norm_mean_percentch = nan(length(xabs), length(channum));
linreg_all = nan(2, length(channum));
pvalues = nan(1,length(channum));
mean_linreg = nan(2,1);
clear i ;
 for i = 1:length(channum)
     
   mean_data = mean(squeeze(data.datafile(i).channel_data.hypo{1,3}.cont3_dMUA_chan(550:1901,:,:)),2);
   norm_mean = (mean_data - min(mean_data))/(max(mean_data)-min(mean_data));
   basedata = norm_mean(25:75);
   mean_bp = mean(basedata,1);
   norm_mean_percentch(:,i) = (norm_mean - mean_bp)*100/mean_bp;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpdMUA      = filtfilt(bwb,bwa, norm_mean_percentch(:,i));
 
   
   %find maxima of the filtered data and plot linear regression
   [pksaMUA, locsaMUA] = findpeaks(lpdMUA(50:1201));
   linreg1 = fitlm(locsaMUA, pksaMUA, 'y ~ x1');
   linreg_all(:,i) = table2array(linreg1.Coefficients(1:2,1));
   pvalues(i) = table2array(linreg1.Coefficients(2,4));
   ylinreg1 = linreg_all(2,i) .* locsaMUA + linreg_all(1,i);
  
plot(locsaMUA, ylinreg1)
%txt1 = [unitnames{i}];
   % text(locsaMUA(1), ylinreg1(1), txt1);
hold on

    if i == 1
    title({'DE50_NDE50_aMUA', 'linreg on maxima: all slopes'}, 'Interpreter', 'none')
    end
 end
 mean_linreg = mean(linreg_all, 2);
 ymeanlin = mean_linreg(2) .*xabs2 + mean_linreg(1);
 hold on
 plot(xabs2,ymeanlin)
 txt1 = ['mean'  ' y = (' num2str(mean_linreg(2)) ')x + (' num2str(mean_linreg(1)) ')'];
 textColor = 'white';
    text(xabs2(200), ymeanlin(200), txt1, 'Color', textColor)
    xlabel('Time from -50ms from stimulus onset (ms)')
    ylabel('% change')
    
%legend('190326','190326','190213','190213','190210','190208','190120', ...
 %   '190119','190119','190124', 'mean')
 
%% mean data time locked on stimulus onset
 
h = figure;
%unitnames = {'190326','190326','190213','190213','190210','190208','190120', ...
    %'190119','190119','190124', 'mean'};
xabs = -50:1301;
nyq = 15000;

channum = 1: file(1) -2;
norm_mean_percentch = nan(length(xabs), length(channum));
mean_resp = nan(length(xabs));
clear i ;
 for i = 1:length(channum)
     
   mean_data = mean(squeeze(data.datafile(i).channel_data.hypo{1,1}.cont1_dMUA_chan(550:1901,:,:)),2);
   norm_mean = (mean_data - min(mean_data))/(max(mean_data)-min(mean_data));
   basedata = norm_mean(25:75);
   mean_bp = mean(basedata,1);
   norm_mean_percentch(:,i) = (norm_mean - mean_bp)*100/mean_bp;
   
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpdMUA      = filtfilt(bwb,bwa, norm_mean_percentch(:,i));
 
   
   %find maxima of the filtered data and plot linear regression
   
plot(xabs, lpdMUA)
%txt1 = [unitnames{i}];
   % text(locsaMUA(1), ylinreg1(1), txt1);
hold on

    if i == 1
    title({'DE0_NDE50_aMUA', 'all responses'}, 'Interpreter', 'none')
    end
 end
 mean_resp = mean(norm_mean_percentch, 2);
 
 hold on
 plot(xabs,mean_resp, 'Color', 'white')
 txt1 = ['mean'];
    text(xabs2(200), ymeanlin(200), txt1, 'Color', 'white')
    xlabel('Time from -50ms from stimulus onset (ms)')
    ylabel('% change')
    
%legend('190326','190326','190213','190213','190210','190208','190120', ...
 %   '190119','190119','190124', 'mean')