if strcmp(getenv('USER'),'maierav')
    npmkdir    = '/Users/alex 1/Desktop/LAB/Loic/NPMK-master/'; 
    nbanadir   = '/Users/alex 1/Desktop/LAB/Loic/nbanalysis/'; 
 
    directory  = '/Users/alex 1/Desktop/LAB/LoicLGNinfo_4LD/';
    BRdatafile = '190119_B_cinterocdrft002';
else
    npmkdir    = 'C:\Users\maier\Documents\MATLAB\NPMK-master\'; 
    nbanadir   = 'C:\Users\maier\Documents\bootcamp-selected\nbanalysis\'; 
 
    directory  = 'C:\Users\maier\Documents\LGNinfo_4LD-20190826T172747Z-001\LGNinfo_4LD\';
    BRdatafile = '190119_B_cinterocdrft002';
end


npmkdir    = 'C:\Users\maier\Documents\MATLAB\NPMK-master\'; 
nbanadir   = 'C:\Users\maier\Documents\bootcamp-selected\nbanalysis\'; 

directory  = 'C:\Users\maier\Documents\LGN_data\190119_B_cinterocdrft002_data\';
BRdatafile = '190119_B_cinterocdrft002';
filename   = [directory BRdatafile]; 

addpath(genpath(directory))
addpath(genpath(npmkdir))
addpath(genpath(nbanadir))

aMUA = load(strcat(filename, '_aMUA_hp15.mat'));  
LFP = load(strcat(filename, '_LFP_hp15.mat')); 
CSD = load(strcat(filename, '_CSD_hp15.mat'));

%% STEP ONE: LOAD STIMULUS CONDITIONS with text file 

patterns   = {'rforidrft','rfsfdrft','posdisparitydrft','disparitydrft','cinterocdrft','coneinterocdrft','conedrft', ...
                'colorflicker','bwflicker','rfori','rfsize','cinteroc','color','rfsf','mcosinteroc','dotmapping'}; 

for p = 1:length(patterns)

   pattern      = patterns{p}; 
  
   if any(strfind(BRdatafile,pattern))
       startlog = strfind(BRdatafile,pattern); 
       if ~isequal(BRdatafile(startlog:end-3),pattern),continue
       else
       match    = patterns{p}; 
       end
   end
   
end

if isequal(match,'dotmapping')
ext  = '.gDotsXY_di';
else
ext  = ['.g' upper(match) 'Grating_di']; 
end

if contains(ext,'DRFT')
      grating     = readgDRFTGrating([filename ext]); % from nbanalysis (or even MLAnalysisOnline--might be out of date)
elseif contains(ext,'Dots')
      grating     = readgDotsXY([filename ext]);
else
      grating     = readgGrating([filename ext]);
end


 

%% STEP TWO: LOAD EVENT TIMES/CODES

NEV             = openNEV([filename '.nev'],'noread','overwrite');
EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;        % we don't know why we have to subtract 128 but we do
EventSamples    = NEV.Data.SerialDigitalIO.TimeStamp;                 % in samples 
EventTimes      = floor(NEV.Data.SerialDigitalIO.TimeStampSec.*1000); % convert to ms 
[pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSamples);          % sorts codes, samps or times into trials


% So far all of these data are from EVERY trial, including trials where
% animal breaks fixation. Lets get rid of those and make a new structure
% with the grating info and the stimulus onsets 

STIM            = sortStimandTimeData(grating,pEvC,pEvT,'stim'); % this is in nbanalysis. definitely double check it before you use it. 
%% LOAD LFP with NS2 file

clear ext 

ext = 'ns2';
el  = 'eC';
lfp = getLFP(filename,ext,el,'ascending');        % if you use this function you need to know the bank + sort direction 
                                                  % this function can be found in MLAnalysisOnline or nbanalysis
                                                  % if your sort direction input is correct, the data come out in order
                                                  % from top--> bottom
%% here's a quick way to get  bank info + sort direction if you don't have it in an Excel sheet or notes: 

% Read in NS Header

NS_Header    = openNSx(strcat(filename,'.',ext),'noread');
banks        = unique({NS_Header.ElectrodesInfo.ConnectorBank}); banks(ismember(banks,'E')) = []; % bank E is BNC cable inputs = eye tracking

for b = 1:length(banks)
    clear neural label 
    neural       = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},banks{b}); 
    firstlabel   = cell2mat({NS_Header.ElectrodesInfo(find(neural,1,'first')).Label}); 
    if str2double(firstlabel(3:4)) < 2
        sortdirection = 'ascending'; 
    else
        sortdirection = 'descending'; 
    end
end

%}                                                                                   
%% LOAD LFP and analog MUA with NS6 file ( you can follow these steps to load with ns2 as well but beware of sampling freq.)
%  let's break the whole thing down. 

clear ext NS_header banks neural 

% Read in NS Header
ext          = 'ns6'; 
NS_Header    = openNSx(strcat(filename,'.',ext),'noread');

% get basic info about recorded data
neural       = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},el(2)); % logicals where contact bank name matches electrode of interest
N.neural     = sum(neural); % number of neural channels 
NeuralLabels = {NS_Header.ElectrodesInfo(neural).Label}; %get labels
Fs           = NS_Header.MetaTags.SamplingFreq; % get sampling frequency
nyq          = Fs/2; 
r            = Fs/1000; 
 %{
data_header = struct();

        eval(['data_header.domsup50' ...
            '_nondom0='...
            'STIM.contrast >= 0.5'...
            ' & STIM.fixedc == 0'])
      
       eval(['data_header.dom0' ...
            '_nondomsup50='...
            'STIM.contrast == 0'...
            ' & STIM.fixedc >= 0.5'])
        
       eval(['data_header.domsup50' ...
            '_nondomsup50='...
            'STIM.contrast >= 0.5'...
            ' & STIM.fixedc >= 0.5'])   
  %}
 %WARNING: the following logical vectors are obtained from 'grating' object,
 %itself obtain from the textfile, that has all the event codes of the
 %whole session. However, the '.ns6' file is of different length, as only
 %a certain number of trials of the session are include in -s00x.ns6' files
 % ==> don't forget to adjust the size accordingly when calling the
 % contrast.
data_header.contrast1 = STIM.contrast ==0 & STIM.fixedc >= 0.5;
data_header.contrast2 = STIM.contrast >=0.5 & STIM.fixedc == 0;
data_header.contrast3 = STIM.contrast >=0.5 & STIM.fixedc >= 0.5;

f = {'_dom0_nondomsup50','_domsup50_nondom0','_domsup50_nondomsup50'};

%f.data_header= {'_dom0_nondomsup50','_domsup50_nondom0','_domsup50_nondomsup50'};
%[data_header.(f)] = f{'_dom0_nondomsup50','_domsup50_nondom0','_domsup50_nondomsup50'};

%contrast = '_domsup50_nondomsup50';

%% if you already have the data stored into a file, write it like
%this:
%{
cont_aMUA = aMUA.aMUA_data(1:1551,:,data_header.contrast2(1:length(aMUA.aMUA_data(1,1,:))));
cont_LFP = LFP.LFP_data(1:1551,:,data_header.contrast2(1:length(LFP.LFP_data(1,1,:))));
cont_CSD = CSD.CSD_data(1:1551,:,data_header.contrast2(1:length(CSD.CSD_data(1,1,:))));
%}

cont_aMUA = aMUA.STIM_aMUA(1:1551,:,data_header.contrast3(1:length(aMUA.STIM_aMUA(1,1,:))));
cont_LFP = LFP.STIM_LFP(1:1551,:,data_header.contrast3(1:length(LFP.STIM_LFP(1,1,:))));
cont_CSD = CSD.STIM_CSD(1:1551,:,data_header.contrast3(1:length(CSD.STIM_CSD(1,1,:))));


 %filter aMUA and CSD
        lpc       = 100; %low pass cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low');
        lpaMUA      = filtfilt(bwb,bwa,cont_aMUA);
        lpCSD      = filtfilt(bwb,bwa,cont_CSD);
        

%% Create plots
%% plot baseline corrected mean normalized data

%data = cat(4,cont_LFP, cont_aMUA, cont_CSD);
%newLimePlotExclusive(data, BRdatafile, filename, f{3})
%newLimePlot(data, BRdatafile, filename, f{3})

%compute mean of the data across trials
mean_aMUA= mean(lpaMUA,3);
mean_CSD = mean(lpCSD,3);

 
%std_aMUA = std(lpaMUA,0,3);
%std_CSD = std(lpCSD,0,3);
% if the data is either LFP or MUA, do a normalization (data - min)/(max-min)
 norm_mean_aMUA = nan(length(mean_aMUA(:,1)),length(mean_aMUA(1,:)));
 
 %norm_std_CSD = nan(length(std_CSD(:,1)),length(std_CSD(1,:)));
 %norm_std_aMUA = nan(length(std_aMUA(:,1)),length(std_aMUA(1,:)));
 
for n=1:length(mean_aMUA(1,:))
    
norm_mean_aMUA(:,n) = (mean_aMUA(:,n) - min(mean_aMUA(:,n)))/(max(mean_aMUA(:,n))-min(mean_aMUA(:,n)));
%norm_std_aMUA(:,n) = (std_aMUA(:,n) - min(std_aMUA(:,n)))/(max(std_aMUA(:,n))-min(std_aMUA(:,n)));

end 

%no normalization for CSD


%baseline corrected data process
baseaMUA = norm_mean_aMUA(25:75,:);
baseCSD = mean_CSD(25:75,:);
mean_bp_aMUA = mean(baseaMUA,1);
mean_bp_CSD = mean(baseCSD,1);

std_bp_aMUA = std(baseaMUA,0,1);
std_bp_CSD = std(baseCSD,0,1);

norm_mean_bscorr_aMUA = norm_mean_aMUA(:,:) - mean_bp_aMUA;
%norm_std_bscorr_aMUA = norm_std_aMUA(:,:) - std_bp_aMUA;

    %Make sure CSD gets inverted
norm_mean_bscorr_CSD = -(mean_CSD(:,:) - mean_bp_CSD); 

%scale CSD to [-1 1]
scalingfactor1 = max(abs(norm_mean_bscorr_CSD(:,18)),[],1);
scalingfactor2 = max(abs(norm_mean_bscorr_CSD(:,18)),[],1);
norm_mean_bscorr_CSD = norm_mean_bscorr_CSD .* 1/scalingfactor1; 
norm_mean_bscorr_aMUA =norm_mean_bscorr_aMUA .* 1/scalingfactor2; 
%norm_std_bscorr_CSD = -(norm_std_CSD(:,:) - std_bp_CSD);

%csd_under_lim = norm_mean_bscorr_CSD - norm_std_bscorr_CSD;
%csd_upper_lim = norm_mean_bscorr_CSD + norm_std_bscorr_CSD;

%amua_under_lim = norm_mean_bscorr_aMUA - norm_std_bscorr_aMUA;
%amua_upper_lim = norm_mean_bscorr_aMUA + norm_std_bscorr_aMUA;

%find local maxima to plot a tendency line
[pksaMUA, locsaMUA] = findpeaks(norm_mean_bscorr_aMUA(:,18));
[pksCSD, locsCSD] = findpeaks(norm_mean_bscorr_CSD(:,18) .*(-1));
linreg1 = polyfit(locsaMUA, pksaMUA, 1);
ylinreg1 = linreg1(1) .* locsaMUA + linreg1(2);
linreg2 = polyfit(locsCSD, pksCSD, 1);
ylinreg2 = linreg2(1) .* locsCSD + linreg2(2);

xabs = -50:1500;

h = figure();
%subplot(1,2,1)
%plot(xabs, amua_upper_lim(:, 18))
%hold on
plot(xabs, norm_mean_bscorr_aMUA(:,18))
hold on
%plot(xabs, amua_under_lim(:, 18))
%hold on
%plot(xabs, csd_upper_lim(:, 18))
%hold on
plot(xabs, norm_mean_bscorr_CSD(:, 18))
hold on
%findpeaks(norm_mean_bscorr_aMUA(:,18))

hold on
%findpeaks(norm_mean_bscorr_CSD(:,18))
hold on
plot(locsaMUA, pksaMUA)
hold on
plot(locsCSD, pksCSD*(-1))
hold on
plot(locsaMUA, ylinreg1)
txt1 = ['y = (' num2str(linreg1(1)) ')x + (' num2str(linreg1(2)) ')'];
text(locsaMUA(1), ylinreg1(1), txt1);
hold on
plot(locsCSD, -ylinreg2)
txt2 = ['y = (' num2str(-linreg2(1)) ')x - (' num2str(linreg2(2)) ')'];
text(locsCSD(1), -ylinreg2(1), txt2);
%plot(xabs, csd_under_lim(:, 18))
%hold on
plot([0 0], ylim,'k')
plot([1150 1150], ylim,'k')
plot(xlim ,[0 0],'k')
title({strcat(sprintf('Scaled LGN %s averaged data',  char(BRdatafile))), ...
    strcat('Normalized, bs corrected aMUA over  bs corrected CSD ', f{3})}, ...
    'Interpreter', 'none', 'fontsize', 20, 'FontWeight', 'bold')

xlabel('\fontsize{9}time (ms)')
ylabel({'\fontsize{9}Channel 9','\fontsize{9}(no unit)'})
legend('\fontsize{9}aMUA', '\fontsize{9}CSD')
hold off


%z-score

aMUA_zscore = norm_mean_bscorr_aMUA .* 1/(std_bp_aMUA);
CSD_zscore = norm_mean_bscorr_CSD .* 1/(std_bp_CSD);

[pksaMUAz, locsaMUAz] = findpeaks(aMUA_zscore);
[pksCSDz, locsCSDz] = findpeaks(CSD_zscore .*(-1));

figure();
plot(xabs, aMUA_zscore(:,18))
hold on 
plot(xabs, norm_mean_bscorr_aMUA(:,18))
%x_width=40; y_width=30;
%set(h, 'PaperPosition', [0 0 x_width y_width]);
%saveas(gcf,strcat(filename, 'all',contrast, '_plots_scaletest2', '.png'));
%saveas(gcf,strcat(filename, 'font_size_test_plots', '.png'));


