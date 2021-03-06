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

dom = unique(STIM.contrast);
nondom = unique(STIM.fixedc);
for temp_c = dom
    for temp_f = nondom
disp(temp_c)
disp(temp_f)
end
end
for temp_c = dom'
    for temp_f = nondom'
        

        eval(['data_header.domsup' num2str(temp_c*1000) ...
            'nondom' num2str(temp_f*1000) '='...
            'STIM.contrast == ' num2str(temp_c)...
            ' & STIM.fixedc ==' num2str(temp_f)])
                
    end
end
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

f = {'_dom0_nondomsup50','_domsup50_nondom0','_domsup50_domsup50'};

%f.data_header= {'_dom0_nondomsup50','_domsup50_nondom0','_domsup50_nondomsup50'};
%[data_header.(f)] = f{'_dom0_nondomsup50','_domsup50_nondom0','_domsup50_nondomsup50'};

%contrast = '_domsup50_nondomsup50';

%% if you already have the data stored into a file, write it like
%this:
cont_aMUA = aMUA.aMUA_data(1:1551,:,data_header.contrast1(1:length(aMUA.aMUA_data(1,1,:))));
cont_LFP = LFP.LFP_data(1:1551,:,data_header.contrast1(1:length(LFP.LFP_data(1,1,:))));
cont_CSD = CSD.CSD_data(1:1551,:,data_header.contrast1(1:length(CSD.CSD_data(1,1,:))));

%{
cont_aMUA = aMUA.STIM_aMUA(1:1551,:,data_header.contrast1(1:length(aMUA.STIM_aMUA(1,1,:))));
cont_LFP = LFP.STIM_LFP(1:1551,:,data_header.contrast1(1:length(LFP.STIM_LFP(1,1,:))));
cont_CSD = CSD.STIM_CSD(1:1551,:,data_header.contrast1(1:length(CSD.STIM_CSD(1,1,:))));
%}

%% Create plots
%% plot baseline corrected mean normalized data

data = cat(4,cont_LFP, cont_aMUA, cont_CSD);
%newLimePlotExclusive(data, BRdatafile, filename, f{3})
newLimePlot(data, BRdatafile, filename, f{1})
 
   %LFP
 LimePlot(cont_LFP, BRdatafile, filename, contrast)
   %aMUA
 LimePlot(cont_aMUA, BRdatafile, filename, contrast)
   %CSD
 LimePlot(cont_CSD, BRdatafile, filename, contrast)
 %%2D plots
 SurfacePlot(cont_LFP, BRdatafile, filename, contrast)
 SurfacePlot(cont_aMUA, BRdatafile, filename, contrast)
 SurfacePlot(cont_CSD, BRdatafile, filename, contrast)

 %% stimuli duration time
onsettimes =EventTimes(EventCodes ==24);
offsettimes = EventTimes(EventCodes ==23);
correctev = EventTimes(EventCodes ==96);

selected_onsettimes= nan(1,length(onsettimes));
selected_offsettimes=nan(1,length(offsettimes));
for n = 1:length(correctev) 
  selected_onsettimes(n)  = EventTimes(find(EventCodes(1:length(correctev)) == 23,1,'last'));
  selected_offsettimes(n) = EventTimes(find(EventCodes(1:length(correctev)) == 24,1,'last'));
end
    
     
        
stimtime1 = selected_onsettimes(1)- selected_offsettimes(1);
stimtime2 = selected_onsettimes(2)- selected_offsettimes(2);
stimtime3 = selected_onsettimes(3)- selected_offsettimes(3);

%%
 %{
 mean_LFP = mean(cont_LFP,3);
 norm_mean = nan(length(mean_data(:,1)),length(mean_data(1,:)));
     for n=1:24
     norm_mean(:,n) = (mean_LFP(:,n) - min(mean_LFP(:,n)))/(max(mean_LFP(:,n))-min(mean_LFP(:,n)));
     end 
 baseline = mean(norm_mean(1:50,:),1);
 bscorr = norm_mean -baseline;
 figure();
 f_ShadedLinePlotbyDepth(bscorr,1:length(bscorr(1,:)),-50:1500,1:length(bscorr(1,:)),1)
}%
%% LOAD discrete MUA with ppnev file
% to get the ppNEV files [post-processed NEV] use the offlineBRAutoSort
% directory under https://github.com/maierav/KiloSortUtils
% input is threshold (std) to apply to envelope to extract spikes

% let me show you where it is...

clear Fs 

ppNEV = load([filename '.ppNEV'],'-MAT','ppnev');

Fs    = double(ppNEV.ppnev.MetaTags.SampleRes);

for i = 1:length(sortedLabels)
    clear elabel 
    elabel               = sortedLabels{i}; 
    eidx                 = find(cell2mat(cellfun(@(x) contains(x',elabel),{ppNEV.ppnev.ElectrodesInfo.ElectrodeLabel},'UniformOutput',0)));
    I                    = ppNEV.ppnev.Data.Spikes.Electrode == eidx;
    SPK                  = double(ppNEV.ppnev.Data.Spikes.TimeStamp(I)); % these are the spike times (in samples)
    sdf                  = spk2sdf(SPK,Fs); % this convolves the spikes. you need jnm_kernel to use the poisson dist 
   
    STIM.dMUA(:,i,:)     = squeeze(trigData(sdf',floor(STIM.onsets./30),-pre,post));  % convolved ppNEV cmua sorted into trials 
    STIM.dMUA_bin(:,i,:) = squeeze(trigBinaryData(SPK,pre,post,floor(STIM.onsets./30))); % binary ppNEV cmua sorted into trials 
end


%% LOOK at kilosorted data 
% use https://github.com/maierav/KiloSortUtils
% use phy to view results! 
% let me show you what it looks like on the MacPro..


%% STEP FOUR: plot the data ! the best part 
mean_dMUA= mean(STIM.dMUA, 3);

%dMUA
meanlimeplotLoic(mean_dMUA, BRdatafile)
SurfacePlot(mean_dMUA, BRdatafile)
%% Really useful Matlab central functions for plotting data to make everything easier and look better: 


% EXPORT FIG -- perfect for making pdfs when you are checking data!


% COLORS
% https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
% see also: http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3


% BAR PLOTS
% https://www.mathworks.com/matlabcentral/fileexchange/57499-superbarsuperbar 




% Read in NS Header


%%
