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

directory  = 'C:\Users\maier\Documents\LGNinfo_4LD-20190826T172747Z-001\LGNinfo_4LD\';
BRdatafile = '190119_B_cinterocdrft002';
filename   = [directory BRdatafile]; 


addpath(genpath(directory))
addpath(genpath(npmkdir))
addpath(genpath(nbanadir))

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




%% STEP THREE: LOAD NEURAL DATA

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

% counters
clear nct
nct = 0;

tic 
% process data electrode by electrode
for e = 1:length(neural)
    
    if neural(e) == 1    % why? because neural is a vector of logicals, 1 = contacts we want

        nct = nct+1;
        
        % open data for this channel. 
        clear NS DAT
        electrode = sprintf('c:%u',e);
        NS        = openNSx(strcat(filename,'.',ext),electrode,'read','uV');
        DAT       = NS.Data; NS.Data = [];  % this is the whole signal on one channel, 30 kHz!
        
        
        % preallocate data matrices 
        if nct == 1
            N.samples = length(DAT); 
            LFP       = zeros(ceil(N.samples/r),N.neural); % preallocating for downsampled data
            MUA       = zeros(ceil(N.samples/r),N.neural);
        end
        
        % extract the LFP. 
        clear lpc lWn bwb bwa lpLFP
        lpc       = 200; %low pass cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low');
        lpLFP      = filtfilt(bwb,bwa,DAT);  %low pass filter 
        
        % extract the MUA:
        clear hpc hWn bwb bwa hpMUA
        hpc       = 750;  %high pass cutoff
        hWn       = hpc/nyq;
        [bwb,bwa] = butter(4,hWn,'high');
        hpMUA     = filtfilt(bwb,bwa,DAT); %high pass filter
        
        % low pass at 5000 Hz and rectify 
        clear lpc lWn bwb bwa 
        lpc       = 5000;  % cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low');
        hpMUA     = abs(filtfilt(bwb,bwa,hpMUA)); %low pass filter &rectify
        
        % low pass filter at x Hz. 
        clear lpc lWn bwb bwa lpMUA
        lpc       = 200; %low pass cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low'); 
        lpMUA     = filtfilt(bwb,bwa,hpMUA);  %low pass filter to smooth
        
      
        % decimate both LFP and analog MUA (aMUA) to get 1kHz samp freq
        MUA(:,nct) = decimate(lpMUA,r); 
        LFP(:,nct) = decimate(lpLFP,r); 
        
        clear DAT 
        
    end
    
end
toc 


% Warning! THESE DATA ARE IN THE SAME ORDER AS THE BR PINS, NOT THE ORDER OF THE PROBE

% sort data from top of electrode to bottom.
idx = zeros(1,length(NeuralLabels));
for i = 1:length(NeuralLabels)
    
   Str  = cell2mat(NeuralLabels(i));
   Key   = 'eC';
   Str(strfind(Str, '%02d')) = [];
   
   Index = strfind(Str, Key);
   idx(1, i) = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
  
end
% let's type this out together 

Srt_MUA = nan(length(MUA(:,1)), length(NeuralLabels));
Srt_LFP = nan(length(MUA(:,1)), length(NeuralLabels));


Srt_MUA(:,idx) = MUA(:,:); 
Srt_LFP(:,idx) = LFP(:,:);

%% calculate CSD 
% calculate CSD before triggering to trials OR on the trial data BUT not on
% the mean LFP. 

CSD = mod_iCSD(Srt_LFP')';  % this function takes LFP in channels x samples so let's transpose LFP and then flip it right back 
                        % feed in units of microV and get back units of
                        % nA/mm^3
% pad array if you want to keep the matrix the same size on the channel
% dimension as the other matrices

CSD = padarray(CSD,[0 1],NaN,'replicate');

%% trigger the neural data to the event codes of interest
pre   = -50;
post  = 300; 

STIM.LFP  = trigData(Srt_LFP,floor(STIM.onsets./30),-pre,post); % this function is MLAnalysisOnline or nbanalysis. pre variable is in absolute units 
STIM.CSD  = trigData(CSD,floor(STIM.onsets./30),-pre,post); 
STIM.aMUA = trigData(Srt_MUA,floor(STIM.onsets./30),-pre,post); 

%% Create plots
mean_LFP = mean(STIM.LFP, 3);
mean_CSD = mean(STIM.CSD, 3);
mean_aMUA= mean(STIM.aMUA, 3);


sortedLabels = 1:length(idx);
xabs = -50:300;
 
h = figure();
[ha, pos] = tight_subplot(24,1,[0.001 .03],[.05 .05],[.1 .01]);
for i = 1:length(idx)
    
    axes(ha(i));
   
    plot(xabs,mean_LFP(1:351,i))
    Aire = area(xabs, mean_LFP(1:351,i));
    Aire.FaceColor =[0.85 0.85 0.85];
    set(h,'position',get(h,'position').*[1 1 1 1.25]);
    xlim([-50 300]);
    ax = gca;
    ax.YAxisLocation = 'origin';
    grid off
    
    if i < 24
    set(ha(i), 'XTick', [])
    
    ax.XAxis.Visible = 'off';
    end
    if i <= 24
        set(ha(i), 'YTick', [])
    end
    if i == 1
    title({'\fontsize{18}LGN 190119_B_cinterocdrft002 LFP'})
    end
     ylh = ylabel(sortedLabels(i));
     %ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1)*5);
   if i == 12
       tit = title(gca, 'Voltage (µV)','interpreter','none','fontsize',17, 'FontWeight', 'Normal');
       tit.Position(1) = tit.Position(1) - abs(tit.Position(1))*1.5;
       %tit.Position(2) = tit.Position(2) + abs(tit.Position(2))*200;
       tit.Rotation =90;
   end
    
end
 set(ha(1:23), 'XTick', [])
 
 %set(ha(1:23), 'XTickLabel',''); 
 set(ha(1:24), 'box', 'off');
 %set(ha, 'YTickLabel', '');
 xlh = xlabel('\fontsize{18}time (ms)');
 xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2)*0.5);
 
 %aMUA
h = figure();
[ha, pos] = tight_subplot(24,1,[0.001 .03],[.05 .05],[.1 .01]);
for i = 1:length(idx)
    
    axes(ha(i));
   
    plot(xabs,mean_aMUA(1:351,i))
    Aire = area(xabs, mean_aMUA(1:351,i));
    Aire.FaceColor =[0.85 0.85 0.85];
    set(h,'position',get(h,'position').*[1 1 1 1.25]);
    xlim([-50 300]);
    ax = gca;
    ax.YAxisLocation = 'origin';
    grid off
    
    if i < 24
    set(ha(i), 'XTick', [])
    
    ax.XAxis.Visible = 'off';
    end
    if i <= 24
        set(ha(i), 'YTick', [])
    end
    if i == 1
    title({'\fontsize{18}LGN 190119_B_cinterocdrft002 aMUA'})
    end
     ylh = ylabel(sortedLabels(i));
     %ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1)*5);
   if i == 12
       tit = title(gca, 'Voltage (µV)','interpreter','none','fontsize',17, 'FontWeight', 'Normal');
       tit.Position(1) = tit.Position(1) - abs(tit.Position(1))*1.5;
       %tit.Position(2) = tit.Position(2) + abs(tit.Position(2))*200;
       tit.Rotation =90;
   end
    
end
 set(ha(1:23), 'XTick', [])
 
 %set(ha(1:23), 'XTickLabel',''); 
 set(ha(1:24), 'box', 'off');
 %set(ha, 'YTickLabel', '');
 xlh = xlabel('\fontsize{18}time (ms)');
 xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2)*0.5);
 
 %CSD
 
 h = figure();
[ha, pos] = tight_subplot(24,1,[0.001 .03],[.05 .05],[.1 .01]);
for i = 1:length(idx)
    
    axes(ha(i));
   
    plot(xabs,mean_CSD(1:351,i))
    Aire = area(xabs, mean_CSD(1:351,i));
    Aire.FaceColor =[0.85 0.85 0.85];
    set(h,'position',get(h,'position').*[1 1 1 1.25]);
    xlim([-50 300]);
    ax = gca;
    ax.YAxisLocation = 'origin';
    grid off
    
    if i < 24
    set(ha(i), 'XTick', [])
    
    ax.XAxis.Visible = 'off';
    end
    if i <= 24
        set(ha(i), 'YTick', [])
    end
    if i == 1
    title({'\fontsize{18}LGN 190119_B_cinterocdrft002 CSD'})
    end
     ylh = ylabel(sortedLabels(i));
     %ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1)*5);
   if i == 12
       tit = title(gca, 'Voltage (µV)','interpreter','none','fontsize',17, 'FontWeight', 'Normal');
       tit.Position(1) = tit.Position(1) - abs(tit.Position(1))*1.5;
       %tit.Position(2) = tit.Position(2) + abs(tit.Position(2))*200;
       tit.Rotation =90;
   end
    
end
 set(ha(1:23), 'XTick', [])
 
 %set(ha(1:23), 'XTickLabel',''); 
 set(ha(1:24), 'box', 'off');
 %set(ha, 'YTickLabel', '');
 xlh = xlabel('\fontsize{18}time (ms)');
 xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2)*0.5);
 
 
%% plot normalized data
norm_LFP = nan(24,351);
norm_aMUA = nan(24,351);
norm_CSD = nan(24,351);
for n=1:24
norm_LFP(n,:) = (mean_LFP(:,n) - min(mean_LFP(:,n)))/(max(mean_LFP(:,n))-min(mean_LFP(:,n)));
norm_aMUA(n,:) = (mean_aMUA(:,n) - min(mean_aMUA(:,n)))/(max(mean_aMUA(:,n))-min(mean_aMUA(:,n)));
norm_CSD(n,:) = (mean_CSD(:,n) - min(mean_CSD(:,n)))/(max(mean_CSD(:,n))-min(mean_CSD(:,n)));
end       
sortedLabels = 1:length(idx);
xabs = -50:300;
 
   %LFP
h = figure();
[ha, pos] = tight_subplot(24,1,[0.001 .03],[.05 .05],[.1 .01]);
for i = 1:length(idx)
    
    axes(ha(i));
   
    plot(xabs,norm_LFP(i,1:351))
    Aire = area(xabs, norm_LFP(i,1:351));
    Aire.FaceColor =[0.85 0.85 0.85];
    set(h,'position',get(h,'position').*[1 1 1 1.25]);
    %xlim([-50 300]);
    ax = gca;
    ax.YAxisLocation = 'origin';
    grid off
    
    if i < length(mean_LFP(1,:))
    set(ha(i), 'XTick', [])
    ax.XAxis.Visible = 'off';
    end
    if i <= length(mean_LFP(1,:))
        set(ha(i), 'YTick', [])
    end
    if i == 1
    title({'\fontsize{18}LGN 190119_B_cinterocdrft002 normalized LFP'})
    end
     ylh = ylabel(sortedLabels(i));
     ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1)*5);
   if i == length(mean_LFP(1,:))/2
       tit = title(gca, 'Voltage (µV)','interpreter','none','fontsize',17, 'FontWeight', 'Normal');
       tit.Position(1) = tit.Position(1) - abs(tit.Position(1))*1.5;
       tit.Rotation =90;
   end
    
end
 set(ha(1:23), 'XTick', [])
 %set(ha(1:23), 'XTickLabel',''); 
 set(ha(1:24), 'box', 'off');
 %set(ha, 'YTickLabel', '');
 xlh = xlabel('\fontsize{18}time (ms)');
 xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2)*0.5);

 %aMUA  
 h = figure();
[ha, pos] = tight_subplot(24,1,[0.001 .03],[.05 .05],[.1 .01]);
for i = 1:length(idx)
    
    axes(ha(i));
   
    plot(xabs,norm_aMUA(i,1:351))
    Aire = area(xabs, norm_aMUA(i,1:351));
    Aire.FaceColor =[0.85 0.85 0.85];
    set(h,'position',get(h,'position').*[1 1 1 1.25]);
    xlim([-50 300]);
    ax = gca;
    ax.YAxisLocation = 'origin';
    grid off
    
    if i < 24
    set(ha(i), 'XTick', [])
    
    ax.XAxis.Visible = 'off';
    end
    if i <= 24
        set(ha(i), 'YTick', [])
    end
    if i == 1
    title({'\fontsize{18}LGN 190119_B_cinterocdrft002 normalized aMUA'})
    end
     ylh = ylabel(sortedLabels(i));
     %ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1)*5);
   if i == 12
       tit = title(gca, 'Voltage (µV)','interpreter','none','fontsize',17, 'FontWeight', 'Normal');
       tit.Position(1) = tit.Position(1) - abs(tit.Position(1))*1.5;
       %tit.Position(2) = tit.Position(2) + abs(tit.Position(2))*200;
       tit.Rotation =90;
   end
    
end
 set(ha(1:23), 'XTick', [])
 
 %set(ha(1:23), 'XTickLabel',''); 
 set(ha(1:24), 'box', 'off');
 %set(ha, 'YTickLabel', '');
 xlh = xlabel('\fontsize{18}time (ms)');
 xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2)*0.5);

 %CSD
 h = figure();
[ha, pos] = tight_subplot(24,1,[0.001 .03],[.05 .05],[.1 .01]);
for i = 1:length(idx)
    
    axes(ha(i));
   
    plot(xabs,norm_CSD(i,1:351))
    Aire = area(xabs, norm_CSD(i,1:351));
    Aire.FaceColor =[0.85 0.85 0.85];
    set(h,'position',get(h,'position').*[1 1 1 1.25]);
    xlim([-50 300]);
    ax = gca;
    ax.YAxisLocation = 'origin';
    grid off
    
    if i < 24
    set(ha(i), 'XTick', [])
    
    ax.XAxis.Visible = 'off';
    end
    if i <= 24
        set(ha(i), 'YTick', [])
    end
    if i == 1
    title({'\fontsize{18}LGN 190119_B_cinterocdrft002 normalized CSD'})
    end
     ylh = ylabel(sortedLabels(i));
     %ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1)*5);
   if i == length(mean_CSD(1,:))/2
       
       tit = title(gca, 'Voltage (µV)','interpreter','none','fontsize',17, 'FontWeight', 'Normal');
       tit.Position(1) = tit.Position(1) - abs(tit.Position(1))*1.5;
       %tit.Position(2) = tit.Position(2) + abs(tit.Position(2))*200;
       tit.Rotation =90;
   end
    
end
 set(ha(1:23), 'XTick', [])
 
 %set(ha(1:23), 'XTickLabel',''); 
 set(ha(1:24), 'box', 'off');
 %set(ha, 'YTickLabel', '');
 xlh = xlabel('\fontsize{18}time (ms)');
 xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2)*0.5);

 %%2D plots
 h = figure();
 xabs = [-50:300];
 yvec = [1:24];
 imagesc(xabs,yvec, norm_LFP)
 title({'\fontsize{18}LGN 190119_B_cinterocdrft002 normalized LFP'})

 h = figure();
 xabs = [-50:300];
 yvec = [1:24];
 imagesc(xabs,yvec, norm_aMUA)
 title({'\fontsize{18}LGN 190119_B_cinterocdrft002 normalized aMUA'})
 
 h = figure();
 xabs = [-50:300];
 yvec = [1:24];
 imagesc(xabs,yvec, norm_CSD)
 title({'\fontsize{18}LGN 190119_B_cinterocdrft002 normalized CSD'})
 
%% LOAD discrete MUA with ppnev file
% to get the ppNEV files [post-processed NEV] use the offlineBRAutoSort
% directory under https://github.com/maierav/KiloSortUtils
% input is threshold (std) to apply to envelope to extract spikes

% let me show you where it is...

clear Fs 

load([filename '.ppNEV'],'-MAT','ppNEV')

Fs    = double(ppNEV.MetaTags.SampleRes);

for i = 1:length(sortedLabels)
    clear elabel 
    elabel               = sortedLabels{i}; 
    eidx                 = find(cell2mat(cellfun(@(x) contains(x',elabel),{ppNEV.ElectrodesInfo.ElectrodeLabel},'UniformOutput',0)));
    I                    = ppNEV.Data.Spikes.Electrode == eidx;
    SPK                  = double(ppNEV.Data.Spikes.TimeStamp(I)); % these are the spike times (in samples)
    sdf                  = spk2sdf(SPK,Fs); % this convolves the spikes. you need jnm_kernel to use the poisson dist 
   
    STIM.dMUA(:,i,:)     = squeeze(trigData(sdf',floor(STIM.onsets./30),-pre,post));  % convolved ppNEV cmua sorted into trials 
    STIM.dMUA_bin(:,i,:) = squeeze(trigBinaryData(SPK,pre,post,floor(STIM.onsets./30))); % binary ppNEV cmua sorted into trials 
end


%% LOOK at kilosorted data 
% use https://github.com/maierav/KiloSortUtils
% use phy to view results! 
% let me show you what it looks like on the MacPro..


%% STEP FOUR: plot the data ! the best part 
fig1 = figure();
subplot(11,1,1);
plot(Srt_LFP(1:400,1))
subplot(11,1,2);
plot(Srt_LFP(1:400,2))
subplot(11,1,3);
plot(Srt_LFP(1:400,3))
subplot(11,1,4);
plot(Srt_LFP(1:400,4))
subplot(11,1,5);
plot(Srt_LFP(1:400,5))
subplot(11,1,6);
plot(Srt_LFP(1:400,6))
subplot(11,1,7);
plot(Srt_LFP(1:400,7))
subplot(11,1,8);
plot(Srt_LFP(1:400,8))
subplot(11,1,9);
plot(Srt_LFP(1:400,9))
subplot(11,1,10);
plot(Srt_LFP(1:400,10))
subplot(11,1,11);
plot(Srt_LFP(1:400,11))

fig2 = figure();
subplot(13,1,1);
plot(Srt_LFP(1:400,12))
subplot(13,1,2);
plot(Srt_LFP(1:400,13))
subplot(13,1,3);
plot(Srt_LFP(1:400,14))
subplot(13,1,4);
plot(Srt_LFP(1:400,15))
subplot(13,1,5);
plot(Srt_LFP(1:400,16))
subplot(13,1,6);
plot(Srt_LFP(1:400,17))
subplot(13,1,7);
plot(Srt_LFP(1:400,18))
subplot(13,1,8);
plot(Srt_LFP(1:400,19))
subplot(13,1,9);
plot(Srt_LFP(1:400,20))
subplot(13,1,10);
plot(Srt_LFP(1:400,21))
subplot(13,1,11);
plot(Srt_LFP(1:400,22))
subplot(13,1,12);
plot(Srt_LFP(1:400,23))
subplot(13,1,13);
plot(Srt_LFP(1:400,24))

%% Really useful Matlab central functions for plotting data to make everything easier and look better: 


% EXPORT FIG -- perfect for making pdfs when you are checking data!


% COLORS
% https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
% see also: http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3


% BAR PLOTS
% https://www.mathworks.com/matlabcentral/fileexchange/57499-superbarsuperbar 




% Read in NS Header


%%
%{
drname       = '/Users/loicdaumail/Documents/Vanderbilt/Maier_lab/Kacie_data/LGNinfo_4LD/'; % ADJUST
filename     = [drname '190119_B_cinterocdrft002'];
extension    = 'ns6';
NS_Header    = openNSx(strcat(filename,'.',extension),'noread');

banks        = unique({NS_Header.ElectrodesInfo.ConnectorBank}); banks(ismember(banks,'E')) = [];

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
%neural       = ~strcmp('E',{NS_Header.ElectrodesInfo.ConnectorBank}); % banks A-D are neural; bank E is the BNCs on the front of the NSP


N.elect      = length(neural);
N.neural     = sum( neural);

%get labels
NeuralLabels = {NS_Header.ElectrodesInfo(neural).Label};
NeuralInfo   = NS_Header.ElectrodesInfo(neural);

% get sampling frequency
Fs           = NS_Header.MetaTags.SamplingFreq;
nyq          = Fs/2;
r            = Fs/1000; % 1000 is the sampling frequency we want after decimation

% counters
clear nct
nct = 0;

% process data electrode by electrode
for e = 1:length(neural)
    
    if neural(e) == 1

        nct = nct+1;
        
        %open data for this chanel
        clear NS DAT
        electrode = sprintf('c:%u',e);
        NS        = openNSx(strcat(filename,'.',extension),electrode,'read','uV');
        if iscell(NS.Data)
            DAT   = cell2mat(NS.Data);
        else
            DAT   = NS.Data;
        end
        NS.Data   = [];
        
        
        % deminate DAT to turn the 30 kHz signal into a 1 kHz one (you can
        % use the r variable to do this) 
        
        
        
            %preallocation -- PREALLOCATE. name the variable you're going
            %to collect the data in
         % preallocate data matrices 
        if nct == 1
            N.samples = length(DAT); 
            LFP       = zeros(ceil(N.samples/r),N.neural); % preallocating for downsampled data
            MUA       = zeros(ceil(N.samples/r),N.neural);
        end
        

        % extract the LFP. 
        clear lpc lWn bwb bwa lpLFP
        lpc       = 200; %low pass cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low');
        lpLFP      = filtfilt(bwb,bwa,DAT);  %low pass filter 
        
        % extract the MUA:
        clear hpc hWn bwb bwa hpMUA
        hpc       = 750;  %high pass cutoff
        hWn       = hpc/nyq;
        [bwb,bwa] = butter(4,hWn,'high');
        hpMUA     = filtfilt(bwb,bwa,DAT); %high pass filter
        
        % low pass at 5000 Hz and rectify 
        clear lpc lWn bwb bwa 
        lpc       = 5000;  % cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low');
        hpMUA     = abs(filtfilt(bwb,bwa,hpMUA)); %low pass filter &rectify
        
        % low pass filter at x Hz. 
        clear lpc lWn bwb bwa lpMUA
        lpc       = 200; %low pass cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low'); 
        lpMUA     = filtfilt(bwb,bwa,hpMUA);  %low pass filter to smooth
        
             
       % collect the data in that varibale
        % decimate both LFP and analog MUA (aMUA) to get 1kHz samp freq
        MUA(:,nct) = decimate(lpMUA,r); 
        LFP(:,nct) = decimate(lpLFP,r); 
         
       % then clear that channels data to keep memory load lower 
       clear DAT
        
    end
    
end

        
% sort electrode contacts in ascending order:
%{
idx         = nan(length(NeuralLabels),1);

for ch      = 1:length(NeuralLabels)
    %chname  = NeuralLabels{ch};
    chname  = strcat(sprintf('%s','eC'),sprintf('%02d',ch)); 
    content = contains(NeuralLabels, chname);
    idx(ch) = sprintf(NeuralLabels(find(content)));
   
end
%}


idx = zeros(1,length(NeuralLabels));
for i = 1:length(NeuralLabels)
    
   Str  = cell2mat(NeuralLabels(i));
   Key   = 'eC';
   Str(strfind(Str, '%02d')) = [];
   
   Index = strfind(Str, Key);
   idx(1, i) = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
  
end
    
% NOW USE IDX (labelnum) TO RE-SORT YOUR DATA VARIABLE

Srt_MUA = nan(length(MUA(:,1)), length(NeuralLabels));
Srt_LFP = nan(length(MUA(:,1)), length(NeuralLabels));


Srt_MUA(:,idx) = MUA(:,:); 
Srt_LFP(:,idx) = LFP(:,:);

fig1 = figure();
subplot(11,1,1);
plot(Srt_LFP(1:400,1))
subplot(11,1,2);
plot(Srt_LFP(1:400,2))
subplot(11,1,3);
plot(Srt_LFP(1:400,3))
subplot(11,1,4);
plot(Srt_LFP(1:400,4))
subplot(11,1,5);
plot(Srt_LFP(1:400,5))
subplot(11,1,6);
plot(Srt_LFP(1:400,6))
subplot(11,1,7);
plot(Srt_LFP(1:400,7))
subplot(11,1,8);
plot(Srt_LFP(1:400,8))
subplot(11,1,9);
plot(Srt_LFP(1:400,9))
subplot(11,1,10);
plot(Srt_LFP(1:400,10))
subplot(11,1,11);
plot(Srt_LFP(1:400,11))

fig2 = figure();
subplot(13,1,1);
plot(Srt_LFP(1:400,12))
subplot(13,1,2);
plot(Srt_LFP(1:400,13))
subplot(13,1,3);
plot(Srt_LFP(1:400,14))
subplot(13,1,4);
plot(Srt_LFP(1:400,15))
subplot(13,1,5);
plot(Srt_LFP(1:400,16))
subplot(13,1,6);
plot(Srt_LFP(1:400,17))
subplot(13,1,7);
plot(Srt_LFP(1:400,18))
subplot(13,1,8);
plot(Srt_LFP(1:400,19))
subplot(13,1,9);
plot(Srt_LFP(1:400,20))
subplot(13,1,10);
plot(Srt_LFP(1:400,21))
subplot(13,1,11);
plot(Srt_LFP(1:400,22))
subplot(13,1,12);
plot(Srt_LFP(1:400,23))
subplot(13,1,13);
plot(Srt_LFP(1:400,24))

%% calculate CSD 
% calculate CSD before triggering to trials OR on the trial data BUT not on
% the mean LFP. 

CSD = mod_iCSD(Srt_LFP')';  % this function takes LFP in channels x samples so let's transpose LFP and then flip it right back 
                        % feed in units of microV and get back units of
                        % nA/mm^3
% pad array if you want to keep the matrix the same size on the channel
% dimension as the other matrices

CSD = padarray(CSD,[0 1],NaN,'replicate');

%% trigger the neural data to the event codes of interest
pre   = -50;
post  = 300; 

STIM.LFP  = trigData(LFP,floor(STIM.onsets./30),-pre,post); % this function is MLAnalysisOnline or nbanalysis. pre variable is in absolute units 
%STIM.CSD  = trigData(CSD,floor(STIM.onsets./30),-pre,post); 
STIM.aMUA = trigData(MUA,floor(STIM.onsets./30),-pre,post); 

%% LOAD discrete MUA with ppnev file
% to get the ppNEV files [post-processed NEV] use the offlineBRAutoSort
% directory under https://github.com/maierav/KiloSortUtils
% input is threshold (std) to apply to envelope to extract spikes

% let me show you where it is...

clear Fs 

load([filename '.ppNEV'],'-MAT','ppNEV')

Fs    = double(ppNEV.MetaTags.SampleRes);

for i = 1:length(sortedLabels)
    clear elabel 
    elabel               = sortedLabels{i}; 
    eidx                 = find(cell2mat(cellfun(@(x) contains(x',elabel),{ppNEV.ElectrodesInfo.ElectrodeLabel},'UniformOutput',0)));
    I                    = ppNEV.Data.Spikes.Electrode == eidx;
    SPK                  = double(ppNEV.Data.Spikes.TimeStamp(I)); % these are the spike times (in samples)
    sdf                  = spk2sdf(SPK,Fs); % this convolves the spikes. you need jnm_kernel to use the poisson dist 
   
    STIM.dMUA(:,i,:)     = squeeze(trigData(sdf',floor(STIM.onsets./30),-pre,post));  % convolved ppNEV cmua sorted into trials 
    STIM.dMUA_bin(:,i,:) = squeeze(trigBinaryData(SPK,pre,post,floor(STIM.onsets./30))); % binary ppNEV cmua sorted into trials 
end


% THEN FIND A WAY TO PLOT ~5 seconds worth of data for each of the
% channels. 
%}
