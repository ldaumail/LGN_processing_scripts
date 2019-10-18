
npmkdir    = 'C:\Users\maier\Documents\MATLAB\NPMK-master\'; 
nbanadir   = 'C:\Users\maier\Documents\bootcamp-selected\nbanalysis\'; 

directory  = 'C:\Users\maier\Documents\LGN_data\190119_B_cinterocdrft002_data\';
BRdatafile = '190119_B_cinterocdrft002';
filename   = [directory BRdatafile]; 
%{
% Set up variables --------------------------------------------------------
folder = fileparts(mfilename('fullpath'));
if ( isempty( folder) )
    folder = pwd;
end
folder = [folder '\data'];
%}
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


eye_info = concatBHV(strcat(filename, '.bhv'));

 session = 'test';
 samplerate = 500;

 trialindex = STIM.trial;
 
 codes                 = eye_info.CodeNumbers{trialindex};
 times                 = eye_info.CodeTimes{trialindex};
 samples = [];
 samples(:,1) =  (-1*times(codes == 23)+1) : 1 : 0 : (length(eye_info.AnalogData{1,1}.EyeSignal(:,1)) - times(codes == 23)); 
 %1:length(eye_info.AnalogData{1,1}.EyeSignal(:,1)); %timestamps of the recording in miliseconds
 samples(:,2) = eye_info.AnalogData{1,1}.EyeSignal(:,1); %horizontal position of the left eye in degrees
 samples(:,3) = eye_info.AnalogData{1,1}.EyeSignal(:,2); %vertical position of the left eye in degrees 
 samples(:,4) = nan();
 samples(:,5) = nan();
 recording = ClusterDetection.EyeMovRecording.Create(directory, session, samples, [], samplerate);

% Runs the saccade detection
[saccades stats] = recording.FindSaccades();

% Plots a main sequence
enum = ClusterDetection.SaccadeDetector.GetEnum;
figure
subplot(2,2,1)
plot(saccades(:,enum.amplitude),saccades(:,enum.peakVelocity),'o')
set(gca,'xlim',[0 1],'ylim',[0 100]);
xlabel('Saccade amplitude (deg)');
ylabel('Saccade peak velocity (deg/s)');


% Plots the traces with the labeled microsaccades
subplot(2,2,[3:4])
plot(samples(:,1), samples(:,2:end));
hold
yl = get(gca,'ylim');
u1= zeros(size(samples(:,1)))+yl(1);
u2= zeros(size(samples(:,1)))+yl(1);
u1((saccades(:,enum.startIndex))) = yl(2);
u2(saccades(:,enum.endIndex)) = yl(2);
u = cumsum(u1)-cumsum(u2);
plot(samples(:,1), u,'k')

xlabel('Time (ms)');
ylabel('Eye Position (deg)');

legend({'Left Horiz', 'Left Vert', 'Right Horiz' , 'Right Vert', 'Microsaccades'})

 
 
 