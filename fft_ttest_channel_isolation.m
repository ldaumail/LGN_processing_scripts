
%load data (loop through the folder)

npmkdir    = 'C:\Users\maier\Documents\MATLAB\NPMK-master\'; 
nbanadir   = 'C:\Users\maier\Documents\bootcamp-selected\nbanalysis\'; 

directory  = 'C:\Users\maier\Documents\LGN_data\kacie_preproc\';
addpath(genpath(directory))
addpath(genpath(npmkdir))
addpath(genpath(nbanadir))

files = dir(strcat(directory, '*cinterocdrft*'));
for file = 1:length(files)
    
STIMBRdatafile = files(file).name;
STIMfilename   = [directory STIMBRdatafile]; 
STIMdMUA = load(strcat(STIMfilename));  

power = nan(length(STIMdMUA.STIM.sdftr(1,:,1)),length(STIMdMUA.STIM.sdftr(1,1,:)),1025);
freq = nan(length(STIMdMUA.STIM.sdftr(1,:,1)),length(STIMdMUA.STIM.sdftr(1,1,:)),1025);
fourhzpower =nan(length(STIMdMUA.STIM.sdftr(1,:,1)),length(STIMdMUA.STIM.sdftr(1,1,:)));
for chan = 1:length(STIMdMUA.STIM.sdftr(1,:,1))
    for tr = 1:length(STIMdMUA.STIM.sdftr(1,1,:))
%compute fft for every trial of a channel, for every channel   
[power(chan,tr,:), freq(chan,tr,:)] = calcFFT(squeeze(STIMdMUA.STIM.sdftr(600:1700,chan,tr)));
%find the index of the frequency vector closest to 4hz and point to the
%power value of this index for every trial
[val,index] = min(abs(4-freq(chan,tr,:)));
fourhzpower(chan,tr) = power(chan,tr,index);
    end
end
%perform t-test accross :
         %- DE50-NDE0 vs DE0-NDE0
         %- NDE50-DE0 vs DE0-NDE0
         %- NDE50-DE50 vs DE0-NDE0

end




% only save all the trials of channels fo which it was significant

         
 
