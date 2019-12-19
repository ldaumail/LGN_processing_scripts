npmkdir    = 'C:\Users\maier\Documents\MATLAB\NPMK-master\'; 
nbanadir   = 'C:\Users\maier\Documents\bootcamp-selected\nbanalysis\'; 

directory  = 'C:\Users\maier\Documents\LGN_data\rfori\';
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


for i = 1:length(data.datafile)
%experiment conditions
data_header.contrast1 = data.datafile(i).STIM.contrast >= 0.5 & data.datafile(i).STIM.fixedc == 0;
f = {'_DE50_NDE0'};
xabs = -100:200;
idx = [1 3 5 7 9 11 13 2 4 6 8 10 12 14];
nyq = 15000;
mean_data = nanmean(squeeze(data.datafile(i).STIM.sdftr(500:800,:)),2);
 sp = subplot(length(1:7), 2, idx(i));
    plot(xabs, mean_data)
    hold on
    plot([0 0], ylim,'k')
    hold on
    %plot([1150 1150], ylim,'k')
end