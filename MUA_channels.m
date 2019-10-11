npmkdir    = 'C:\Users\maier\Documents\MATLAB\NPMK-master\'; 
nbanadir   = 'C:\Users\maier\Documents\bootcamp-selected\nbanalysis\'; 
addpath(genpath(npmkdir))
addpath(genpath(nbanadir))

directory1  = 'C:\Users\maier\Documents\LGN_data\190326_B_cinterocdrft001_data\';
BRdatafile1 = '190326_B_cinterocdrft001';
filename1   = [directory1 BRdatafile1]; 
addpath(genpath(directory1))

directory2  = 'C:\Users\maier\Documents\LGN_data\190213_B_cinterocdrft001_data\';
BRdatafile2 = '190213_B_cinterocdrft001';
filename2   = [directory2 BRdatafile2]; 
addpath(genpath(directory2))

directory3  = 'C:\Users\maier\Documents\LGN_data\190210_B_cinterocdrft001_data\';
BRdatafile3 = '190210_B_cinterocdrft001';
filename3   = [directory3 BRdatafile3]; 
addpath(genpath(directory3))

directory4  = 'C:\Users\maier\Documents\LGN_data\190208_B_cinterocdrft_data\';
BRdatafile4 = '190208_B_cinterocdrft001';
filename4   = [directory4 BRdatafile4]; 
addpath(genpath(directory4))

directory5  = 'C:\Users\maier\Documents\LGN_data\190120_B_cinterocdrft001_data\';
BRdatafile5 = '190120_B_cinterocdrft001';
filename5   = [directory5 BRdatafile5]; 
addpath(genpath(directory5))

directory6  = 'C:\Users\maier\Documents\LGN_data\190119_B_cinterocdrft002_data\';
BRdatafile6 = '190119_B_cinterocdrft002';
filename6   = [directory6 BRdatafile6]; 
addpath(genpath(directory6))

directory7  = 'C:\Users\maier\Documents\LGN_data\190124_B_cinterocdrft_data\';
BRdatafile7 = '190124_B_cinterocdrft001';
filename7   = [directory7 BRdatafile7]; 
addpath(genpath(directory7))

f = {'_dom0_nondomsup50','_domsup50_nondom0','_domsup50_domsup50'};


aMUA11 = load(strcat(filename1, 'channel17',f{1}, '_aMUA.mat')); 
aMUA12 = load(strcat(filename1, 'channel18',f{1}, '_aMUA.mat')); 
aMUA13 = load(strcat(filename2, 'channel16',f{1}, '_aMUA.mat')); 
aMUA14 = load(strcat(filename2, 'channel17',f{1}, '_aMUA.mat')); 
aMUA15 = load(strcat(filename3, 'channel22',f{1}, '_aMUA.mat'));
aMUA16 = load(strcat(filename4, 'channel14',f{1}, '_aMUA.mat')); 
aMUA17 = load(strcat(filename5, 'channel9',f{1}, '_aMUA.mat')); 
aMUA18 = load(strcat(filename6, 'channel18',f{1}, '_aMUA.mat')); 
aMUA19 = load(strcat(filename6, 'channel20',f{1}, '_aMUA.mat'));
aMUA110 = load(strcat(filename7, 'channel23',f{1}, '_aMUA.mat'));


aMUA21 = load(strcat(filename1, 'channel17',f{2}, '_aMUA.mat')); 
aMUA22 = load(strcat(filename1, 'channel18',f{2}, '_aMUA.mat')); 
aMUA23 = load(strcat(filename2, 'channel16',f{2}, '_aMUA.mat')); 
aMUA24 = load(strcat(filename2, 'channel17',f{2}, '_aMUA.mat')); 
aMUA25 = load(strcat(filename3, 'channel22',f{2}, '_aMUA.mat'));
aMUA26 = load(strcat(filename4, 'channel14',f{2}, '_aMUA.mat')); 
aMUA27 = load(strcat(filename5, 'channel9',f{2}, '_aMUA.mat')); 
aMUA28 = load(strcat(filename6, 'channel18',f{2}, '_aMUA.mat')); 
aMUA29 = load(strcat(filename6, 'channel20',f{2}, '_aMUA.mat'));
aMUA210 = load(strcat(filename7, 'channel23',f{2}, '_aMUA.mat'));

aMUA31 = load(strcat(filename1, 'channel17',f{3}, '_aMUA.mat')); 
aMUA32 = load(strcat(filename1, 'channel18',f{3}, '_aMUA.mat')); 
aMUA33 = load(strcat(filename2, 'channel16',f{3}, '_aMUA.mat')); 
aMUA34 = load(strcat(filename2, 'channel17',f{3}, '_aMUA.mat')); 
aMUA35 = load(strcat(filename3, 'channel22',f{3}, '_aMUA.mat'));
aMUA36 = load(strcat(filename4, 'channel14',f{3}, '_aMUA.mat')); 
aMUA37 = load(strcat(filename5, 'channel9',f{3}, '_aMUA.mat')); 
aMUA38 = load(strcat(filename6, 'channel18',f{3}, '_aMUA.mat')); 
aMUA39 = load(strcat(filename6, 'channel20',f{3}, '_aMUA.mat'));
aMUA310 = load(strcat(filename7, 'channel23',f{3}, '_aMUA.mat'));


DE0_NDE50_aMUA1 = aMUA11.chaMUA;
DE0_NDE50_aMUA2 = aMUA12.chaMUA;
DE0_NDE50_aMUA3 = aMUA13.chaMUA;
DE0_NDE50_aMUA4 = aMUA14.chaMUA;
DE0_NDE50_aMUA5 = aMUA15.chaMUA;
DE0_NDE50_aMUA6 = aMUA16.chaMUA;
DE0_NDE50_aMUA7 = aMUA17.chaMUA;
DE0_NDE50_aMUA8 = aMUA18.chaMUA;
DE0_NDE50_aMUA9 = aMUA19.chaMUA;
DE0_NDE50_aMUA10 = aMUA110.chaMUA;

DE0_NDE50_aMUA = struct();
DE0_NDE50_aMUA.aMUA1 = DE0_NDE50_aMUA1;
DE0_NDE50_aMUA.aMUA2 = DE0_NDE50_aMUA2;
DE0_NDE50_aMUA.aMUA3 = DE0_NDE50_aMUA3;
DE0_NDE50_aMUA.aMUA4 = DE0_NDE50_aMUA4;
DE0_NDE50_aMUA.aMUA5 = DE0_NDE50_aMUA5;
DE0_NDE50_aMUA.aMUA6 = DE0_NDE50_aMUA6;
DE0_NDE50_aMUA.aMUA7 = DE0_NDE50_aMUA7;
DE0_NDE50_aMUA.aMUA8 = DE0_NDE50_aMUA8;
DE0_NDE50_aMUA.aMUA9 = DE0_NDE50_aMUA9;
DE0_NDE50_aMUA.aMUA10 = DE0_NDE50_aMUA10;


DE50_NDE0_aMUA1 = aMUA21.chaMUA;
DE50_NDE0_aMUA2 = aMUA22.chaMUA;
DE50_NDE0_aMUA3 = aMUA23.chaMUA;
DE50_NDE0_aMUA4 = aMUA24.chaMUA;
DE50_NDE0_aMUA5 = aMUA25.chaMUA;
DE50_NDE0_aMUA6 = aMUA26.chaMUA;
DE50_NDE0_aMUA7 = aMUA27.chaMUA;
DE50_NDE0_aMUA8 = aMUA28.chaMUA;
DE50_NDE0_aMUA9 = aMUA29.chaMUA;
DE50_NDE0_aMUA10 = aMUA210.chaMUA;

DE50_NDE0_aMUA = struct();
DE50_NDE0_aMUA.aMUA1 = DE50_NDE0_aMUA1;
DE50_NDE0_aMUA.aMUA2 = DE50_NDE0_aMUA2;
DE50_NDE0_aMUA.aMUA3 = DE50_NDE0_aMUA3;
DE50_NDE0_aMUA.aMUA4 = DE50_NDE0_aMUA4;
DE50_NDE0_aMUA.aMUA5 = DE50_NDE0_aMUA5;
DE50_NDE0_aMUA.aMUA6 = DE50_NDE0_aMUA6;
DE50_NDE0_aMUA.aMUA7 = DE50_NDE0_aMUA7;
DE50_NDE0_aMUA.aMUA8 = DE50_NDE0_aMUA8;
DE50_NDE0_aMUA.aMUA9 = DE50_NDE0_aMUA9;
DE50_NDE0_aMUA.aMUA10 = DE50_NDE0_aMUA10;


DE50_NDE50_aMUA1 = aMUA31.chaMUA;
DE50_NDE50_aMUA2 = aMUA32.chaMUA;
DE50_NDE50_aMUA3 = aMUA33.chaMUA;
DE50_NDE50_aMUA4 = aMUA34.chaMUA;
DE50_NDE50_aMUA5 = aMUA35.chaMUA;
DE50_NDE50_aMUA6 = aMUA36.chaMUA;
DE50_NDE50_aMUA7 = aMUA37.chaMUA;
DE50_NDE50_aMUA8 = aMUA38.chaMUA;
DE50_NDE50_aMUA9 = aMUA39.chaMUA;
DE50_NDE50_aMUA10 = aMUA310.chaMUA;

DE50_NDE50_aMUA = struct();
DE50_NDE50_aMUA.aMUA1 = DE50_NDE50_aMUA1;
DE50_NDE50_aMUA.aMUA2 = DE50_NDE50_aMUA2;
DE50_NDE50_aMUA.aMUA3 = DE50_NDE50_aMUA3;
DE50_NDE50_aMUA.aMUA4 = DE50_NDE50_aMUA4;
DE50_NDE50_aMUA.aMUA5 = DE50_NDE50_aMUA5;
DE50_NDE50_aMUA.aMUA6 = DE50_NDE50_aMUA6;
DE50_NDE50_aMUA.aMUA7 = DE50_NDE50_aMUA7;
DE50_NDE50_aMUA.aMUA8 = DE50_NDE50_aMUA8;
DE50_NDE50_aMUA.aMUA9 = DE50_NDE50_aMUA9;
DE50_NDE50_aMUA.aMUA10 = DE50_NDE50_aMUA10;

fns = fieldnames(DE0_NDE50_aMUA);
xabs = -50:1500;
h = figure;
subplot(length(fns),3,1)
for i = 1:length(fns)
    data = DE0_NDE50_aMUA.(fns{i});
    mean_data = mean(squeeze(data),3);
    norm_mean = (mean_data - min(mean_data))/(max(mean_data)-min(mean_data));
    basedata = norm_mean(25:75);
    mean_bp = mean(basedata,1);
    norm_mean_bscorr = norm_mean - mean_bp;
    subplot(i,3,1)
    plot(xabs, norm_mean_bscorr)
end