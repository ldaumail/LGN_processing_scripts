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


data = struct();
data.DE0_NDE50_aMUA = DE0_NDE50_aMUA;
data.DE50_NDE0_aMUA = DE50_NDE0_aMUA; 
data.DE50_NDE50_aMUA =DE50_NDE50_aMUA;

%fns1 = fieldnames(data);

%% plot normalized and baseline corrected data
  

h = figure;

xabs = -50:1500;
scale_factor_y = 1.2;
fns = fieldnames(DE0_NDE50_aMUA);
norm_mean_bscorr = nan(length(xabs), length(fns));

 for i = 1:length(fns)
     
    data = DE0_NDE50_aMUA.(fns{i});
    mean_data = mean(squeeze(data),3);
    norm_mean = (mean_data - min(mean_data))/(max(mean_data)-min(mean_data));
    basedata = norm_mean(25:75);
    mean_bp = mean(basedata,1);
    norm_mean_bscorr(:,i) = norm_mean - mean_bp;
    
    sp = subplot(length(fns), 1, i);
    plot(xabs, norm_mean_bscorr(:, i))
    %set(sp,'position',get(sp,'position').*[1 1 1 scale_factor_y]);
    hold on
    plot([0 0], ylim,'k')
    hold on
    plot([1150 1150], ylim,'k')

    xlim([-50 1500]);
    ylim([-0.2 0.2]);
    
    if i == 1
    title('DE0_NDE50_aMUA', 'Interpreter', 'none')
    end
    if i == length(fns)/2
        ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}(microV)'});
    end
    if i <= length(fns)
        set(subplot(length(fns),1,i), 'XTick', [])
    end
   
 end
    xlabel('\fontsize{9}time (ms)')

 %% plot percent change
 
 h = figure;

xabs = -50:1500;
scale_factor_y = 1.2;
fns = fieldnames(DE0_NDE50_aMUA);
norm_mean_percentch = nan(length(xabs), length(fns));

 for i = 1:length(fns)
     
    data = DE0_NDE50_aMUA.(fns{i});
    mean_data = mean(squeeze(data),2);
    norm_mean = (mean_data - min(mean_data))/(max(mean_data)-min(mean_data));
    basedata = norm_mean(25:75);
    mean_bp = mean(basedata,1);
    norm_mean_percentch(:,i) = (norm_mean - mean_bp)/mean_bp;
    
    sp = subplot(length(fns), 1, i);
    plot(xabs, norm_mean_percentch(:, i))
    %set(sp,'position',get(sp,'position').*[1 1 1 scale_factor_y]);
    hold on
    plot([0 0], ylim,'k')
    hold on
    plot([1150 1150], ylim,'k')

    xlim([-50 1500]);
    ylim([-1 3]);
    
    if i == 1
    title({'DE0_NDE50_aMUA', 'percent change'}, 'Interpreter', 'none')
    end
    if i == length(fns)/2
        ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}(no unit)'});
    end
    if i <= length(fns)
        set(subplot(length(fns),1,i), 'XTick', [])
    end
   
 end
    xlabel('\fontsize{9}time (ms)')
    
    %% filtered data percent change and linear regression

h = figure;

xabs = -50:1500;
nyq = 15000;
fns = fieldnames(DE50_NDE50_aMUA);
norm_mean_percentch = nan(length(xabs), length(fns));

clear i ;
 for i = 1:length(fns)
     
   mean_data = mean(squeeze(DE50_NDE50_aMUA.(fns{i})),2);
   norm_mean = (mean_data - min(mean_data))/(max(mean_data)-min(mean_data));
   basedata = norm_mean(25:75);
   mean_bp = mean(basedata,1);
   norm_mean_percentch(:,i) = (norm_mean - mean_bp)/mean_bp;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpaMUA      = filtfilt(bwb,bwa, norm_mean_percentch(:,i));
 
   
   %find maxima of the filtered data and plot linear regression
   [pksaMUA, locsaMUA] = findpeaks(lpaMUA(50:1201));
   linreg1 = polyfit(locsaMUA, pksaMUA, 1);
   
   
   ylinreg1 = linreg1(1) .* locsaMUA + linreg1(2);
   sp = subplot(length(fns), 1, i);
    plot(xabs,lpaMUA)
    hold on
    %plot(locsaMUA, pksaMUA)
    hold on
    plot(locsaMUA, ylinreg1)
    txt1 = ['y = (' num2str(linreg1(1)) ')x + (' num2str(linreg1(2)) ')'];
    text(locsaMUA(1), ylinreg1(1), txt1);
    %set(sp,'position',get(sp,'position').*[1 1 1 scale_factor_y]);
    hold on
    plot([0 0], ylim,'k')
    hold on
    plot([1150 1150], ylim,'k')

    xlim([-50 1500]);
    ylim([-1 2]);
    
    if i == 1
    title({'DE50_NDE50_aMUA', 'percent change'}, 'Interpreter', 'none')
    end
    if i == length(fns)/2
        ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}(no unit)'});
    end
    if i <= length(fns)
        set(subplot(length(fns),1,i), 'XTick', [])
    end
   
 end
    xlabel('\fontsize{9}time (ms)')
    
%% plot all slopes together (filtered data)

h = figure;

xabs = -50:1500;
nyq = 15000;
fns = fieldnames(DE50_NDE0_aMUA);
norm_mean_percentch = nan(length(xabs), length(fns));
linreg_all = nan(2, length(fns));
clear i ;
 for i = 1:length(fns)
     
   mean_data = mean(squeeze(DE50_NDE0_aMUA.(fns{i})),2);
   norm_mean = (mean_data - min(mean_data))/(max(mean_data)-min(mean_data));
   basedata = norm_mean(25:75);
   mean_bp = mean(basedata,1);
   norm_mean_percentch(:,i) = (norm_mean - mean_bp)/mean_bp;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpaMUA      = filtfilt(bwb,bwa, norm_mean_percentch(:,i));
 
   
   %find maxima of the filtered data and plot linear regression
   [pksaMUA, locsaMUA] = findpeaks(lpaMUA(50:1201));
   linreg1 = polyfit(locsaMUA, pksaMUA, 1);
   linreg_all(:,i) = linreg1;
   
   ylinreg1 = linreg1(1) .* locsaMUA + linreg1(2);
plot(locsaMUA, ylinreg1)
hold on
    if i == 1
    title({'DE50_NDE0_aMUA', 'linreg on maxima: all slopes'}, 'Interpreter', 'none')
    end
 end
 
 
 %% percent change and linear regression (no filter)
    
h = figure;

xabs = -50:1500;
nyq = 15000;
fns = fieldnames(DE0_NDE50_aMUA);
norm_mean_percentch = nan(length(xabs), length(fns));
linreg_all = nan(2, length(fns));
clear i ;
 for i = 1:length(fns)
     
   mean_data = mean(squeeze(DE0_NDE50_aMUA.(fns{i})),2);
   norm_mean = (mean_data - min(mean_data))/(max(mean_data)-min(mean_data));
   basedata = norm_mean(25:75);
   mean_bp = mean(basedata,1);
   norm_mean_percentch(:,i) = (norm_mean - mean_bp)/mean_bp;
   %{
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpaMUA      = filtfilt(bwb,bwa, norm_mean_percentch(:,i));
 %}
   
   %find maxima of the filtered data and plot linear regression
   [pksaMUA, locsaMUA] = findpeaks(norm_mean_percentch(50:1201,i));
   linreg1 = polyfit(locsaMUA, pksaMUA, 1);
   linreg_all(:,i) = linreg1;
   
   ylinreg1 = linreg1(1) .* locsaMUA + linreg1(2);
   sp = subplot(length(fns), 1, i);
    plot(xabs, norm_mean_percentch(:,i))
    hold on
    %plot(locsaMUA, pksaMUA)
    hold on
    plot(locsaMUA, ylinreg1)
    txt1 = ['y = (' num2str(linreg1(1)) ')x + (' num2str(linreg1(2)) ')'];
    text(locsaMUA(1), ylinreg1(1), txt1);
    %set(sp,'position',get(sp,'position').*[1 1 1 scale_factor_y]);
    hold on
    plot([0 0], ylim,'k')
    hold on
    plot([1150 1150], ylim,'k')

    xlim([-50 1500]);
    ylim([-1 2]);
    
    if i == 1
    title({'DE0_NDE50_aMUA', 'percent change'}, 'Interpreter', 'none')
    end
    if i == length(fns)/2
        ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}(no unit)'});
    end
    if i <= length(fns)
        set(subplot(length(fns),1,i), 'XTick', [])
    end
   
 end
    xlabel('\fontsize{9}time (ms)')
    
    %% plot all slopes together (no filter)

h = figure;
unitnames = {'190326','190326','190213','190213','190210','190208','190120', ...
    '190119','190119','190124', 'mean'};
xabs = -50:1500;
xabs2 = 0:1150;
nyq = 15000;
fns = fieldnames(DE50_NDE0_aMUA);
norm_mean_percentch = nan(length(xabs), length(fns));
linreg_all = nan(2, length(fns));
mean_linreg = nan(2,1);
clear i ;
 for i = 1:length(fns)
     
   mean_data = mean(squeeze(DE50_NDE0_aMUA.(fns{i})),2);
   norm_mean = (mean_data - min(mean_data))/(max(mean_data)-min(mean_data));
   basedata = norm_mean(25:75);
   mean_bp = mean(basedata,1);
   norm_mean_percentch(:,i) = (norm_mean - mean_bp)/mean_bp;
   %{
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpaMUA      = filtfilt(bwb,bwa, norm_mean_percentch(:,i));
 %}
   
   %find maxima of the filtered data and plot linear regression
   [pksaMUA, locsaMUA] = findpeaks(norm_mean_percentch(50:1201, i));
   linreg1 = polyfit(locsaMUA, pksaMUA, 1);
   linreg_all(:,i) = linreg1;
   
   ylinreg1 = linreg1(1) .* locsaMUA + linreg1(2);
  
plot(locsaMUA, ylinreg1)
txt1 = [unitnames{i}];
    text(locsaMUA(1), ylinreg1(1), txt1);
hold on

    if i == 1
    title({'DE50_NDE0_aMUA', 'linreg on maxima: all slopes'}, 'Interpreter', 'none')
    end
 end
 mean_linreg = mean(linreg_all, 2);
 ymeanlin = mean_linreg(1) .*xabs2 + mean_linreg(2);
 hold on
 plot(xabs2,ymeanlin)
 txt1 = [unitnames{11}  'y = (' num2str(mean_linreg(1)) ')x + (' num2str(mean_linreg(2)) ')'];
    text(xabs2(200), ymeanlin(200), txt1)
legend('190326','190326','190213','190213','190210','190208','190120', ...
    '190119','190119','190124', 'mean')
%% same with filtered data
h = figure;
unitnames = {'190326','190326','190213','190213','190210','190208','190120', ...
    '190119','190119','190124', 'mean'};
xabs = -50:1500;
xabs2 = 0:1150;
nyq = 15000;
fns = fieldnames(DE50_NDE50_aMUA);
norm_mean_percentch = nan(length(xabs), length(fns));
linreg_all = nan(2, length(fns));
pvalues = nan(1,length(fns));
mean_linreg = nan(2,1);
clear i ;
 for i = 1:length(fns)
     
   mean_data = mean(squeeze(DE50_NDE50_aMUA.(fns{i})),2);
   norm_mean = (mean_data - min(mean_data))/(max(mean_data)-min(mean_data));
   basedata = norm_mean(25:75);
   mean_bp = mean(basedata,1);
   norm_mean_percentch(:,i) = (norm_mean - mean_bp)/mean_bp;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpaMUA      = filtfilt(bwb,bwa, norm_mean_percentch(:,i));
 
   
   %find maxima of the filtered data and plot linear regression
   [pksaMUA, locsaMUA] = findpeaks(lpaMUA(50:1201));
   linreg1 = fitlm(locsaMUA, pksaMUA, 'y ~ x1');
   linreg_all(:,i) = table2array(linreg1.Coefficients(1:2,1));
   pvalues(i) = table2array(linreg1.Coefficients(2,4));
   ylinreg1 = linreg_all(2,i) .* locsaMUA + linreg_all(1,i);
  
plot(locsaMUA, ylinreg1)
txt1 = [unitnames{i}];
    text(locsaMUA(1), ylinreg1(1), txt1);
hold on

    if i == 1
    title({'DE50_NDE50_aMUA', 'linreg on maxima: all slopes'}, 'Interpreter', 'none')
    end
 end
 mean_linreg = mean(linreg_all, 2);
 ymeanlin = mean_linreg(2) .*xabs2 + mean_linreg(1);
 hold on
 plot(xabs2,ymeanlin)
 txt1 = [unitnames{11}  'y = (' num2str(mean_linreg(2)) ')x + (' num2str(mean_linreg(1)) ')'];
    text(xabs2(200), ymeanlin(200), txt1)
legend('190326','190326','190213','190213','190210','190208','190120', ...
    '190119','190119','190124', 'mean')

