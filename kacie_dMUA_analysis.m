

npmkdir    = 'C:\Users\maier\Documents\MATLAB\NPMK-master\'; 
nbanadir   = 'C:\Users\maier\Documents\bootcamp-selected\nbanalysis\'; 

directory  = 'C:\Users\maier\Documents\LGN_data\kacie_preproc\';
STIMBRdatafile = '160616_I_cinterocdrft004_p02_trigAUTO';
STIMfilename   = [directory STIMBRdatafile]; 

directoryc  = 'C:\Users\maier\Documents\LGN_data\kacie_preproc\kacie_channels\';
chanfilename = [directoryc STIMBRdatafile];

directoryp = 'C:\Users\maier\Documents\LGN_data\kacie_preproc\plots';
plotfilename = [directoryp STIMBRdatafile];
addpath(genpath(directory))
addpath(genpath(npmkdir))
addpath(genpath(nbanadir))

STIMdMUA = load(strcat(STIMfilename, '.mat'));  


data_header.contrast1 = STIMdMUA.STIM.contrast ==0 & STIMdMUA.STIM.fixedc >= 0.5;
data_header.contrast2 = STIMdMUA.STIM.contrast >=0.5 & STIMdMUA.STIM.fixedc == 0;
data_header.contrast3 = STIMdMUA.STIM.contrast >=0.5 & STIMdMUA.STIM.fixedc >= 0.5;

f = {'_DE0_NDE50','_DE50_NDE0','_DE50_NDE50'};

dMUA = STIMdMUA.STIM.sdftr;

cont_dMUA = STIMdMUA.STIM.sdftr(1:length(STIMdMUA.STIM.sdftr(:,1,1)),:,data_header.contrast3(1:length(STIMdMUA.STIM.sdftr(1,1,:))));
chdMUA = cont_dMUA(:,24,:);
%save(strcat(chanfilename,'channel24',f{3}, '_dMUA.mat'), 'chdMUA');

 %}


%% plot percent change
 
 h = figure;

xabs = -50:1250;
idx = [1 3 5 7 9 11 13 15 17 19 21 23 2 4 6 8 10 12 14 16 18 20 22 24];
norm_mean = nan(length(xabs), length(cont_dMUA(1,:,1)));
norm_mean_percentch = nan(length(xabs), length(cont_dMUA(1,:,1)));

 for i = 1:length(cont_dMUA(1,:,1))
     
    mean_data = mean(cont_dMUA,3);
    norm_mean = (mean_data(550:1850,i) - min(mean_data(550:1850,i)))/(max(mean_data(550:1850,i))-min(mean_data(550:1850,i)));
    basedata = norm_mean(25:75);
    mean_bp = mean(basedata,1);
    norm_mean_percentch(:,i) = (norm_mean - mean_bp)*100/mean_bp;

    sp = subplot(length(cont_dMUA(1,:,1))/2, 2, idx(i));
    
    plot(xabs, norm_mean_percentch(:, i))
    %set(sp,'position',get(sp,'position').*[1 1 1 scale_factor_y]);
    hold on
    plot([0 0], ylim,'k')
    hold on
    plot([1150 1150], ylim,'k')

    xlim([-50 1250]);
    ylim([-500 1300]);
    
    if i == length(cont_dMUA(1,:,1))/4
        ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}(% change)'});
    end
   
    if i < 12 || (i >= 13)&&(i < 24)
        set(sp, 'XTick', [])
    end
    ylabelh = text(max(xabs), 0, num2str(i),'HorizontalAlignment','left','FontName', 'Arial','FontSize', 13);
           
 end
  sgtitle({STIMBRdatafile, strcat( f{3}, ' percent change')}, 'Interpreter', 'none')
    xlabel('\fontsize{9}time (ms)')
 saveas(gcf,strcat(plotfilename, f{3}, '_plots', '.png'));
    
    
     %% stimuli duration time
  
selected_onsettimes = STIMdMUA.STIM.onsets;
selected_offsettimes = STIMdMUA.STIM.offsets;
     
stimtime1 = selected_onsettimes(1)- selected_offsettimes(1);
stimtime2 = selected_onsettimes(2)- selected_offsettimes(2);
stimtime3 = selected_onsettimes(3)- selected_offsettimes(3);
