%after saving the data with new_data_set.m, we isolate the peaks for the
%new analysis of the refined data
%we also save the data we want to plot (only the clean data)

newdatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
channelfilename = [newdatadir 'refined_dataset']; 
data_file = load(channelfilename);

%exclude 160517, (first unit, left empty, it is a K neuron)
%Reject 180806 p1 uclust17, M cell, as doesn't seem well triggered (46)
%Reject 181207 (B) uclust22, M cell, as doesn't seem well triggered (55)
 layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = [];
 f = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
 
 %% find peaks in order to analyze them on R and fit a LMER
 
channum = 1: length(data_file.new_data);
xabs = -199:1300;
nyq = 500;
mean_filtered_dSUA = struct();

data_peaks = struct();
clean_high_SUA = struct();
 
 for i = channum
 
   filename = [data_file.new_data(i).channel_data.filename, f{2}];
   filename = erase(filename, '.mat');

   blankcontrast = data_file.new_data(i).channel_data.contrast ==  0 & data_file.new_data(i).channel_data.fixedc ==  0;
   highcontrast = data_file.new_data(i).channel_data.contrast >=  0.5 & data_file.new_data(i).channel_data.fixedc ==  0; 
    
    trialidx = 1:length(data_file.new_data(i).channel_data.sdftr_chan(1,:));
    raw_bs = nan(length(xabs), length(trialidx));
    filtered_dSUA = nan(length(xabs), length(trialidx));
    all_norm_lpdSUA= nan(length(xabs),length(trialidx));
       
    
    powerstim = nan(length(trialidx),1025);
    freqstim = nan(length(trialidx),1025);
    fourhzpowerstim =nan(length(trialidx),1);
    mean_wnd1 = nan(1,length(trialidx));
    
     all_norm_pks = nan(4,length(data_file.new_data(i).channel_data.sdftr_chan(1,highcontrast)));
     
     for tridx = trialidx
      
    all_data = data_file.new_data(i).channel_data.sdftr_chan(401:1900,tridx);
 
    raw_bs(:,tridx) = all_data(1:end)- mean(all_data(1:200));
   
    lpc       = 4.5; %low pass cutoff
    lWn       = lpc/nyq;
    [bwb,bwa] = butter(4,lWn,'low');
    lpdSUA      = filtfilt(bwb,bwa, raw_bs(:,tridx));
  

    filtered_dSUA(:,tridx) = lpdSUA;
    all_norm_lpdSUA(:,tridx) = (lpdSUA - min(lpdSUA))/(max(lpdSUA)- min(lpdSUA));
    mean_wnd1(tridx) = mean(lpdSUA(231:480));
    
    %%% power

          
    [powerstim(tridx,:), freqstim(tridx,:)] = calcFFT(all_data(200:1350));
        
    %find the index of the frequency vector closest to 4hz and point to the
    %power value of this index for every trial, and store the value in
    %fourhzpower
    [val,index] = min(abs(4-freqstim(tridx,:)));
    fourhzpowerstim(tridx,1) = powerstim(tridx,index);
    
     end
     
   %%%%%%%%%%% %reject trials below the 95%CI in the blank condition %%%%%%
   %power related variables
   power0 = fourhzpowerstim(blankcontrast);
   power5 = fourhzpowerstim(highcontrast);
  
   %spiking activity related variables

   filtered_dSUA_high = filtered_dSUA(:, highcontrast);
   all_norm_lpdSUA_high = all_norm_lpdSUA(:, highcontrast);
   %first peak location related variables
   mean_wnd1_5 = mean_wnd1(highcontrast); 
   sua_bsl =  mean(filtered_dSUA_high(1:200,:),1);
    
   for tr = 1:length(power5)
       if power5(tr) > mean(power0)+1.96*std(power0) && mean_wnd1_5(tr) > mean(sua_bsl) + 1.96*std(sua_bsl)
     
           all_norm_lpdSUA_high(:,tr) = all_norm_lpdSUA_high(:,tr);
           
       else
   
           all_norm_lpdSUA_high(:,tr) = nan(length(all_norm_lpdSUA_high(:,tr)),1);
       end
   end
   all_norm_lpdSUA_high = all_norm_lpdSUA_high(:,~all(isnan(all_norm_lpdSUA_high))); % for nan - cols
   
   
   if length(all_norm_lpdSUA_high(1,:)) >=10
   
    %determine the first peak location for each trial of a given single
    %unit
    
  
    for trial = 1:length(all_norm_lpdSUA_high(1,:))
        
         for ln = 30:550
             if all_norm_lpdSUA_high(200+ln,trial) < all_norm_lpdSUA_high(200+ln+1,trial)
              
            
             locsdSUA_trial = findpeaks(all_norm_lpdSUA_high(200+ln:1350,trial));
         
             break 
             end 
         end
         
        if length(locsdSUA_trial.loc) >= 4
             %adjust location to the first data point of lpsu (+ln),
    lpsulocs = locsdSUA_trial.loc(1:4) + 200+ ln;
    all_norm_pks(:,trial) = all_norm_lpdSUA_high(lpsulocs(1:4), trial);
  
        end 
    end
   namelist(1,1:length(sprintf('chan_%d',i))) = sprintf('chan_%d',i);
   clean_high_SUA(i).namelist =  all_norm_lpdSUA_high;
    
    elseif length(all_norm_lpdSUA_high(1,:)) <10  || length(locsdSUA_trial.loc) <= 4
       all_norm_pks(:,:) = nan(length(all_norm_pks(:,1)),length(all_norm_pks(1,:)));
       clean_high_SUA(i).namelist =  nan(size(all_norm_lpdSUA_high));
   end
   
 namelist(1,1:length(sprintf('chan_%d',i))) = sprintf('chan_%d',i);
 data_peaks(i).namelist = all_norm_pks(:,~all(isnan(all_norm_pks)));
 all_norm_pks = all_norm_pks(:,~all(isnan(all_norm_pks)));
channelfilename = [newdatadir 'norm_su_peaks\' filename];
%save(strcat(channelfilename, '.mat'), 'all_norm_pks');
 end
  allfilename = [newdatadir 'all_norm_data_peaks'];
 %save(strcat(allfilename, '.mat'), 'data_peaks');
 allfilename = [newdatadir 'clean_norm_SUA_sup_50'];
 %save(strcat(allfilename, '.mat'), 'clean_high_SUA');
 