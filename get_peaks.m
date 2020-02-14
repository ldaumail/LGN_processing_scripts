%after saving the data with new_data_set.m, we isolate the peaks for the
%new analysis of the refined data

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

mean_blank_power = nan(length(channum),1);
mean_fourhzpowerhighcontrast = nan(length(channum),1);

all_locsdSUA_filtered = nan(1,length(channum));
norm_data_peaks = struct();
 
 for i = channum
 data = squeeze(data_file.new_data(i).channel_data.hypo{1,2}.cont_su(401:1900,:,:));
 
 filename = [data_file.new_data(i).channel_data.filename, f{2}];
 filename = erase(filename, '.mat');

 all_norm_pks = nan(4,length(data(1,:)));
 
 
   blankcontrast = data_file.new_data(i).channel_data.contrast ==  0 & data_file.new_data(i).channel_data.fixedc ==  0;
   highcontrast = data_file.new_data(i).channel_data.contrast >=  0.5 & data_file.new_data(i).channel_data.fixedc ==  0; 
    
    trialidx = 1:length(data_file.new_data(i).channel_data.sdftr_chan(1,:));
    raw_bs = nan(length(xabs), length(trialidx));
    filtered_dSUA = nan(length(xabs), length(trialidx));
    all_norm_lpdSUA= nan(length(xabs),length(trialidx));

    
    powerstim = nan(length(trialidx),1025);
    freqstim = nan(length(trialidx),1025);
    fourhzpowerstim =nan(length(trialidx),1);
    bsl = nan(1, length(trialidx));
    mean_wnd1 = nan(1,length(trialidx));
     for tridx = trialidx
      
    all_data = data_file.new_data(i).channel_data.sdftr_chan(401:1900,tridx);
      
    bsl(tridx) = mean(all_data(1:200));
    raw_bs(:,tridx) = all_data(1:end)- bsl(tridx);
   
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
   %power related variables
   power0 = fourhzpowerstim(blankcontrast);
   power5 = fourhzpowerstim(highcontrast);
   mean_blank_power(i,1) = mean(power0);
   mean_fourhzpowerhighcontrast(i,1) = mean(power5);
   
   %spiking activity related variables
   mean_wnd1_5 =mean_wnd1(highcontrast);
   filtered_dSUA_high = filtered_dSUA(:, highcontrast);
   all_norm_lpdSUA_high = all_norm_lpdSUA(:, highcontrast);
   %first peak location related variables 
    sua_bsl =  mean(filtered_dSUA_high(1:200,:),1);
    

  %%%%%%%%%%% %reject trials below the 95%CI in the blank condition %%%%%%
   for tr = 1:length(power5)
       if power5(tr) > mean(power0)+1.96*std(power0) && mean_wnd1_5(tr) > mean(sua_bsl) + 1.96*std(sua_bsl)
     
           filtered_dSUA_high(:,tr) = filtered_dSUA_high(:,tr);
           
       else
   
           filtered_dSUA_high(:,tr) = nan(length(filtered_dSUA_high(:,tr)),1);
       end
   end
   filtered_dSUA_high = filtered_dSUA_high(:,~all(isnan(filtered_dSUA_high))); % for nan - cols
   
   
   all_locsdSUA_trials = nan(6,length(filtered_dSUA_high(1,:)));
   if length(filtered_dSUA_high(1,:)) >=10
   
    %determine the first peak location for each trial of a given single
    %unit
    
  
    for trial = 1:length(filtered_dSUA_high(1,:))
        
         for ln = 30:550
             if filtered_dSUA_high(200+ln,trial) < filtered_dSUA_high(200+ln+1,trial)
              
            
             locsdSUA_trial = findpeaks(filtered_dSUA_high(200+ln:1350,trial));
             %store 4 peaks locations 
             all_locsdSUA_trials(1:length(locsdSUA_trial.loc),trial) = locsdSUA_trial.loc+200+ln;
             break 
             end 
         end
         
        if length(locsdSUA_trial.loc) >= 4
             %adjust location to the first data point of lpsu (+len),
    lpsulocs = locsdSUA_trial.loc(1:4) + 200+ ln;
    all_norm_pks(:,trial) = all_norm_lpdSUA_high(lpsulocs(1:4), trial);
  
         end 
    end
    elseif length(filtered_dSUA_high(1,:)) <=10  || length(locsdSUA_trial.loc) <= 4
      all_locsdSUA_trials(:,:) = nan(length(all_locsdSUA_trials(:,1)),length(all_locsdSUA_trials(1,:)));
      all_norm_pks(:,:) = nan(length(all_norm_pks(:,1)),length(all_norm_pks(1,:)));
   end
 namelist(1,1:length(sprintf('chan_%d',i))) = sprintf('chan_%d',i);
 norm_data_peaks(i).namelist = all_norm_pks(:,~all(isnan(all_norm_pks)));
 all_norm_pks = all_norm_pks(:,~all(isnan(all_norm_pks)));
channelfilename = [newdatadir 'su_peaks\' filename];
save(strcat(channelfilename, '.mat'), 'all_norm_pks');
 end
  allfilename = [newdatadir 'all_norm_data_peaks'];
 save(strcat(allfilename, '.mat'), 'norm_data_peaks');
 
    %{
 %all_locs =nan(4,length(channel_data(1,1,:)));
 for n =1:length(data(1,:))
   lpc       = 4.5; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   lpsu      = filtfilt(bwb,bwa, data(:,n));
   
   norm_lpsu = (lpsu - min(lpsu))/(max(lpsu)- min(lpsu));
   %all_lpsu(:,n) = lpsu;
  % plot(lpsu)
 if ~all(isnan(norm_lpsu)) 
  for len = 30:550
            if norm_lpsu(len) < norm_lpsu(len+1)
   locs = findpeaks(norm_lpsu(len:1200));
        break
            end
  end
         
         if length(locs.loc) >= 4
             %adjust location to the first data point of lpsu (+len),
    lpsulocs = locs.loc(1:4) + len;
   all_norm_pks(:,n) = norm_lpsu(lpsulocs(1:4));
  
         end 
 end
% plot(locs.loc +len, lpsu(locs.loc +len))
 end
 
 namelist(1,1:length(sprintf('chan_%d',i))) = sprintf('chan_%d',i);
 norm_data_peaks(i).namelist = all_norm_pks;
channelfilename = [newdatadir filename];
%save(strcat(channelfilename, '.mat'), 'all_pks');
 end
    %}
 allfilename = [newdatadir 'all_norm_data_peaks'];
 save(strcat(allfilename, '.mat'), 'norm_data_peaks');