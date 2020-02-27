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

clean_high_SUA = struct();
data_peaks = struct();
peaks_locs = struct();
 
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
    bsl = nan(1, length(trialidx));
    mean_wnd1 = nan(1,length(trialidx));
    
    all_pks = nan(4,length(data_file.new_data(i).channel_data.sdftr_chan(1,highcontrast)));
     
     for tridx = trialidx
      
    all_data = data_file.new_data(i).channel_data.sdftr_chan(401:1900,tridx);
    raw_bs(:,tridx) = all_data(1:end)- mean(all_data(1:200));
   
    lpc       = 4.5; %low pass cutoff
    lWn       = lpc/nyq;
    [bwb,bwa] = butter(4,lWn,'low');
    lpdSUA      = filtfilt(bwb,bwa, raw_bs(:,tridx));
  

    filtered_dSUA(:,tridx) = lpdSUA;
    %all_norm_lpdSUA(:,tridx) = (lpdSUA - min(lpdSUA))/(max(lpdSUA)- min(lpdSUA));
    mean_wnd1(tridx) = mean(lpdSUA(201:480));
    
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
   mean_wnd1_5 =mean_wnd1(highcontrast);
   filtered_dSUA_high = filtered_dSUA(:, highcontrast);
   %all_norm_lpdSUA_high = all_norm_lpdSUA(:, highcontrast);
   %first peak location related variables 
   sua_bsl =  mean(filtered_dSUA_high(1:200,:),1);
  
   for tr = 1:length(power5)
      if mean_wnd1_5(tr) > mean(sua_bsl)+1.96*std(sua_bsl)  && power5(tr) > mean(power0)+1.96*std(power0) %
     
           filtered_dSUA_high(:,tr) = filtered_dSUA_high(:,tr);
           
       else
   
           filtered_dSUA_high(:,tr) = nan(length(filtered_dSUA_high(:,tr)),1);
       end
   end
  
      %determine the first peak location for each trial of a given single
    %unit
    all_locsdSUA_trials = nan(6,length(filtered_dSUA_high(1,:)));
    %peak1 = nan(1, length(filtered_dSUA_high(1,:)));
    for trial = 1:length(filtered_dSUA_high(1,:))
        
         for ln = 1:550
             if filtered_dSUA_high(200+ln,trial) < filtered_dSUA_high(200+ln+1,trial) && ~all(isnan(filtered_dSUA_high(:,trial)))
  
             locsdSUA_trial = findpeaks(filtered_dSUA_high(200+ln:1350,trial));
             %if peak1 is too small, peak2 becomes peak1
                     if filtered_dSUA_high(locsdSUA_trial.loc(1)+200+ln,trial) >= 0.4*filtered_dSUA_high(locsdSUA_trial.loc(2)+200+ln)
             %store first peak location 
                     all_locsdSUA_trials(1:length(locsdSUA_trial.loc),trial) = locsdSUA_trial.loc(1:end)+200+ln;
                     else %if filtered_dSUA_high(locsdSUA_trial.loc(1)+200+ln,trial) < 0.4*filtered_dSUA_high(locsdSUA_trial.loc(2)+200+ln) && filtered_dSUA_high(locsdSUA_trial.loc(2)+200+ln,trial) >= 0.4*filtered_dSUA_high(locsdSUA_trial.loc(3)+200+ln)
                     all_locsdSUA_trials(1:length(locsdSUA_trial.loc(2:end)),trial) = locsdSUA_trial.loc(2:end)+200+ln; 
                     
                     end
                     
              break 
             end 
         end
         
        if nnz(~isnan(all_locsdSUA_trials(:,trial))) >= 4 && ~all(isnan(all_locsdSUA_trials(:,trial))) 
             %adjust location to the first data point of lpsu (+ln),
    %lpsulocs = ~isnan(all_locsdSUA_trials(:,trial));
    all_pks(:,trial) = filtered_dSUA_high(all_locsdSUA_trials(1:4,trial), trial);
    filtered_dSUA_high(:,trial) = filtered_dSUA_high(:,trial); 
    all_locsdSUA_trials(:,trial) = all_locsdSUA_trials(:,trial);
   % peak1(1,trial) = filtered_dSUA_high(all_locsdSUA_trials(1,trial), trial);
        else
    filtered_dSUA_high(:,trial) = nan(length(filtered_dSUA_high(:,trial)),1);
    all_locsdSUA_trials(:,trial) = nan(size(all_locsdSUA_trials(:,trial)));
        end 
    end
    %%% reject outlier peaks and the corresponding trials in
    %%% filtered_dSUA_high
   
   
   %reject if there is a peak 1 outlier, if the max peak value in the
   %baseline is an outlier
   
    % First find peaks before stimulus onset 
    
    bsl_peaks = nan(1, length(filtered_dSUA_high(1,:)));
    clear tr
    for tr = 1:length(filtered_dSUA_high(1,:))
        
         for loc = 1:200
            if filtered_dSUA_high(loc,tr) < filtered_dSUA_high(loc+1,tr) && ~all(isnan(filtered_dSUA_high(:,tr)))
           bsl_peak_locs = findpeaks(filtered_dSUA_high(loc:200,tr));
           bsl_peaks(1,tr) = max(filtered_dSUA_high(bsl_peak_locs.loc+loc,tr));
          break
            end 
         end
    end
    
    out_bsl_peaks = isoutlier(bsl_peaks);
    
    p1outliers = isoutlier(all_pks(1,:));
    clear tr
    for tr = 1:length(filtered_dSUA_high(1,:))
        %exclude trials
        if p1outliers(tr) == 0 && ~all(isnan(all_pks(:,tr))) && out_bsl_peaks(tr) ==0 
            
    filtered_dSUA_high(:,tr) = filtered_dSUA_high(:, tr);
    all_pks(:, tr) = all_pks(:,tr);
    all_locsdSUA_trials(:,tr) = all_locsdSUA_trials(:,tr);
   % peak1(tr) = peak1(tr);
        else 
    filtered_dSUA_high(:,tr) = nan(length(filtered_dSUA_high(:,tr)),1);
    all_pks(:,tr) = nan(length(all_pks(:,tr)),1);
    all_locsdSUA_trials(:,tr) = nan(size(all_locsdSUA_trials(:,tr)));
   % peak1(tr) = nan(1,1);
        end
    end
    filtered_dSUA_high = filtered_dSUA_high(:,~all(isnan(filtered_dSUA_high))); % for nan - cols
   all_locsdSUA_trials =  all_locsdSUA_trials(:,~all(isnan(all_locsdSUA_trials)));
   all_pks = all_pks(:, ~all(isnan(all_pks)));
   %peak1 = peak1(:,~isnan(peak1));
   
 %%% perform the peak1 outlier exclusion a second time
 %{
    p1outliers = isoutlier(all_pks(1,:));
    clear tr
    for tr = 1:length(filtered_dSUA_high(1,:))
        %exclude trials
        if  ~all(isnan(all_pks(:,tr))) && p1outliers(tr) == 0 %&& out_bsl_peaks(tr) ==0 
            
    filtered_dSUA_high(:,tr) = filtered_dSUA_high(:, tr);
    all_pks(:, tr) = all_pks(:,tr);
    all_locsdSUA_trials(:,tr) = all_locsdSUA_trials(:,tr);
        else 
    filtered_dSUA_high(:,tr) = nan(length(filtered_dSUA_high(:,tr)),1);
    all_pks(:,tr) = nan(length(all_pks(:,tr)),1);
    all_locsdSUA_trials(:,tr) = nan(size(all_locsdSUA_trials(:,tr)));
        end
    end
 
    
    filtered_dSUA_high = filtered_dSUA_high(:,~all(isnan(filtered_dSUA_high))); % for nan - cols
   all_locsdSUA_trials =  all_locsdSUA_trials(:,~all(isnan(all_locsdSUA_trials)));
   all_pks = all_pks(:, ~all(isnan(all_pks)));
 %}
   %%%%
  
   if length(filtered_dSUA_high(1,:)) >=10
   clean_high_SUA(i).namelist =  filtered_dSUA_high;
   peaks_locs(i).locs = all_locsdSUA_trials; 
   elseif length(filtered_dSUA_high(1,:)) <10  
    all_pks(:,:) = nan(length(all_pks(:,1)),length(all_pks(1,:)));
    clean_high_SUA(i).namelist =  [];
    peaks_locs(i).locs = []; 
   end
   
 data_peaks(i).namelist = all_pks(:,~all(isnan(all_pks)));
 all_pks = all_pks(:,~all(isnan(all_pks)));
channelfilename = [newdatadir 'su_peaks_02262020_2\' filename];
save(strcat(channelfilename, '.mat'), 'all_pks');
 end
 allfilename = [newdatadir 'su_peaks_02262020_2\all_units\all_data_peaks'];
 save(strcat(allfilename, '.mat'), 'data_peaks');
 allfilename = [newdatadir 'su_peaks_02262020_2\all_units\clean_SUA_sup_50'];
 save(strcat(allfilename, '.mat'), 'clean_high_SUA');
 allfilename = [newdatadir 'su_peaks_02262020_2\all_units\clean_SUA_locs'];
 save(strcat(allfilename, '.mat'), 'peaks_locs');
 
   cnt =0;
   for i =1:length(data_peaks)
       if ~isempty(data_peaks(i).namelist)
           cnt = cnt+1;
       end
        
   end
   