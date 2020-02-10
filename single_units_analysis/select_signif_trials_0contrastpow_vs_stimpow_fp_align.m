gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);

xabs = -199:1300;

nyq = 15000;
%starting from 2 as the first one should be excluded
channum = 1: length(data_file.good_data);
mean_filtered_dSUA = struct();
signi = nan(length(channum),1);
mean_blank_power = nan(length(channum),1);
mean_fourhzpowerstim = nan(length(channum),1);

up_dist = nan(1, length(channum));
all_locsdSUA_filtered = nan(1,length(channum));
mean_sua_bsl = nan(length(channum),1);
mean_sua_wnd1 = nan(length(channum),1);
signi_fploc = nan(length(channum),1);



for i = channum  
    contrast = data_file.good_data(i).channel_data.contrast ==  0 & data_file.good_data(i).channel_data.fixedc ==  0;
    
    
    trialidx = 1:length(data_file.good_data(i).channel_data.sdftr_chan(1,:));
    raw_bs = nan(length(xabs), length(trialidx));
    filtered_dSUA = nan(length(xabs), length(trialidx));
    
    powerstim = nan(length(trialidx),1025);
    freqstim = nan(length(trialidx),1025);
    fourhzpowerstim =nan(length(trialidx),1);
    bsl = nan(1, length(trialidx));
    mean_wnd1 = nan(1,length(trialidx));
     for tridx = trialidx
      
    data = data_file.good_data(i).channel_data.sdftr_chan(401:1900,tridx);
      
    bsl(tridx) = mean(data(1:200));
    raw_bs(:,tridx) = data(1:end)- bsl(tridx);
   
    lpc       = 100; %low pass cutoff
    lWn       = lpc/nyq;
    [bwb,bwa] = butter(4,lWn,'low');
    lpdSUA      = filtfilt(bwb,bwa, raw_bs(:,tridx));
  

    filtered_dSUA(:,tridx) = lpdSUA;
    mean_wnd1(tridx) = mean(lpdSUA(231:480));
    
    %%% power

          
    [powerstim(tridx,:), freqstim(tridx,:)] = calcFFT(data(200:1350));
        
    %find the index of the frequency vector closest to 4hz and point to the
    %power value of this index for every trial, and store the value in
    %fourhzpower
    [val,index] = min(abs(4-freqstim(tridx,:)));
    fourhzpowerstim(tridx,1) = powerstim(tridx,index);
    
     end
   %power related variables
   power0 = fourhzpowerstim(contrast);

   mean_blank_power(i,1) = mean(power0);
   mean_fourhzpowerstim(i,1) = mean(fourhzpowerstim);
   
   %first peak location related variables 
    sua_bsl =  mean(filtered_dSUA(1:200,:),1);
    [signi_fploc(i,1), ~] = ttest2(sua_bsl, mean_wnd1);
    mean_sua_bsl(i,1) = mean(sua_bsl);
    mean_sua_wnd1(i,1) = mean(mean_wnd1);
  %%%%%%%%%%% %reject trials below the 95%CI in the blank condition %%%%%%
   for tr = 1:length(fourhzpowerstim)
       if fourhzpowerstim(tr) > mean(power0)+1.96*std(power0) && signi_fploc(i,1) ==1 && mean_sua_bsl(i,1) < mean_sua_wnd1(i,1)
           %fourhzpowerstim(tr) = fourhzpowerstim(tr);
           filtered_dSUA(:,tr) = filtered_dSUA(:,tr);
           
       else
            %fourhzpowerstim(tr) = nan(1,1);
             filtered_dSUA(:,tr) = nan(length(filtered_dSUA(:,tr)),1);
       end
   end
   filtered_dSUA = filtered_dSUA(:,~all(isnan(filtered_dSUA))); % for nan - cols
   [signi(i,1), pvalue] = ttest2(power0,fourhzpowerstim);
   
   if length(filtered_dSUA(1,:)) >=10
   
    %determine the first peak location for each trial of a given single
    %unit
    all_locsdSUA_trials = nan(1,length(filtered_dSUA(1,:)));
    up_dist_trials = nan(1, length(filtered_dSUA(1,:)));
    for trial = 1:length(filtered_dSUA(1,:))
        
         for ln = 30:550
             if filtered_dSUA(200+ln,trial) < filtered_dSUA(200+ln+1,trial)
              
            
             locsdSUA_trial = findpeaks(filtered_dSUA(200+ln:1350,trial));
             %store first peak location 
             all_locsdSUA_trials(:,trial) = locsdSUA_trial.loc(1)+200+ln;
         
             break 
             end 
         end
         up_dist_trials(:,trial)= length(xabs)- all_locsdSUA_trials(trial);
    end
    max_low_dist_unit = max(all_locsdSUA_trials(1,:));
    %create new matrix with the length(max(d)+max(xabs - d))
    new_dist_unit = max_low_dist_unit + max(up_dist_trials); 
    fp_locked_trials = nan(new_dist_unit,length(filtered_dSUA(1,:)));
   
     clear n
     for n = 1:length(filtered_dSUA(1,:))
         if ~isnan(all_locsdSUA_trials(n))
         lower_unit_bound =max_low_dist_unit-all_locsdSUA_trials(n)+1;
         upper_unit_bound =max_low_dist_unit-all_locsdSUA_trials(n)+length(xabs);
     
         fp_locked_trials(lower_unit_bound:upper_unit_bound,n) = filtered_dSUA(:,n);
         end
     end
     fp_locked_trials_out = fp_locked_trials(:,~all(isnan(fp_locked_trials))); % for nan - cols
      %compute the mean single unit activity if more than 10 trials
      mean_filtered_dSUA(i).mean_unit = mean(fp_locked_trials_out,2); % for nan - cols
 
   elseif length(filtered_dSUA(1,:)) <=10
      mean_filtered_dSUA(i).mean_unit = nan(length(fp_locked_trials(:,1)),1);  
   end
   %{
   figure();
   x = 1:length(fp_locked_trials_out(:,1));
   plot(x,fp_locked_trials)
   hold on
   plot(x, mean(fp_locked_trials_out,2),'LineWidth',1, 'Color', 'black')
   
   mean_aligned = mean(fp_locked_trials_out,2);
   nanmean_aligned = nanmean(fp_locked_trials_out,2);
   %}
   %align the mean single units that have significant power at 4hz to the first peak
    if signi(i,1) == 1 &&  mean(power0) < mean(fourhzpowerstim) %find peaks for every unit
    for len = 30:1550
        
        if ~all(isnan(mean_filtered_dSUA(i).mean_unit)) && mean_filtered_dSUA(i).mean_unit(len) < mean_filtered_dSUA(i).mean_unit(len+1)
              
            
             locsdSUA_filtered = findpeaks(mean_filtered_dSUA(i).mean_unit(len:end));
             %store first peak location 
             all_locsdSUA_filtered(:,i) = locsdSUA_filtered.loc(1)+len;
         
             break 
        end 
    end
   else
     %if stimulation power at 4Hz is below baseline
     % generate 4 arbitrary linearly spaced peaks locations 
     locsdSUA_filtered = linspace(1,1150,4); 
     all_locsdSUA_filtered(:,i) = locsdSUA_filtered(2);
                     
    end

    %compute the distance between the first peak and the last datapoint and store  
   %in a matrix
   up_dist(:,i)= length(mean_filtered_dSUA(i).mean_unit)- all_locsdSUA_filtered(i);
  
end        

%get the max distance between the first peak and the stimulus onset (or
%first data point)
 max_low_dist = max(all_locsdSUA_filtered(1,:));
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist);
 clear fp_locked_data norm_fp_locked
 fp_locked_data = nan(round(new_dist),length(channum));
 norm_fp_locked = nan(round(new_dist),length(channum));
 clear n
 for n = 1:length(channum)
    
 lower_bound =round(max_low_dist)-round(all_locsdSUA_filtered(n))+1;
 upper_bound =round(max_low_dist)-round(all_locsdSUA_filtered(n))+length(mean_filtered_dSUA(n).mean_unit);
 % need to change this :
       if signi(n,1) ==1 && mean_blank_power(n,1) < mean_fourhzpowerstim(n,1)
    
 fp_locked_data(lower_bound:upper_bound,n) = mean_filtered_dSUA(n).mean_unit;
  %{
 for len = 30:1500
            if fp_locked_data(len,n) < fp_locked_data(len+1,n)
locsdSUA_mean = findpeaks(fp_locked_data(len:end, n));
        break
            end 
 end
 %x1 = -locsdSUA_mean.loc(1)-len+1:length(fp_locked_data(:,n))-locsdSUA_mean.loc(1)-len;
%}
 

  norm_fp_locked(lower_bound:upper_bound,n) = (fp_locked_data(lower_bound:upper_bound,n)-min(fp_locked_data(lower_bound:upper_bound,n)))/(max(fp_locked_data(lower_bound:upper_bound,n))-min(fp_locked_data(lower_bound:upper_bound,n)));
  
        else
       fp_locked_data(:,n) = nan(length(fp_locked_data(:,n)),1);
       norm_fp_locked(:,n) = nan(length(fp_locked_data(:,n)),1);
       
       end
        %}
  
 end
 
 fp_locked_data_out = fp_locked_data(:,~all(isnan(fp_locked_data))); % for nan - cols
 

figure();
x = 1:length(fp_locked_data_out(:,1));
 plot(x, fp_locked_data_out(:,:),'HandleVisibility','off')
 hold on
 plot(x, mean(fp_locked_data_out(:,:),2),'LineWidth',1, 'Color', 'black')
 plot(xlim, [0 0],'k','HandleVisibility','off')
% plot(xlim, [mean(data_file_bsl) mean(data_file_bsl)],'k')
 %plot(xlim, [mean(data_file_bsl)+1.96*data_file_std_bsl mean(data_file_bsl)+1.96*data_file_std_bsl],'k')
 
  %txtmeanBSL = ['mean baseline'];
  %text(0, mean(data_file_bsl)+3, txtmeanBSL)
  %txt95ci = ['mean baseline +1.96std(bsl)'];
  %text(0, data_file_bsl+1.96*data_file_std_bsl+3, txt95ci)
  title('Single Units Data All Contrasts Involved' , 'Interpreter', 'none')
