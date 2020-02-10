gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);

xabs = -199:1300;

nyq = 15000;
%starting from 2 as the first one should be excluded
channum = 1: length(data_file.good_data);
mean_filtered_dSUA = nan(length(xabs), length(channum));
signi = nan(length(channum),1);
mean_bs_power = nan(length(channum),1);
mean_fourhzpowerstim = nan(length(channum),1);

up_dist = nan(1, length(channum));
all_locsdSUA_filtered = nan(1,length(channum));





for i = channum  
    
    trialidx = 1:length(data_file.good_data(i).channel_data.sdftr_chan(1,:));
    raw_bs = nan(length(xabs), length(trialidx));
    filtered_dSUA = nan(length(xabs), length(trialidx));
    
    powerbs = nan(length(trialidx),513);
    freqbs = nan(length(trialidx),513);
    fourhzpowerbs =nan(length(trialidx),1);
    
    powerstim = nan(length(trialidx),1025);
    freqstim = nan(length(trialidx),1025);
    fourhzpowerstim =nan(length(trialidx),1);
    bsl = nan(1, length(trialidx));

     for tridx = trialidx
      
    data = data_file.good_data(i).channel_data.sdftr_chan(401:1900,tridx);
    bs_data = data_file.good_data(i).channel_data.sdftr_chan(1:600,tridx);
    
    bsl(tridx) = mean(data(1:200));
    raw_bs(:,tridx) = data(1:end)- bsl(tridx);
   
    lpc       = 100; %low pass cutoff
    lWn       = lpc/nyq;
    [bwb,bwa] = butter(4,lWn,'low');
    lpdSUA      = filtfilt(bwb,bwa, raw_bs(:,tridx));
  

    filtered_dSUA(:,tridx) = lpdSUA;
 
    
    %%% power
    % baseline
        [powerbs(tridx,:), freqbs(tridx,:)] = calcFFT(squeeze(bs_data));
        
         %find the index of the frequency vector closest to 4hz and point to the
         %power value of this index for every trial, and store the value in
         %fourhzpower
          [val,index] = min(abs(4-freqbs(tridx,:)));
          fourhzpowerbs(tridx,1) = powerbs(tridx,index);
          
       % stimulation data 
          
         [powerstim(tridx,:), freqstim(tridx,:)] = calcFFT(data(200:1350));
        
    %find the index of the frequency vector closest to 4hz and point to the
    %power value of this index for every trial, and store the value in
    %fourhzpower
    [val,index] = min(abs(4-freqstim(tridx,:)));
    fourhzpowerstim(tridx,1) = powerstim(tridx,index);
    
     end
   
     bs_power = fourhzpowerbs;
    [signi(i,1), pvalue] = ttest2(fourhzpowerbs,fourhzpowerstim);
   mean_bs_power(i,1) = mean(fourhzpowerbs);
   mean_fourhzpowerstim(i,1) = mean(fourhzpowerstim);
    
    mean_filtered_dSUA(:,i) = mean(filtered_dSUA,2);
       if signi(i,1) == 1 &&  mean(bs_power) < mean(fourhzpowerstim)
 %find peaks for every unit
    for len = 30:550
        
        if mean_filtered_dSUA(200+len,i) < mean_filtered_dSUA(200+len+1,i)
              
            
             locsdSUA_filtered = findpeaks(mean_filtered_dSUA(200+len:1350,i));
             %store first peak location 
             all_locsdSUA_filtered(:,i) = locsdSUA_filtered.loc(1)+200+len;
         
             break 
        end 
    end
   else
     %if stimulation power at 4Hz is below baseline
     % generate 4 arbitrary linearly spaced peaks locations 
     locsdSUA_filtered = linspace(230,1350,4); 
     all_locsdSUA_filtered(:,i) = locsdSUA_filtered(2);
                     
      end

    %compute the distance between the first peak and the last datapoint and store  
   %in a matrix
   up_dist(:,i)= length(xabs)- all_locsdSUA_filtered(i);
  
end        

%get the max distance between the first peak and the stimulus onset
 max_low_dist = max(all_locsdSUA_filtered(1,:));
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist);
 clear fp_locked_data norm_fp_locked
 fp_locked_data = nan(round(new_dist),length(channum));
 norm_fp_locked = nan(round(new_dist),length(channum));
 for n = 1:length(channum)
    
 lower_bound =round(max_low_dist)-round(all_locsdSUA_filtered(n))+1;
 upper_bound =round(max_low_dist)-round(all_locsdSUA_filtered(n))+length(xabs);
 % need to change this :
       if signi(n,1) ==1 && mean_bs_power(n,1) < mean_fourhzpowerstim(n,1)
    
 fp_locked_data(lower_bound:upper_bound,n) = mean_filtered_dSUA(:,n);
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
