gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);

xabs = -199:1300;

nyq = 15000;
mean_filtered_dSUA = nan(length(xabs), length(channum));
signi = nan(length(channum),1);
mean_data_file_bsl = nan(length(channum),1);
mean_bsl = nan(length(channum),1);
up_dist = nan(1, length(channum));
all_locsdSUA_filtered = nan(1,length(channum));

%starting from 2 as the first one should be excluded
channum = 1: length(data_file.good_data);


for i = channum  
    data_file.good_data(i).channel_data.bsls =mean(data_file.good_data(i).channel_data.sdftr_chan(401:600,:),1);
    trialidx = 1:length(data_file.good_data(i).channel_data.sdftr_chan(1,:));
    raw_bs = nan(length(xabs), length(trialidx));
    filtered_dSUA = nan(length(xabs), length(trialidx));
    mean_wnd1 = nan(1,length(trialidx));
    
    bsl = nan(1, length(trialidx));

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
     end
   
    data_file_bsl =  mean(filtered_dSUA(1:200,:),1);
    [signi(i,1), pvalue] = ttest2(data_file_bsl, mean_wnd1);
    mean_data_file_bsl(i,1) = mean(data_file_bsl);
    mean_bsl(i,1) = mean(mean_wnd1);
    
    mean_filtered_dSUA(:,i) = mean(filtered_dSUA,2);
       if signi(i,1) == 1 &&  mean(data_file_bsl) < mean(mean_wnd1)
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
     %if there is no peak in the first window where first peak
     %is expected, generate 4 arbitrary linearly spaced peaks locations 
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
        %if signi(n,1) ==1 && mean_data_file_bsl(n,1) < mean_bsl(n,1)
    
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
  %{
         else
       fp_locked_data(:,n) = nan(length(fp_locked_data(:,n)),1);
       norm_fp_locked(:,n) = nan(length(fp_locked_data(:,n)),1);
       
         end
        %}
  
 end
 
 fp_locked_data_out = fp_locked_data(:,~all(isnan(fp_locked_data))); % for nan - cols
 
for i = 1:length(fp_locked_data(1,:))
figure();
x = 1:length(fp_locked_data(:,1));
 plot(x, fp_locked_data(:,i),'HandleVisibility','off')
 hold on
 plot(x, mean(fp_locked_data(i,:),2),'LineWidth',1, 'Color', 'black')
 plot(xlim, [0 0],'k','HandleVisibility','off')
 plot(xlim, [mean(data_file_bsl) mean(data_file_bsl)],'k')
 %plot(xlim, [mean(data_file_bsl)+1.96*data_file_std_bsl mean(data_file_bsl)+1.96*data_file_std_bsl],'k')
 
  txtmeanBSL = ['mean baseline'];
  text(0, mean(data_file_bsl)+3, txtmeanBSL)
  %txt95ci = ['mean baseline +1.96std(bsl)'];
  %text(0, data_file_bsl+1.96*data_file_std_bsl+3, txt95ci)
  title('Single Units Data All Contrasts Involved' , 'Interpreter', 'none')
end