

gooddatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\';
channelfilename = [gooddatadir 'good_single_units_data_4bmpmore']; 
data_file = load(channelfilename);


channum = 1: length(data_file.good_data); 

C_DE = nan(6, length(channum));
FC_NDE = nan(6, length(channum));

%create logical vectors necessary to select the trial corresponding to
%specific contrasts
data_header = [];
for i = 1:length(channum)
C_DE(1:length(unique(data_file.good_data(i).channel_data.contrast)),i) = unique(data_file.good_data(i).channel_data.contrast);
%FC_NDE(1:length(unique(data_file.good_data(i).channel_data.fixedc)),i) = unique(data_file.good_data(i).channel_data.fixedc);
%let's see the cases where contrast NDE =0
FC_NDE(1:length(unique(data_file.good_data(i).channel_data.contrast)),i) = zeros(length(unique(data_file.good_data(i).channel_data.contrast)),1);

 filename = data_file.good_data(i).channel_data.filename;
 filename= erase(filename, "cinterocdrft_stab_fft_sigmat.mat") ;
 
for temp_c = C_DE(:,i)'
    for temp_f = 0
   %tridx = STIM.contrast >= 0.5 & STIM.fixedc == 0;
   %data_header.DE0_NDE0 = data_file.good_data(i).channel_data.contrast==0 & data_file.good_data(i).channel_data.fixedc==0;
        if ~isnan(temp_c) && ~isnan(temp_f)
           
           eval(['data_header.DE' num2str(temp_c*100000) '_NDE' num2str(temp_f*100000) '_' filename  '=' ...
               'data_file.good_data(i).channel_data.contrast == ' num2str(temp_c) ...
            ' & data_file.good_data(i).channel_data.fixedc == ' num2str(temp_f)])
        end
    end
end

end


 %select logical vectors of contrasts of interest
    %store them with the data_header fieldname as a name
    %get the corresponding dataset (using the filename ==fieldname)
    %store the layer class in the data stored for each contrast
    
 layer = {'','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
 
 logical_fieldnames = fieldnames(data_header);
 Key = cell2mat(strfind(logical_fieldnames, '_'));
% cont1 = struct();
cont_data = struct();
for nb = 1:numel(fieldnames(data_header))
   fieldname = char(logical_fieldnames(nb));
   %for a more generalized version, use regexp, but more complexe like regexp(fieldname(1:Key(1)), '([0-9\- '']|\. )+', 'match'))
   dom_contrast = str2double(fieldname(3:Key(nb,1)-1));
   
   data_file_filename = [fieldname(Key(nb,2)+1:end) 'cinterocdrft_stab_fft_sigmat.mat'];
   if (0 <= dom_contrast)&&(dom_contrast<= 33333)
       %isolate logicals (not very important)
   %eval(['cont1.' fieldname '='...
      % 'data_header.' fieldname])
       %get the corresponding dataset
       clear i
       for i = 1:length(channum)
           if ~isempty(strfind(data_file.good_data(i).channel_data.filename, data_file_filename))
               eval(['cont_data.cont1_data.' fieldname '.sdftr_chan' '=' ...
                  'data_file.good_data(i).channel_data.sdftr_chan(:,data_header.' fieldname ');'])
              if ~isempty(layer(i))
               eval(['cont_data.cont1_data.' fieldname '.cell_class' '='...
                  'layer(i);'])
              end
           end
       end
   
   elseif (33333 < dom_contrast)&&(dom_contrast<= 66666)
       clear i
       for i = 1:length(channum)
           if ~isempty(strfind(data_file.good_data(i).channel_data.filename, data_file_filename))
               eval(['cont_data.cont2_data.' fieldname '.sdftr_chan' '=' ...
                  'data_file.good_data(i).channel_data.sdftr_chan(:,data_header.' fieldname ');'])
             if ~isempty(layer(i))
               eval(['cont_data.cont2_data.' fieldname '.cell_class' '='...
                  'layer(i);'])
              end
           end
       end
    elseif (66666 < dom_contrast)&&(dom_contrast<= 100000)
       clear i
       for i = 1:length(channum)
           if ~isempty(strfind(data_file.good_data(i).channel_data.filename, data_file_filename))
               eval(['cont_data.cont3_data.' fieldname '.sdftr_chan' '=' ...
                  'data_file.good_data(i).channel_data.sdftr_chan(:,data_header.' fieldname ');'])
             if ~isempty(layer(i))
              eval(['cont_data.cont3_data.' fieldname '.cell_class' '='...
                  'layer(i);'])
              end
           end
       end
   end
end


norm_data = struct();
org_data = struct();
clear channum layer_idx
for c =1:3

eval(['cont_data_fieldnames = fieldnames(cont_data.cont' num2str(c) '_data)';])   


  log_layer = nan(numel(cont_data_fieldnames),1);
  for nb = 1:numel(cont_data_fieldnames)
      
fieldname = char(cont_data_fieldnames(nb));
cell = char('M');
     if eval(['strcmp(cont_data.cont' num2str(c) '_data.' fieldname '.cell_class' ',' 'cell' ')== 1'])
       log_layer(nb) = 1;
     elseif eval(['strcmp(cont_data.cont' num2str(c) '_data.' fieldname '.cell_class' ',' 'cell' ')==0'])
       log_layer(nb) = 0;
     end
  end

layer_idx = find(log_layer);
contrast = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
xabs = -199:1300;
nyq = 15000;
channum = 1:numel(cont_data_fieldnames);
raw_mean_bs = nan(length(xabs), length(channum));
filtered_dSUA = nan(length(xabs), length(channum));
all_locsdSUA_filtered = nan(length(channum),1);
up_dist = nan(1, length(channum));

clear i ;
 for i = 1:length(channum)
  fieldname = char(cont_data_fieldnames(i));  
     if log_layer(i) == 1
   eval(['mean_data = mean(squeeze(cont_data.cont' num2str(c) '_data.' fieldname '.sdftr_chan(401:1900,:,:)),2);']);
   bsl = mean(mean_data(1:200));
   raw_mean_bs(:,i) = mean_data(1:end)- bsl;
   
   lpc       = 100; %low pass cutoff
   lWn       = lpc/nyq;
   [bwb,bwa] = butter(4,lWn,'low');
   %cdata = ;
   lpdSUA      = filtfilt(bwb,bwa, raw_mean_bs(:,i));
  
   %{
   locssu = findpeaks(lpdSUA(50:1201));
x = xabs - locssu.loc(1);
plot(x, lpdSUA)
hold on
%}
    
 filtered_dSUA(:,i) = lpdSUA;
 %find peaks for every channel
  for len = 30:550
            if filtered_dSUA(200+len,i) < filtered_dSUA(200+len+1,i)
   locsdSUA_filtered = findpeaks(filtered_dSUA(200+len:1350,i));
       break
            end
           %break  
  end
         

 %store first peak location
 all_locsdSUA_filtered(i) = locsdSUA_filtered.loc(1)+200+len;
 
  %compute the distance between the first peak and the last datapoint and store  
 %in a matrix
 up_dist(:,i)= length(xabs)- all_locsdSUA_filtered(i);
     end
 end
 
 %get the max distance between the first peak and the stimulus onset
 max_low_dist = max(all_locsdSUA_filtered);
 %create new matrix with the length(max(d)+max(xabs - d))
 new_dist = max_low_dist + max(up_dist);
 clear fp_locked_data norm_fp_locked
 fp_locked_data = nan(new_dist,length(layer_idx));
 norm_fp_locked = nan(new_dist,length(layer_idx));
 for n = 1:length(layer_idx)
    
 lower_bound =max_low_dist-all_locsdSUA_filtered(layer_idx(n))+1;
 upper_bound =max_low_dist-all_locsdSUA_filtered(layer_idx(n))+length(xabs);
 
 fp_locked_data(lower_bound:upper_bound,n) = filtered_dSUA(:,layer_idx(n));
  %{
 for len = 30:1500
            if fp_locked_data(len,n) < fp_locked_data(len+1,n)
locsdSUA_mean = findpeaks(fp_locked_data(len:end, n));
        break
            end 
   end
 x1 = -locsdSUA_mean.loc(1)-len+1:length(fp_locked_data(:,n))-locsdSUA_mean.loc(1)-len;

 %x = 1:length(fp_locked_data(:,1));
 plot(x1, fp_locked_data(:,n))
 hold on
 %}
  norm_fp_locked(lower_bound:upper_bound,n) = (fp_locked_data(lower_bound:upper_bound,n)-min(fp_locked_data(lower_bound:upper_bound,n)))/(max(fp_locked_data(lower_bound:upper_bound,n))-min(fp_locked_data(lower_bound:upper_bound,n)));
 end
 eval(['norm_data.norm_fp_locked' num2str(c) '= norm_fp_locked;'])
 eval(['org_data.fp_locked_data' num2str(c) '= fp_locked_data;'])
 end

h = figure();
clear c
for c= 1:3
    %get rid out of the nans to find out the first peak location
 eval(['mean_filtered = mean(org_data.fp_locked_data' num2str(c) ', 2)']);
 %avoiding the wrong peak is not really important for the mean, but let's
 %do it to make sure it's all aligned well
  for len = 30:1500
            if mean_filtered(len) < mean_filtered(len+1)
 locsdSUA_mean = findpeaks(mean_filtered(len:end));
  break
            end 
   end
 x1 = -locsdSUA_mean.loc(1)-len+1:length(mean_filtered)-locsdSUA_mean.loc(1)-len;
 %include nans back in for the plotting
eval(['mean_with_nans = nanmean(org_data.fp_locked_data' num2str(c) ', 2)']);
 plot(x1,mean_with_nans,'LineWidth',1 )
 hold on
 %ylim([-.1 1.1])
 %green= [167/255 185/255 54/255])
 %black = [24/255 23/255 23/255] )
 %pink = [229/255, 49/255, 90/255])

  hold on
  plot([0 0], ylim,'k','HandleVisibility','off')

      set(gca, 'linewidth',2)
      set(gca,'box','off')
   title({'M class cells mean spiking activity'}, 'Interpreter', 'none')
    xlabel('Time from first{\bf peak} (ms)')
    ylabel('Spike rate (Normalized)')
     xlim([-300 1100])
     
end
legend( 'Mean contrasts1 (0-0.33)','Mean contrasts2 (0.33-0.66)','Mean contrasts3 (0.66-1)','Location', 'bestoutside')
    
  filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\plots\', strcat(contrast{2},'M_class_fpaligned_nanmean_contrasts_xlim'));  
  %saveas(gcf, strcat(filename, '.svg'));
  saveas(gcf, strcat(filename, '.png'));
  
  
  %% only plot contrasts 1
  
  h = figure();

for i = 1:length(org_data.fp_locked_data1(1,:))
 mean_filtered = org_data.fp_locked_data1(:,i);
 for len = 1:1500
            if mean_filtered(len) < filtered_dSUA(len+1)
 locsdSUA_mean = findpeaks(mean_filtered(1:end));
      
            end
              break
 end
 x1 = -locsdSUA_mean.loc(1)-len+1:length(mean_filtered)-locsdSUA_mean.loc(1)-len;
 plot(x1,mean_filtered) %'LineWidth',1 )
 hold on
 %ylim([-.1 1.1])
 %green= [167/255 185/255 54/255])
 %black = [24/255 23/255 23/255] )
 %pink = [229/255, 49/255, 90/255])

  hold on
  plot([0 0], ylim,'k')

      set(gca, 'linewidth',2)
      set(gca,'box','off')
   title({'M class cells mean spiking activity'}, 'Interpreter', 'none')
    xlabel('Time from first{\bf peak} (ms)')
    ylabel('Spike rate (Normalized)')
      xlim([-300 1100])   
end
filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\plots\', strcat(contrast{2},'mean_M_class_fpaligned_norm_mean_contrast1'));  
  %saveas(gcf, strcat(filename, '.svg'));
  saveas(gcf, strcat(filename, '.png'));
  