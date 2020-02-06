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

%% determine overall baseline of a given unit across all trials and contrasts


%{
all_bsls = struct();
 for c =1:3
     eval(['cont_data_fieldnames = fieldnames(cont_data.cont' num2str(c) '_data);'])
     channum = 1:numel(cont_data_fieldnames);
     for i = channum
         fieldname1 = char(cont_data_fieldnames(i));
         underscores = strfind(fieldname1, '_');
         su_name1 = char(fieldname1(underscores(2)+1:end-1));
         eval(['all_bsls.all_bsls' su_name1 '_' num2str(c) '= [];']);
       for nb = channum
          fieldname2 = char(cont_data_fieldnames(nb));
          su_name2 = char(fieldname2(underscores(2)+1:end-1));
          
          if strcmp(su_name1, su_name2) == 1
              cell = char('M');
              if eval(['strcmp(cont_data.cont' num2str(c) '_data.' fieldname2 '.cell_class' ',' 'cell' ')== 1'])
                  eval(['trialidx = 1:length(cont_data.cont' num2str(c) '_data.' fieldname2 '.sdftr_chan(1,:));']);
                  eval('bsl = nan(length(trialidx),1);')
                  for tridx = trialidx
                  eval(['data = cont_data.cont' num2str(c) '_data.' fieldname2 '.sdftr_chan(401:1900,tridx);']);
                  eval('bsl(tridx,1) = mean(data(1:200));');
                  end
          eval(['all_bsls_mat = all_bsls.all_bsls' su_name1 '_' num2str(c) ';']) 
          all_bsls_mat = [all_bsls_mat; bsl];
          eval(['all_bsls.all_bsls' su_name1 '_' num2str(c) '= all_bsls_mat;'])
          eval(['clear bsl' su_name2 num2str(nb)])
              end
               
          end
         
       end

      end
 end
 
all_bsls_fieldnames = fieldnames(all_bsls);
all_bsls_stack = struct();
 for c =1:3
     channum2 = 1:numel(all_bsls_fieldnames);
     for i = channum2
         fieldname3 = char(all_bsls_fieldnames(i));
         underscores = strfind(fieldname3, '_');
         su_name3 = char(fieldname3(1:underscores(5)-1));
       for nb = channum2
         fieldname4 = char(all_bsls_fieldnames(i));
         underscores = strfind(fieldname4, '_');
        
         su_name4 = char(fieldname4(1:underscores(5)-1));
            if strcmp(su_name3, su_name4) == 1
                testfield =char(strcat(su_name3, '_', num2str(c)));
                if eval(['isfield(all_bsls,' 'testfield' ') ==1;']) 
          eval(['bsls = all_bsls.' su_name3 '_' num2str(c) ';']) 
          all_bsls_mat = [all_bsls_mat; bsls];
          eval(['all_bsls_stack.' su_name3 '= all_bsls_mat;'])
                end
            end
       end
     end
 
 end
 %}
channum = 1: length(data_file.good_data); 

all_bsls = struct();
for i = channum
    filefieldname = data_file.good_data(i).channel_data.filename;
    Key = strfind(filefieldname, '_');
    filefieldname = filefieldname(1:Key(4)-1);   
  eval(['all_bsls.dat_' filefieldname '_tad.bsls =mean(data_file.good_data(i).channel_data.sdftr_chan(401:600,:),1);'])
  eval(['all_bsls.dat_' filefieldname '_tad.std_bsls = std(all_bsls.dat_' filefieldname '_tad.bsls);'])
  eval(['all_bsls.dat_' filefieldname '_tad.mean_bsls = mean(all_bsls.dat_' filefieldname '_tad.bsls);'])
  
  eval(['all_bsls.dat_' filefieldname '_tad.bslcorr_data = data_file.good_data(i).channel_data.sdftr_chan(401:600,:) - all_bsls.dat_' filefieldname '_tad.bsls ;'])
  eval(['all_bsls.dat_' filefieldname '_tad.bsls_bslcorr = mean(all_bsls.dat_' filefieldname '_tad.bslcorr_data(1:200,:),1);'])
  eval(['all_bsls.dat_' filefieldname '_tad.std_bsls_bslcorr = std(all_bsls.dat_' filefieldname '_tad.bsls_bslcorr);'])
  eval(['all_bsls.dat_' filefieldname '_tad.mean_bsls_bslcorr = mean(all_bsls.dat_' filefieldname '_tad.bsls_bslcorr);'])
  
end        
%% look at the power at 4Hz for the non responses
%{

[power(tr,:), freq(tr,:)] = calcFFT(squeeze(STIMdMUA.STIM.sdftr(600:1700,tr)));

%find the index of the frequency vector closest to 4hz and point to the
%power value of this index for every trial, and store the value in
%fourhzpower
[val,index] = min(abs(4-freq(tr,:)));
fourhzpower(tr) = power(tr,index);
%}

%%

clear channum layer_idx
bsl_fieldnames = fieldnames(all_bsls);


%contrast = {'DE0_NDE50','DE50_NDE0','DE50_NDE50'};
xabs = -199:1300;
%xabs2 =-599:1300;
nyq = 15000;
c =1;
eval(['cont_data_fieldnames = fieldnames(cont_data.cont' num2str(c) '_data);']) 
channum = 1:numel(cont_data_fieldnames);
up_dist = nan(1, length(channum));
all_locsdSUA_filtered = nan(1,length(channum));
mean_filtered_dSUA = nan(length(xabs), length(channum));
signi = nan(length(channum),1);
mean_data_file_bsl = nan(length(channum),1);
mean_bsl = nan(length(channum),1);

for i = channum

fieldname = char(cont_data_fieldnames(i)); 
eval(['trialidx = 1:length(cont_data.cont' num2str(c) '_data.' fieldname '.sdftr_chan(1,:));']);

raw_bs = nan(length(xabs), length(trialidx));
filtered_dSUA = nan(length(xabs), length(trialidx));

mean_wnd1 = nan(1,length(trialidx));
%max_wnd1 = nan(1,length(trialidx));

bsl = nan(1, length(trialidx));
   for tridx = trialidx
    
         
    eval(['data = cont_data.cont' num2str(c) '_data.' fieldname '.sdftr_chan(401:1900,tridx);']);
    bsl(tridx) = mean(data(1:200));
    raw_bs(:,tridx) = data(1:end)- bsl(tridx);
   
    lpc       = 100; %low pass cutoff
    lWn       = lpc/nyq;
    [bwb,bwa] = butter(4,lWn,'low');
    lpdSUA      = filtfilt(bwb,bwa, raw_bs(:,tridx));
  
   %{
  %locssu = findpeaks(lpdSUA(50:1201));
%x = xabs - locssu.loc(1);
%x = 1:length(lpdSUA);
plot(xabs2, lpdSUA)
hold on
%}

    filtered_dSUA(:,tridx) = lpdSUA;
    mean_wnd1(tridx) = mean(lpdSUA(231:480));
   end
    %{
plot(xabs2, mean(filtered_dSUA,2),'LineWidth',1, 'Color', 'black')  
hold on
    %}
   mean_filtered_dSUA(:,i) = mean(filtered_dSUA,2);
   underscores = strfind(fieldname, '_');
   splitfieldname = fieldname(underscores(2)+1:end-1);
    
   for nb = 1:numel(bsl_fieldnames)
      bsl_fieldname = char(bsl_fieldnames(nb));
      underscores2 = strfind(bsl_fieldname, '_');
      splitfieldname2 = bsl_fieldname(underscores2(1)+1:underscores2(5)-1);

       if strcmp(splitfieldname, splitfieldname2) ==1
         eval(['data_file_bsl = all_bsls.dat_' splitfieldname2 '_tad.bsls_bslcorr;'])
         %eval(['data_file_std_bsl = all_bsls.dat_' splitfieldname2 '_tad.std_bsls_bslcorr;'])
       
       end
   end 

   [signi(i,1), pvalue] = ttest2(data_file_bsl, mean_wnd1);
   mean_data_file_bsl(i,1) = mean(data_file_bsl);
   mean_bsl(i,1) = mean(mean_wnd1);
   if signi(i,1) == 1 &&  mean(data_file_bsl) < mean(mean_wnd1)
 %find peaks for every unit
    for len = 30:550
        
        if mean_filtered_dSUA(200+len,i) < mean_filtered_dSUA(200+len+1,i)
              
            
             locsdSUA_filtered = findpeaks(mean_filtered_dSUA(200+len:1350,i));
             %store first peak location 
             all_locsdSUA_filtered(:,i) = locsdSUA_filtered.loc(1)+200+len;
             %locspk1(tridx) =locsdSUA_filtered.loc(1)+200+len;
            % locspk2(tridx) =locsdSUA_filtered.loc(2)+200+len;
             break 
        end 
    end
   else
     %if there is no peak in the first window where first peak
     %is expected, generate 4 arbitrary linearly spaced peaks locations 
     locsdSUA_filtered = linspace(200,1350,4); 
     all_locsdSUA_filtered(:,i) = locsdSUA_filtered(2);
                     
  end

     
  %{  
     if filtered_dSUA(all_locsdSUA_filtered(:,tridx),tridx) < data_file_bsl+ 1.96*data_file_std_bsl
              
         locsdSUA_filtered = findpeaks(filtered_dSUA(locspk2(tridx):1350,tridx));
         all_locsdSUA_filtered(:,tridx) = locsdSUA_filtered.loc(1)+locspk2(tridx)-1;
             
     end
    %}
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
 if signi(n,1) ==1 && mean_data_file_bsl(n,1) < mean_bsl(n,1)
    
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
  % 
         else
       fp_locked_data(:,n) = nan(length(fp_locked_data(:,n)),1);
       norm_fp_locked(:,n) = nan(length(fp_locked_data(:,n)),1);
       
         end
        %}
  
 end
 
fp_locked_data_out = fp_locked_data(:,~all(isnan(fp_locked_data))); % for nan - cols
 
underscore = strfind(fieldname, '_');
domcont = str2double(fieldname(3:underscore(1)-1))/100000;
figure();
x = 1:length(fp_locked_data_out(:,1));
 plot(x, fp_locked_data_out(:,:),'HandleVisibility','off')
 hold on
 plot(x, mean(fp_locked_data_out(:,:),2),'LineWidth',1, 'Color', 'black')
 plot(xlim, [0 0],'k','HandleVisibility','off')
 plot(xlim, [mean(data_file_bsl) mean(data_file_bsl)],'k')
 %plot(xlim, [mean(data_file_bsl)+1.96*data_file_std_bsl mean(data_file_bsl)+1.96*data_file_std_bsl],'k')
 
  txtmeanBSL = ['mean baseline'];
  text(0, mean(data_file_bsl)+3, txtmeanBSL)
  %txt95ci = ['mean baseline +1.96std(bsl)'];
  %text(0, data_file_bsl+1.96*data_file_std_bsl+3, txt95ci)
  title(strcat({'Contrast =' domcont}), 'Interpreter', 'none')
  
 