newdatadir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\su_peaks_02252020_2\all_units\';
channelfilename = [newdatadir 'clean_SUA_sup_50']; 
data_file = load(channelfilename);

xabs = -199:1300;

nyq = 500;
channum = 1: length(data_file.clean_high_SUA);
mean_filtered_dSUA = struct();
suas_trials = struct();
up_dist = nan(1, length(channum));
all_locsdSUA_filtered = nan(1,length(channum));


for i = channum  
    if ~isempty(data_file.clean_high_SUA(i).namelist)
    trialidx = 1:length(data_file.clean_high_SUA(i).namelist(1,:));
   
   filtered_dSUA = data_file.clean_high_SUA(i).namelist;
  % filtered_dSUA = filtered_dSUA(:,~all(isnan(filtered_dSUA))); % for nan - cols
   
   
   bsl_peaks = nan(1, length(filtered_dSUA(1,:)));
   for tr = trialidx
      for loc = 1:200
          if filtered_dSUA(loc,tr) < filtered_dSUA(loc+1,tr)
      bsl_peak_locs = findpeaks(filtered_dSUA(loc:200,tr));
      bsl_peaks(1,tr) = max(filtered_dSUA(bsl_peak_locs.loc+loc,tr));
      break
          end 
      end
   end
   
   
   %if length(filtered_dSUA(1,:)) >=10
       %store unaligned data
    suas_trials(i).unaligned = filtered_dSUA;
    %determine the first peak location for each trial of a given single
    %unit
    all_locsdSUA_trials = nan(1,length(filtered_dSUA(1,:)));
    peak1 = nan(1, length(filtered_dSUA(1,:)));
  
    for trial = 1:length(filtered_dSUA(1,:))
         for ln = 1:550
             if ~all(isnan(filtered_dSUA(:,trial))) && filtered_dSUA(200+ln,trial) < filtered_dSUA(200+ln+1,trial)
            
                 locsdSUA_trial = findpeaks(filtered_dSUA(200+ln:1350,trial));
             
                if filtered_dSUA(locsdSUA_trial.loc(1)+200+ln,trial) > 0.4*filtered_dSUA(locsdSUA_trial.loc(2)+200+ln)
         %store first peak location 
             all_locsdSUA_trials(:,trial) = locsdSUA_trial.loc(1)+200+ln;
                else
             all_locsdSUA_trials(:,trial) = locsdSUA_trial.loc(2)+200+ln;    
            
                end
                 %get distribution of peak 1 across trials
             peak1(1,trial) = filtered_dSUA(all_locsdSUA_trials(1,trial),trial);
             break 
             end 
         end
         
    end
    up_dist_trials= length(xabs)- all_locsdSUA_trials;
    max_low_dist_unit = max(all_locsdSUA_trials(1,:));
    %create new matrix with the length(max(d)+max(xabs - d))
    new_dist_unit = max_low_dist_unit + max(up_dist_trials); 
    fp_locked_trials = nan(new_dist_unit,length(filtered_dSUA(1,:)));
   
    %find peak1 outliers
    outlier = isoutlier(peak1);
    %find bsl peaks outliers
    out_bsl_peaks = isoutlier(bsl_peaks);
     clear n
     for n = 1:length(filtered_dSUA(1,:))
        if ~isnan(all_locsdSUA_trials(n)) && outlier(n) == 0 %&& out_bsl_peaks(n) == 0
         lower_unit_bound =max_low_dist_unit-all_locsdSUA_trials(n)+1;
         upper_unit_bound =max_low_dist_unit-all_locsdSUA_trials(n)+length(xabs);
     
         fp_locked_trials(lower_unit_bound:upper_unit_bound,n) = filtered_dSUA(:,n);
        end
     end
     fp_locked_trials_out = fp_locked_trials(:,~all(isnan(fp_locked_trials))); % for nan - cols
      %compute the mean single unit activity if more than 10 trials
      mean_filtered_dSUA(i).mean_unit = mean(fp_locked_trials_out,2); % for nan - cols
     %get the aligned data if it exists for the unit 
     suas_trials(i).aligned = fp_locked_trials_out; 
  
  % elseif length(filtered_dSUA(1,:)) <10
  %    mean_filtered_dSUA(i).mean_unit = nan(length(fp_locked_trials(:,1)),1);  
  %    suas_trials(i).unaligned = nan(1,1);
  % end
  end 
   %{
   figure();
   x = 1:length(fp_locked_trials_out(:,1));
   plot(x,fp_locked_trials)
   hold on
   plot(x, mean(fp_locked_trials_out,2),'LineWidth',1, 'Color', 'black')
   
   %mean_aligned = mean(fp_locked_trials_out,2);
   %nanmean_aligned = nanmean(fp_locked_trials_out,2);
   %}
   %{
   %%%%%%% align mean single units to first peak %%%%
      %start finding the peaks at the first non NaN location 
        for len1 =1:800
            if ~isnan(mean_filtered_dSUA(i).mean_unit(len1))
                break
            end
        end
        for len2 = len1+30:1550
            if ~all(isnan(mean_filtered_dSUA(i).mean_unit)) && mean_filtered_dSUA(i).mean_unit(len2) < mean_filtered_dSUA(i).mean_unit(len2+1)

             locsdSUA_filtered = findpeaks(mean_filtered_dSUA(i).mean_unit(len2:end));
             
                if mean_filtered_dSUA(i).mean_unit(locsdSUA_filtered.loc(1)+len2) > 0.6*mean_filtered_dSUA(i).mean_unit(locsdSUA_filtered.loc(2)+len2)
         %store first peak location 
             all_locsdSUA_filtered(:,i) = locsdSUA_filtered.loc(1)+len2;
                else
             all_locsdSUA_filtered(:,i) = locsdSUA_filtered.loc(2)+len2;    
            
                end
             break 
            end 
        end

    %compute the distance between the first peak and the last datapoint and store  
   %in a matrix
   up_dist(:,i)= length(mean_filtered_dSUA(i).mean_unit)- all_locsdSUA_filtered(i);
  %}
end

pvaluesdir = 'C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\lmer_results_peaks\';
 pvalfilename = [pvaluesdir 'lmer_results_rawbs.csv'];
 pvalues = dlmread(pvalfilename, ',', 1,1);
 
%% Plots of peaks with pvalues 
%exclude 160517, (first unit, left empty, it is a K neuron)
%Reject 180806 p1 uclust17, M cell, as doesn't seem well triggered (46)
%Reject 181207 (B) uclust22, M cell, as doesn't seem well triggered (55)
 layer = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
layer([1,46,55]) = []; 
channum = 1: length(data_file.clean_high_SUA);

for chan = 1:15
 
     mean_unaligned = mean(suas_trials(chan).unaligned,2);
     unaligned = suas_trials(chan).unaligned;
     mean_aligned = mean(suas_trials(chan).aligned,2);
     aligned = suas_trials(chan).aligned;
     if ~isempty(mean_unaligned) && ~isempty(mean_aligned)
    h = figure;
    xabs = -199:1300;
    subplot(2, 1,1);
    plot(xabs, unaligned)
    hold on
    plot(xabs, mean_unaligned,'LineWidth',1, 'Color', 'black')
    hold on
    plot([0 0], ylim,'k')
    hold on
    plot([1150 1150], ylim,'k')
 
   
    xalign = 1:length(suas_trials(chan).aligned(:,1));
    subplot(2, 1,2);
    plot(xalign, aligned)
    hold on
    plot(xalign, mean_aligned,'LineWidth',1, 'Color', 'black')
  
    ylh = ylabel({'\fontsize{9}Contacts','\fontsize{9}Spike Rate (spikes/s)'});
   
 
      set(gca, 'linewidth',2)
      set(gca,'box','off')

    sgtitle({sprintf('%s | %s', num2str(chan), char(layer(chan))), 'responses, p<0.05, associated to adaptation pvalues'}, 'Interpreter', 'none')
    xlabel('Time from -50ms from stimulus onset (ms)')
   set(gcf,'Units','inches') 
   set(gcf,'position',[1 1 8.5 11])
   % filename = strcat('C:\Users\maier\Documents\LGN_data\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\plots\',strcat( sprintf('bscorr_data_aligned_unaligned_P1040P2_1ms_stimonset_unit_%d_exclude4peaks', chan)));
    %saveas(gcf, strcat(filename, '.png'));
     end
end
