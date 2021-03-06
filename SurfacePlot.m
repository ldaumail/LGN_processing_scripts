function [outputplot] = SurfacePlot(data, filetitle, filename, contrast)
%this function creates surface plots using imagesc function, of either LFP,
%MUA, CSD. During this process, either LFP, aMUA/dMUA,or CSD is baseline adjusted(corrected);
%-baseline corrected LFP mean, baseline corrected MUA mean, or baseline corrected CSD mean is then computed; 
%-LFP or MUA: normalized;
%-filtered and interpolated with filterCSD function and plotted.

%compute mean of the data across trials
mean_data = mean(data,3);
% if the data is either LFP or MUA, do a normalization (data - min)/(max-min)
if any(strfind(inputname(1), 'LFP')) || any(strfind(inputname(1), 'MUA'))
 norm_mean = nan(length(mean_data(:,1)),length(mean_data(1,:)));
     for n=1:24
     norm_mean(:,n) = (mean_data(:,n) - min(mean_data(:,n)))/(max(mean_data(:,n))-min(mean_data(:,n)));
     end 
else 
   norm_mean = mean_data; 
end
%baseline adjusted data process
basedata = norm_mean(25:75,:);
mean_bp = mean(basedata,1);

norm_mean_adj = norm_mean(:,:) - mean_bp;

fint_data = filterCSD(norm_mean_adj')';

h = figure();
 xabs = [-50:1500];
 yvec = [1:14,16:24];
 
 imagesc(xabs,yvec, fint_data')
 
 if any(strfind(inputname(1), 'CSD'))
 colormap(flipud(jet));
 end 
 colorbar;
 %climit = max(abs(get(gca,'CLim'))*.8);
 %{
 if any(strfind(inputname(1), 'CSD'))
 set(gca,'CLim',[-800 800],'Box','off','TickDir','out')
 else
 end
 %}
 hold on;
plot([0 0], ylim,'k')
plot([1100 1100], ylim,'k')
clrbar = colorbar;
if any(strfind(inputname(1), 'LFP')) || any(strfind(inputname(1), 'MUA'))
 title({sprintf('LGN %s ',  char(filetitle)), 'baseline-corrected, averaged,', 'then normalized' inputname(1), contrast}, 'Interpreter', 'none', 'fontsize', 18)
 clrbar.Label.String = 'no unit';
else
 title({sprintf('LGN %s ',  char(filetitle)), 'baseline-corrected, averaged', inputname(1), contrast}, 'Interpreter', 'none', 'fontsize', 18)
 clrbar.Label.String = 'nA/mm^3';
end
 xlabel('time (ms)')
 ylabel('Contact index')

 if any(strfind(inputname(1), 'LFP')) 
saveas(gcf,strcat(filename, 'baselinecorr_mean_normalized_LFP_surfaceplot',contrast,'.png'));
 else
     if any(strfind(inputname(1), 'aMUA'))
  saveas(gcf,strcat(filename, 'baselinecorr_mean_normalized_aMUA_surfaceplot',contrast,'.png')); 
     else
         if any(strfind(inputname(1), 'dMUA'))
  saveas(gcf,strcat(filename, 'baselinecorr_mean_normalized_dMUA_surfaceplot',contrast,'.png')); 
         else
             if any(strfind(inputname(1), 'CSD'))
  saveas(gcf,strcat(filename, 'baselinecorr_mean_CSD_surfaceplot',contrast,'.png'));        
             end
         end
     end
 end
outputplot = h;
end

