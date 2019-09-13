function [outputplot] = SurfacePlot(mean_data, filetitle)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% if the data is either LFP or MUA, do a normalization (data -
% min)/(max-min)
if any(strfind(inputname(1), 'LFP')) || any(strfind(inputname(1), 'MUA'))
 norm_data = nan(length(mean_data(1,:)),length(mean_data(:,1)));
     for n=1:24
     norm_data(n,:) = (mean_data(:,n) - min(mean_data(:,n)))/(max(mean_data(:,n))-min(mean_data(:,n)));
     end 
   fcsd = filterCSD(norm_data)';
else 
   fcsd = filterCSD(mean_data')'; % filter baseline corrected CSD
end

h = figure();
 xabs = [-50:2000];
 yvec = [1:14,16:24];
 
 imagesc(xabs,yvec, fcsd')
 
 if any(strfind(inputname(1), 'CSD'))
 colormap(flipud(jet));
 end 
 colorbar;
 %climit = max(abs(get(gca,'CLim'))*.8);
 %set(gca,'CLim',[3 13],'Box','off','TickDir','out')
 hold on;
plot([0 0], ylim,'k')
plot([1000 1000], ylim,'k')
clrbar = colorbar;
if any(strfind(inputname(1), 'LFP')) || any(strfind(inputname(1), 'MUA'))
 title({sprintf('LGN %s normalized', char(filetitle)), inputname(1)}, 'Interpreter', 'none', 'fontsize', 18)
 clrbar.Label.String = 'µV';
else
 title({sprintf('LGN %s baseline-corrected', char(filetitle)), inputname(1)}, 'Interpreter', 'none', 'fontsize', 18)
 clrbar.Label.String = 'nA/mm^3';
end
 xlabel('time (ms)')
 ylabel('Contact index')

outputplot = h;
end

