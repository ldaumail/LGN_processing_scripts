function [outputplot] = LimePlot(data, filetitle, filename, contrast)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%baseline adjusted data process

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
if any(strfind(inputname(1), 'LFP')) || any(strfind(inputname(1), 'MUA'))
norm_mean_adj = norm_mean(:,:) - mean_bp;
else
 norm_mean_adj = -(norm_mean(:,:) - mean_bp); 
end

sortedLabels = 1:length(data(1,:,1));
xabs = -50:1500;
h = figure();
[ha, pos] = tight_subplot(24,1,[0.001 .03],[.05 .05],[.1 .01]);
for i = 1:length(norm_mean_adj(1,:))
    
    axes(ha(i));
   
    plot(xabs,norm_mean_adj(:,i))
    Aire = area(xabs, norm_mean_adj(:,i));
    Aire.FaceColor =[0.85 0.85 0.85];
    set(h,'position',get(h,'position').*[1 1 1 1.25]);
    ylh = ylabel(sortedLabels(i));
    ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1)*0.09);
    %xlim([-50 300]);
    ax = gca;
    ax.YAxisLocation = 'origin';
    grid off
    
    if i < length(norm_mean_adj(1,:))
    set(ha(i), 'XTick', [])
    ax.XAxis.Visible = 'off';
    end
    if i <= length(norm_mean_adj(1,:))
        set(ha(i), 'YTick', [])
    end
    if i == 1
   if any(strfind(inputname(1), 'LFP')) || any(strfind(inputname(1), 'MUA'))
 title({sprintf('LGN %s ',  char(filetitle)), 'baseline-corrected, averaged, then normalized', inputname(1) contrast}, 'Interpreter', 'none', 'fontsize', 12)
 clrbar.Label.String = 'no unit';
else
 title({sprintf('LGN %s ',  char(filetitle)), 'baseline-corrected, averaged', inputname(1) contrast}, 'Interpreter', 'none', 'fontsize', 12)
 clrbar.Label.String = 'nA/mm^3';
   end
    end
     
   if i == length(norm_mean_adj(1,:))
       tit = title(gca, 'Normalized (no unit)','interpreter','none','fontsize',15, 'FontWeight', 'Normal');
       tit.Position(1) = tit.Position(1) - abs(tit.Position(1)*1.6);
       tit.Position(2) = tit.Position(2) + abs(tit.Position(2))*12;
       tit.Rotation =90;
   end
    
end
 set(ha(1:23), 'XTick', [])
 set(ha(1:24), 'box', 'off');
 
 xlh = xlabel('\fontsize{18}time (ms)');
 xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2)*0.7);
 if any(strfind(inputname(1), 'LFP')) 
saveas(gcf,strcat(filename, 'baselinecorr_mean_normalized_LFP_limeplot',contrast,'.png'));
 else
     if any(strfind(inputname(1), 'aMUA'))
  saveas(gcf,strcat(filename, 'baselinecorr_mean_normalized_aMUA_limeplot',contrast,'.png')); 
     else
         if any(strfind(inputname(1), 'dMUA'))
  saveas(gcf,strcat(filename, 'baselinecorr_mean_normalized_dMUA_limeplot',contrast,'.png')); 
         else
             if any(strfind(inputname(1), 'CSD'))
  saveas(gcf,strcat(filename, 'baselinecorr_mean_CSD_limeplot',contrast,'.png'));        
             end
         end
     end
 end
outputplot = h;
end

