function [outputplot] = meanlimeplotLoic(mean_data, filetitle)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

sortedLabels = 1:length(mean_data(1,:));
xabs = -50:300;
h = figure();
[ha, pos] = tight_subplot(24,1,[0.001 .03],[.05 .05],[.1 .01]);
for i = 1:length(mean_data(1,:))
    
    axes(ha(i));
   
    plot(xabs,mean_data(1:end,i))
    Aire = area(xabs, mean_data(1:end,i));
    Aire.FaceColor =[0.85 0.85 0.85];
    set(h,'position',get(h,'position').*[1 1 1 1.25]);
    ylh = ylabel(sortedLabels(i));
    ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1)*0.09);
    %xlim([-50 300]);
    ax = gca;
    ax.YAxisLocation = 'origin';
    grid off
    
    if i < length(mean_data(1,:))
    set(ha(i), 'XTick', [])
    ax.XAxis.Visible = 'off';
    end
    if i <= length(mean_data(1,:))
        set(ha(i), 'YTick', [])
    end
    if i == 1
    title(sprintf('LGN %s normalized ', char(filetitle), inputname(1)), 'fontsize', 18)
    end
     
   if i == length(mean_data(1,:))
       tit = title(gca, 'Voltage (�V)','interpreter','none','fontsize',15, 'FontWeight', 'Normal');
       tit.Position(1) = tit.Position(1) - abs(tit.Position(1)*1.6);
       tit.Position(2) = tit.Position(2) + abs(tit.Position(2))*12;
       tit.Rotation =90;
   end
    
end
 set(ha(1:23), 'XTick', [])
 set(ha(1:24), 'box', 'off');
 
 xlh = xlabel('\fontsize{18}time (ms)');
 xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2)*0.7);

outputplot = h;
end

