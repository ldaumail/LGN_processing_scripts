function [outputplot] = newLimePlot(data, BRdatafile, filename, contrast)
%newLimePlot plots all the data (LFP, aMUA, CSD) stored in the 4-D vector
%'data', for which the dimensions are as follows:
%dim1: time points; dim2: channel indices, dim3: trials, dim4: data type
%(either LFP, aMUA, or CSD)
%   Detailed explanation goes here
%baseline adjusted data process

h = figure();
for inp = 1:3
%compute mean of the data across trials
mean_data = mean(squeeze(data(:,:,:,inp)),3);
% if the data is either LFP or MUA, do a normalization (data - min)/(max-min)
if any(strfind(inputname(1), 'LFP')) || any(strfind(inputname(1), 'MUA'))
 norm_mean = nan(length(mean_data(:,1)),length(mean_data(1,:)));
     for n=1:24
     norm_mean(:,n) = (mean_data(:,n) - min(mean_data(:,n)))/(max(mean_data(:,n))-min(mean_data(:,n)));
     end 
else 
    %no normalization for CSD
   norm_mean = mean_data; 
end
%baseline corrected data process
basedata = norm_mean(25:75,:);
mean_bp = mean(basedata,1);

if inp == 1 || inp == 2
norm_mean_bscorr = norm_mean(:,:) - mean_bp;
else
    %Make sure CSD gets inverted
 norm_mean_bscorr = -(norm_mean(:,:) - mean_bp); 
end
xabs = -50:1500;

subplot(1,4,inp)
%f_ShadedLinePlotbyDepth(norm_mean_bscorr, 0:(1/length(norm_mean_bscorr(1,:))):1,-50:1500,1:length(norm_mean_bscorr(1,:)),1);
f_ShadedLinePlotbyDepth(norm_mean_bscorr, 1:length(norm_mean_bscorr(1,:)),xabs,1:length(norm_mean_bscorr(1,:)),1);
hold on
plot([0 0], ylim,'k')
plot([1100 1100], ylim,'k')
sgtitle({strcat(sprintf('Scaled LGN %s averaged data',  char(BRdatafile)))}, 'Interpreter', 'none', 'fontsize', 30, 'FontWeight', 'bold')

xlabel('\fontsize{9}time (ms)')
if inp == 1
title({strcat('LFP ',contrast), 'normalized, and baseline-corrected'}, 'Interpreter', 'none', 'fontsize', 30)
ylabel({'\fontsize{9}Contacts indexed down from surface','\fontsize{9}(no unit)'})
else
if inp == 2
  title({strcat('aMUA ',contrast), 'normalized, and baseline-corrected'}, 'Interpreter', 'none', 'fontsize', 30)
ylabel({'\fontsize{9}Contacts indexed down from surface','\fontsize{9}(no unit)'})
else
title({strcat('CSD ', contrast), 'baseline-corrected'}, 'Interpreter', 'none', 'fontsize', 30)
ylabel({'\fontsize{9}Contacts indexed down from surface','\fontsize{9}(nA/mm^3)'})
end
end
 if inp ==3
     subplot(1,4,4)
     fint_data = filterCSD(norm_mean_bscorr')';
     scalingfactor = max(max(abs(fint_data),[],1)); 
     scalingrange  = [-1 1];
     scaled_data    = fint_data .* 1/scalingfactor;
     yvec = [1:24];
     imagesc(xabs,yvec, fint_data')
     colormap(jet);
     colorbar;
     hold on;
     plot([1100 1100], ylim,'k')
     title({strcat('CSD ', contrast), 'baseline-corrected'}, 'Interpreter', 'none', 'fontsize', 30)
     xlabel('\fontsize{9}time (ms)')
     ylabel({'\fontsize{9}Contacts indexed down from surface','\fontsize{9}(nA/mm^3)'})
 end
hold off
end
x_width=40; y_width=30;
set(h, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf,strcat(filename, 'all',contrast, '_plots', '.png'));

