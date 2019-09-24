function f_ShadedLinePlotbyDepth(ASAbyDepth,corticaldepth,tm,N,offset,fill_deflection)

if nargin < 6
    fill_deflection = -1; 
end

if size(ASAbyDepth,1) > size(ASAbyDepth,2)
    ASAbyDepth = ASAbyDepth';
end


%    scale means to values between 1 and -1 with basline at 0
%baslined_ASA  = bsxfun(@minus,ASAbyDepth, mean(ASAbyDepth(:,tm<0),2));
baslined_ASA = ASAbyDepth;
scalingfactor = max(max(abs(baslined_ASA),[],2)); 
scalingrange  = [-1 1];
scaled_ASA    = baslined_ASA .* 1/scalingfactor;

supra_color   = [248 202 161]./255; 
gran_color    = [232 182 212]./255; 
infra_color   = [170 200 211]./255; 

% plot!
for depth = size(ASAbyDepth,1):-1:1
    
    x  = tm;
    bl = -offset*depth ;
    y  = scaled_ASA(depth,:) + bl;
    
    % basline
    baseline = repmat(bl,1,length(x));
    plot(x,baseline,'Color',[.8 .8 .8])
 
    
    % fill positive and negative curves
    sink = y;
    sink(sink>baseline) = bl;
    source = y;
    source(source<baseline) = bl;

    if fill_deflection == 0 
        
     Hl=plot(x,y,'w'); hold on; 

        
    elseif fill_deflection < 0 
    mmfill(x,baseline,source,[0.95 0.95 0.95]), hold on
    
    if ~isempty(corticaldepth)
        if corticaldepth(depth) >= 0.6
            mmfill(x,baseline,sink,supra_color),hold on
        elseif corticaldepth(depth) <= 0.5 &&  corticaldepth(depth) >= 0.0
            mmfill(x,baseline,sink,gran_color),hold on
        else
            mmfill(x,baseline,sink,infra_color),hold on
        end
    else
        mmfill(x,baseline,sink,[0 0 0]), hold on
    end
    
    else
           mmfill(x,baseline,sink,[0.95 0.95 0.95]), hold on
    if ~isempty(corticaldepth)
        if corticaldepth(depth) >= 0.6
            mmfill(x,baseline,source,supra_color),hold on
        elseif corticaldepth(depth) <= 0.5 &&  corticaldepth(depth) >= 0.0
            mmfill(x,baseline,source,gran_color),hold on
        else
            mmfill(x,baseline,source,infra_color),hold on
        end
    else
        mmfill(x,baseline,source,[0 0 0]), hold on
    end 
    end
   
    
    % timecourse
    if ~isempty(corticaldepth)
        if corticaldepth(depth) >= 0.6
            plot(x,y,'Color',supra_color,'LineWidth',1); hold on;
        elseif corticaldepth(depth) <= 0.5 &&  corticaldepth(depth) >= 0.0
            plot(x,y,'Color',gran_color,'LineWidth',1); hold on;
        else
            plot(x,y,'Color',infra_color,'LineWidth',1); hold on;
        end
        hold on;
    else
        
        plot(x,y,'Color',[.8 .8 .8],'LineWidth',1); hold on;
        hold on; 
        
    end

    % label depth and n
    if ~isempty(N)
        if rem(N(depth),5) == 0
            ylabelh = text(min(x),bl,num2str(N(depth)),'HorizontalAlignment','right','FontName', 'Arial','FontSize', 8);
            pos  = get(ylabelh , 'position' );
            pos1 = pos - [10 0 0];
            set(ylabelh , 'position' , pos1 );        
        end
    end
    %{
    if ~isempty(N)
        text(max(x),bl,num2str(N(depth)),'HorizontalAlignment','left','FontName', 'Arial','FontSize', 12); hold on;
    end
%}
end

%    add scale bar and label
ystr = sprintf('Scale Bar = %.2f %s\n', scalingfactor/2*offset);
xlabel(ystr,'fontsize',10); 

% make pretty
xlim([min(tm) max(tm)]);
ylim([-1*offset*size(ASAbyDepth,1) -1*offset] + [-1*offset offset]);
plot([0 0], ylim,'Color',[0 148 68]./255,'linewidth',.8);
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [0 0 0], ...
    'YColor'      , [0 0 0], ...
    'YTick'       , [], ...
    'LineWidth'   , 1         );

% axies labels and title
hYLabel = ylabel(ystr);
set( hYLabel                ,...
    'FontName'   , 'Arial' ,...
    'FontSize'   , 8           );
%{
titleHandle = get( gca ,'ylabel' );
pos  = get( titleHandle , 'position' );
pos1 = pos - [0 0.2 0];
set( titleHandle , 'position' , pos1 );
    %}
%hXLabel= xlabel('Time');
% set( hXLabel                ,...
%     'FontName'   , 'Arial' ,...
%     'FontSize'   , 12      );
end

