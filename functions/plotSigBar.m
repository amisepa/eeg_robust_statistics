%% Plots bars of statistical significance (h) at the botom of the plot
% 
% h must be a logical vector, with true as significant and false as non-significant.
% f = xAxis values (e.g., frequencies)
% 
% h and f must have the same size
% 
% Usage: plotSigBar(h, f)
% 
% Cedric Cannard, 2021

function plotSigBar(h, f)

%Find segments with significant statistical values
seg = [];
count = 1;
for i = 1:length(h)
    if i == 1 && h(i) == 0 
        seg(count,1) = 0;
    elseif h(i) == 1 && i == 1
        seg(count,1) = 1;
    end
    if i>1
        if h(i-1) == 0 && h(i) == 1
        seg(count,1) = i;
        end
    end
    if h(i) == 1 && i == length(h)
        seg(count,2) = i; break
    end
    if h(i) == 1 && h(i+1) == 0
%         if h(i) == 1 && i == 1 && h(i+1) == 0 %if only first time point is significant
%             seg(count,2) = i+0.5;
%         else    %else, end of the significant segment
            seg(count,2) = i;
%         end
        count = count+1;
    end
end

if sum(h) ~= 0
    seg = f(seg);
end

%plot bars
if size(seg,2) == 2
    xSize = diff(f(1:2));
    for i = 1:size(seg,1)
%         hold on; plot([seg(i,1)-xSize/10 seg(i,2)+xSize/10],[min(ylim)+0.00001 min(ylim)+0.00001],'k','LineWidth',6);
        % hold on; plot([seg(i,1)-xSize/2 seg(i,2)+xSize/2],[min(ylim)+0.001 min(ylim)+0.001],'k','LineWidth',6);
        hold on; plot([seg(i,1)-xSize seg(i,2)+xSize], [min(ylim)+0.00001 min(ylim)+0.00001],'k','LineWidth',6);
    end
end
