% Version 1.3: Last modified on 5/15/18
function plotResPages(series, targetRes, loadCapacities, stdDevPercentCapacity, loadlist)

n_data = size(targetRes,1);
n_dim = size(targetRes,2);

targetRes;

targetResPct = 100 * targetRes ./ (ones(n_data,1)*loadCapacities);

thr = 0.25; %residual threshold

srs = unique(series);
n_series = length(srs);
ind_s = zeros(n_series,2);
for i = 1:n_series
    ind_s(i,1) = min(find(series==srs(i)));
    ind_s(i,2) = max(find(series==srs(i)));
end

sub = 0; 
r = min(n_dim,6);
for i = 1:n_dim
    if i>1 && rem(i,6) == 1
        h1=gcf;
        figure('Name',h1.Name,'NumberTitle','off','WindowState','maximized');
        sub = 6;
%         r = n_dim - 6;
        r = 6;
    end
    subplot(r, 1, i-sub); hold on
    axis([1 n_data -1 1])
    
    for j = 1:n_series
        if mod(j,2) == 0
            patch([ind_s(j,1)-0.5 ind_s(j,2)+0.5 ind_s(j,2)+0.5 ind_s(j,1)-0.5], [-1 -1 1 1], [0.8 0.8 0.8])
        end
    end

    plot([1:n_data],targetResPct(:,i))
    plot([1:n_data],thr*ones(1,n_data),'k--',[1:n_data],-thr*ones(1,n_data),'k--')
    xlabel('Data point index');
    ylabel(strcat('\Delta',loadlist{i}));
    title(sprintf('Residual; %% of Load Capacity; Standard Deviation = %0.4f%%',stdDevPercentCapacity(i)));
%     set(gcf,'Position',[100 50 825 min(r*175.5,900)]);
    hold off
end

end