function [taretal,tarestdevsub] = meantare(series,checkit)

% SOLVE TARES BY TAKING THE MEAN
[s,~,s_id] = unique(series);
nseries = length(s);
tares=zeros(nseries,size(checkit,2));
tarestd=zeros(nseries,size(checkit,2));

for i = 1:nseries
    series_Res = checkit(s_id==i,:);
    tares(i,:) = mean(series_Res,1);
    tarestd(i,:) = std(series_Res,0,1);
end
taretal = tares(s_id,:);
tarestdevsub = tarestd(s_id,:); 

