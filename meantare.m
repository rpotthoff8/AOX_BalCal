function [taretal,tarestdevsub] = meantare(series,checkit)
%Function solves for tare in each series by taking mean of residual

%INPUTS:
%  series = Labels of series number for each point
%  checkit  =  Residual input.  (Non tare-corrected load approximation) -(target load)

%OUTPUTS:
%  taretal = tall matrix with tare value for each datapoint
%  tarestdevsub = tall matrix with standard deviation of series checkit for each datapoint

% SOLVE TARES BY TAKING THE MEAN
[s,~,s_id] = unique(series); %Determine series labels and break points
nseries = length(s); %Count total number of series

%Initialize variables
tares=zeros(nseries,size(checkit,2));
tarestd=zeros(nseries,size(checkit,2));
for i = 1:nseries %Loop through for each series
    series_Res = checkit(s_id==i,:); %Take subset of checkit for specific series
    tares(i,:) = mean(series_Res,1); %Tares defined as mean of checkit (residual) in series
    tarestd(i,:) = std(series_Res,0,1); %Standard deviation of residuals in series
end
taretal = tares(s_id,:); %Expand out to number of datapoints
tarestdevsub = tarestd(s_id,:); %Expand out to number of datapoints

