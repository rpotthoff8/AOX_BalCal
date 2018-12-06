function taretal = meantare(series,checkit)
% SOLVE TARES BY TAKING THE MEAN
[~,~,s_id] = unique(series);
for i = 1:nseries
    zoop = checkit(s_id==i,:);
    zap(i,:) = mean(zoop);
end
taretal = zap(s_id,:);