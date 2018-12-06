function [localZeros,localZerosAllPoints] = localzeros(series,excessVec0)
numpts = size(series,1);
[~,s_1st,s_id] = unique(series);

localZeros = excessVec0(s_1st,:);
localZerosAllPoints = localZeros(s_id,:);