function [series_lhs,targetMatrix_lhs,excessVec_lhs, sample] = AOX_LHS(series,targetMatrix,excessVec,p)
thr = 0.05;

[n, dim] = size(excessVec);
XI = [1:n]';

[srs, i_s1st, i_s] = unique(series);
n_s = length(srs);

tM_tmp = excessVec - ones(n,1)*min(excessVec);
x_norm = tM_tmp./( ones(n,1)*max(tM_tmp) + eps );

LZ = XI(i_s1st);

XI_tmp  = XI;        i_s_tmp = i_s;
XI_tmp(i_s1st) = []; i_s_tmp(i_s1st) = [];
x_norm_tmp = x_norm(XI_tmp,:);

boundary = (x_norm_tmp<thr) + (x_norm_tmp>(1-thr));
bnd = sum(boundary,2) > 0;
BD = XI_tmp(bnd);

XI_tmp(bnd) = []; i_s_tmp(bnd) = [];
x_norm_tmp = x_norm(XI_tmp,:);

x_tmp = x_norm_tmp - ones(size(x_norm_tmp,1),1)*min(x_norm_tmp);
x_renorm = x_tmp./( ones(size(x_tmp,1),1)*max(x_tmp) + eps );

n_sample = round(p*size(x_renorm,1));
x_lhs = lhsdesign(n_sample,dim);

D = pdist2(x_renorm,x_lhs);
XI_sample = [];
for i = 1:n_sample
    [D_min1,DI] = min(D,[],2);
    [~,DI_2] = min(D_min1);
    DI_1 = DI(DI_2);
    XI_sample = [XI_sample; XI_tmp(DI_2)];
    
    D(DI_2,:) = []; D(:,DI_1) = [];
    XI_tmp(DI_2) = [];
end

X_cal = sort([LZ;BD;XI_sample]);
%X_val = sort([LZ;XI_tmp]);

series_lhs = series(X_cal,:);
excessVec_lhs = excessVec(X_cal,:);
targetMatrix_lhs = targetMatrix(X_cal,:);
sample = X_cal;