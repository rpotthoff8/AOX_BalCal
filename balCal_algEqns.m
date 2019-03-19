% Copyright ©2017 Andrew Meade, Javier Villarreal.  All Rights Reserved.
function in_comb = balCal_algEqns(model_FLAG,in,series)

% Detect the size of the input
n = size(in,1); %number of data points
d = size(in,2); %data dimensionality.
if model_FLAG == 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LINEAR MODEL
    in_comb = in;%                                                        3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    in_comb = interceptTerms(in_comb,series);
    return 
end

in_sq = in.^2;
j = 1;
ini_inj = zeros(n,(d^2-d)/2);
for k = 1:d-1
    for m = k+1:d
        ini_inj(:,j) = in(:,k).*in(:,m);
        j = j+1;
    end
end
if model_FLAG == 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TRUNCATED MODEL
    in_comb = [in,      ...                                               2
               in_sq,   ...
               ini_inj];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    in_comb = interceptTerms(in_comb,series);
    return 
end

abs_in = abs(in);
in_absin = in.*abs_in;
abs_iniinj = abs(ini_inj);
j = 1;
ini_absinj = zeros(n,(d^2-d)/2);
absini_inj = zeros(n,(d^2-d)/2);
for k = 1:d-1
    for m = k+1:d
        ini_absinj(:,j) = in(:,k).*abs_in(:,m);
        absini_inj(:,j) = abs_in(:,k).*in(:,m);
        j = j+1;
    end
end
in_cu = in.^3;
abs_incu = abs(in_cu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FULL MODEL
in_comb = [in,        ...                                                 1
           abs_in,    ...
           in_sq,     ...
           in_absin,  ...
           ini_inj,   ...
           abs_iniinj,...
           ini_absinj,...
           absini_inj,...
           in_cu,     ...
           abs_incu];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_comb = interceptTerms(in_comb,series);
end

function in_comb = interceptTerms(in_comb,series)
n = size(in_comb,1);
[~,s_1st,s_id] = unique(series);
nseries = length(s_1st);
ints = zeros(n,nseries);
ids = sub2ind(size(ints),[1:n]',s_id);
ints(ids) = 1;
in_comb = [in_comb, ints];
end