% Copyright �2017 Andrew Meade, Javier Villarreal.  All Rights Reserved.
function in_comb = balCal_algEqns(model_FLAG,in,series,intercept_FLAG,normFLAG)
%validation_FLAG==1 if producing in_comb for the validation section.  The
%only difference is the series intercept terms are not included for
%validation

if nargin <5
    normFLAG = 0;
end

if normFLAG == 1
    range = max(in) - min(in);
    shift = min(in) + range/2;
    in_n = (in - shift)./(range/2);
else
    in_n = in;
end

% Detect the size of the input
n = size(in_n,1); %number of data points
d = size(in_n,2); %data dimensionality.
if model_FLAG == 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LINEAR MODEL
    in_comb = in_n;%                                                        3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    in_comb = interceptTerms(in_comb,series,intercept_FLAG);
    return
end

in_sq = in_n.^2;
j = 1;
ini_inj = zeros(n,(d^2-d)/2);
for k = 1:d-1
    for m = k+1:d
        ini_inj(:,j) = in_n(:,k).*in_n(:,m);
        j = j+1;
    end
end
if model_FLAG == 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TRUNCATED MODEL
    in_comb = [in_n,      ...                                               2
               in_sq,   ...
               ini_inj];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    in_comb = interceptTerms(in_comb,series,intercept_FLAG);
    return
end


abs_in = abs(in);
if normFLAG == 1
    range = max(abs_in) - min(abs_in);
    shift = min(abs_in) + range/2;
    abs_in = (abs_in - shift)./(range/2);
end

in_absin = in_n.*abs_in;
j = 1;
abs_iniinj = zeros(n,(d^2-d)/2);
ini_absinj = zeros(n,(d^2-d)/2);
absini_inj = zeros(n,(d^2-d)/2);
for k = 1:d-1
    for m = k+1:d
        abs_iniinj(:,j) = abs_in(:,k).*abs_in(:,m);
        ini_absinj(:,j) = in_n(:,k).*abs_in(:,m);
        absini_inj(:,j) = abs_in(:,k).*in_n(:,m);
        j = j+1;
    end
end
in_cu = in_n.^3;
abs_incu = abs_in.^3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FULL MODEL
in_comb = [in_n,        ...                                                 1
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
in_comb = interceptTerms(in_comb,series,intercept_FLAG);
end

function in_comb = interceptTerms(in_comb,series,intercept_FLAG)
n = size(in_comb,1);
[~,s_1st,s_id] = unique(series);
nseries = length(s_1st);
ints = zeros(n,nseries);
ids = sub2ind(size(ints),[1:n]',s_id);
ints(ids) = 1;

if intercept_FLAG==1 %If preference for series intercepts is on, add 1's for intercept terms: (These are not included for validation)
    in_comb = [in_comb, ints];
end
end
