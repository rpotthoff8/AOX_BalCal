function customMatrix=customMatrix_builder(voltdimFlag,termInclude,loaddimFlag, interceptFlag)
%Function creates a custom matrix for which algebraic terms should be
%included in the model.  This matrix of 1s and 0s denotes the terms
%selected.  The possible algebraic terms are listed in order below.
%Function works for any data dimension (any # of channels)

%INPUTS:
%  voltdimFlag = Number of voltage data channels
%  termInclude = Boolean 1x12 vector for which algebraic term combinations should be included in model
%  loaddimFlag = Number of load data channels
%  interceptFlag = Flag for if global intercept should be included in model 

%OUTPUTS:
%  customMatrix = Custom Equation Matrix for calculations

%First determine the index numbers for each term in the full equation set
termIndex=cell(12,1);
%Terms are listed in following order:
%  F, |F|, F*F, F*|F|, F*G, |F*G|, F*|G|, |F|*G, F*F*F, |F*F*F|, F*G*G, F*G*H
termIndex{1}=(1:voltdimFlag)'; % F
termIndex{2}=(1:voltdimFlag)'+max(termIndex{1}); % |F|
termIndex{3}=(1:voltdimFlag)'+max(termIndex{2}); % F*F
termIndex{4}=(1:voltdimFlag)'+max(termIndex{3}); % F*|F|
termIndex{5}=(1:(voltdimFlag^2-voltdimFlag)/2)'+max(termIndex{4}); % F*G
termIndex{6}=(1:(voltdimFlag^2-voltdimFlag)/2)'+max(termIndex{5}); % |F*G|
termIndex{7}=(1:(voltdimFlag^2-voltdimFlag)/2)'+max(termIndex{6}); % F*|G|
termIndex{8}=(1:(voltdimFlag^2-voltdimFlag)/2)'+max(termIndex{7}); % |F|*G
termIndex{9}=(1:voltdimFlag)'+max(termIndex{8}); % F*F*F
termIndex{10}=(1:voltdimFlag)'+max(termIndex{9}); % |F*F*F|
termIndex{11}=(1:(factorial(voltdimFlag)/factorial(voltdimFlag-2)))'+max(termIndex{10}); % F*G*G
termIndex{12}=(1:(factorial(voltdimFlag)/(factorial(3)*factorial(voltdimFlag-3))))'+max(termIndex{11}); % F*G*H

%Build a custom equation matrix based on the terms selected and the data
%dimension
nterms = 2*voltdimFlag*(voltdimFlag+2)+factorial(voltdimFlag)/factorial(voltdimFlag-2)+factorial(voltdimFlag)/(factorial(3)*factorial(voltdimFlag-3));
customMatrix=zeros(nterms,loaddimFlag);
customMatrix(cell2mat(termIndex(boolean(termInclude))),:)=1;

%Add intercept flag to customMatrix
if interceptFlag==1
    intercept=ones(1,loaddimFlag);
else
    intercept=zeros(1,loaddimFlag);
end
customMatrix=[intercept; customMatrix];
end