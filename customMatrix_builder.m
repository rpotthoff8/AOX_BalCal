function customMatrix=customMatrix_builder(dimFlag,termInclude)
%Function creates a custom matrix for which algebraic terms should be
%included in the model.  This matrix of 1s and 0s denotes the terms
%selected.  The possible algebraic terms are listed in order below.
%Function works for any data dimension (any # of channels)

%INPUTS:
%  dimFlag = Number of data channels
%  termInclude = Boolean 1x10 vector for which algebraic term combinations should be included in model

%OUTPUTS:
%  customMatrix = Custom Equation Matrix for calculations

%First determine the index numbers for each term in the full equation set
termIndex=cell(10,1);
%Terms are listed in following order:
% F, |F|, F*F, F*|F|, F*G, |F*G|, F*|G|, |F|*G, F*F*F, |F*F*F|
termIndex{1}=(1:dimFlag)'; % F
termIndex{2}=(1:dimFlag)'+max(termIndex{1}); % |F|
termIndex{3}=(1:dimFlag)'+max(termIndex{2}); % F*F
termIndex{4}=(1:dimFlag)'+max(termIndex{3}); % F*|F|
termIndex{5}=(1:(dimFlag^2-dimFlag)/2)'+max(termIndex{4}); % F*G
termIndex{6}=(1:(dimFlag^2-dimFlag)/2)'+max(termIndex{5}); % |F*G|
termIndex{7}=(1:(dimFlag^2-dimFlag)/2)'+max(termIndex{6}); % F*|G|
termIndex{8}=(1:(dimFlag^2-dimFlag)/2)'+max(termIndex{7}); % |F|*G
termIndex{9}=(1:dimFlag)'+max(termIndex{8}); % F*F*F
termIndex{10}=(1:dimFlag)'+max(termIndex{9}); % |F*F*F|

%Build a custom equation matrix based on the terms selected and the data
%dimension
customMatrix=zeros(2*dimFlag*(dimFlag+2),dimFlag);
customMatrix(cell2mat(termIndex(boolean(termInclude))),:)=1;

end