function V=vif_dl2(X)
% V=vif(X)
%
% Variance inflaction factor in Regression Analysis
% the variance inflation factor (VIF) quantifies the severity of 
% multicollinearity in an ordinary least squares regression analysis. 
%It provides an index that measures how much the variance of an estimated 
% regression coefficient is increased because of collinearity.
%
% INPUT:
% 
% X is the matrix n onservation x p variables  
%
% OUTPUT:
%
% V is a column vector of vif indices

V=diag(pinv(corr(X)));
end
