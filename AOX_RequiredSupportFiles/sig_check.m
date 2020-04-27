function [sig,p_T]=sig_check(X,y,pct)
% Function returns vector of significance for each predictor term
% Reference: http://reliawiki.org/index.php/Multiple_Linear_Regression_Analysis

%INPUTS:
%  X = Matrix of predictor variables. Each row is different observation, each column is different predictor
%  y  =  Target variable for regression.  Each row is observation
%  pct  =  Percent confidence level for ANOVA calculations

%OUTPUTS:
% sig = vector for if terms are significant
% p_T = P-value of coefficients

%% Multiple Linear Regression
% Instead of calculating the coefficient matrix using pinv or backslash,
% the coeffciients will be solved for the normal way, that is:
% b = (X' X)^-1 X y
% NOTE: IF THIS METHOD THROWS UP WARNINGS ABOUT ILL-CONDITIONED, DO NOT
% SUPPRESS THEM, IT IS IMPORTANT INFORMATION TO HAVE.
invXtX = pinv(X'*X);
beta = invXtX*X'*y;

% The approximation can then be made with
% y_hat = X b = X (X' X)^-1 X' y = H y,
% where H is called the "hat" matrix, H = X(X'X)^-1X'
H = X*invXtX*X';
y_hat = H*y;

% The error e = y - y_hat
e = y - y_hat;

%% Covariance matrix and standard error
% The covariance matrix is calculated as C = sigma_hat^2 (X' X)^-1
% sigma_hat is an estimate of the variance (std dev), based on the MSE
%%
% (mean square error).
% MSE = SSE / dof(SSE)
% SSE is the sum square error, and dof() is the degrees of freedom of that
% error. The degrees of freedom will equal the number of points minus the
% number of coefficients that the regression is prediction.

%First, we'll calculate the SSE,
SSE = sum(e.^2);
% then, to calculate the degrees of freedom, we'll need the size
% information of the predictor matrix X
[n,k] = size(X);
dof_e = n-k;
if dof_e<0
    fprintf('WARNING: Negative Degrees of Freedom in ANOVA \n')
end
MSE = SSE/dof_e;
sigma_hat_sq = MSE;

C = sigma_hat_sq*invXtX;

% The standard error of each variable is then the square root of it's
% variance (the diagonal terms in the covariance matrix C.
se = sqrt(diag(C));

%% T test for individual coefficients b_i
% The t statistic is the ration of a coefficient to its standard error
% The input to tinv is 0.975 is because it is a two-sided test. so .025 is
% subtracted from either end, meaning subtracting 0.05 total from both ends
% -> 95% confidence.
T = beta./se;
p_T = 2*(1 - tcdf(abs(T),dof_e)); % P-value of coefficients
T_cr = tinv(1-((1-(pct/100))/2),dof_e);
sig = ~((T<T_cr) & (T>-T_cr));
% If sig == 1, then that coefficient rejects the null hypothesis that b =
% 0, so the data suggests that the terms matters to this regression.
% NOTE: This is a test of statistical significance of a term IN A
% MODEL THAT INCLUDES ALL THE TERMS. In other words, the t-test might have
% different results depending on which all terms were included in the
% regression. Obviously multicollinearity will greatly affect how
% statistically significant or insignificant a term appears to be, given
% different data sets.

end