% Statistical function collection.
% Reference: http://reliawiki.org/index.php/Multiple_Linear_Regression_Analysis
% These equations are not yet directly used by the code, and they are not
% really even functions, but scripts that hold different statistical
% analysis equations.

% A breakpoint can be included in AOX_BalCal.m at approximately line 198,
% (right before xcalib_k is calculated). The scripts in this file will work
% once the workspace has been populated to that point.
% (For notation convenience, we'll store the loop iteration variable in
% iter.}
iter = k; clear k;

% Before beginning, the predictors and the responses will be re-named with
% new variables for clarity.
X = comIN_k;
y = targetMatrix(:,iter);

%% Multiple Linear Regression
% Instead of calculating the coefficient matrix using pinv or backslash,
% the coeffciients will be solved for the normal way, that is:
% b = (X' X)^-1 X y
% NOTE: IF THIS METHOD THROWS UP WARNINGS ABOUT ILL-CONDITIONED, DO NOT
% SUPPRESS THEM, IT IS IMPORTANT INFORMATION TO HAVE.
beta = inv(X'*X)*X'*y;

% The approximation can then be made with
% y_hat = X b = X (X' X)^-1 X' y = H y,
% where H is called the "hat" matrix, H = X(X'X)^-1X'
H = X*inv(X'*X)*X';
y_hat = H*y;

% The error e = y - y_hat
e = y - y_hat;

%% Covariance matrix and standard error
% The covariance amtrix is calculated as C = sigma_hat^2 (X' X)^-1
% sigma_hat is an estimate of the variance (std dev), based on the MSE
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
MSE = SSE/dof_e;
sigma_hat_sq = MSE;

C = sigma_hat_sq*inv(X'*X);

% The standard error of each variable is then the square root of it's
% variance (the diagonal terms in the covariance matrix C.
se = sqrt(diag(C));

%% Hypothesis testing (t-testing and F-testing)

% The F statistic is defined as F_0 = MSR/MSE, where MSR is the mean square
% of the regression.
%First we'll need the mean value of y.
y_bar = sum(y)/n;  %I'm aware mean() exists, writing equations out for clarity)
SSR = sum((y_hat-y_bar).^2);

% NOTE: for MSR, the dof = k (not including intercept term)
dof_r = k-1;
MSR = SSR/dof_r;

%F-statistic, F = MSR/MSE
F = MSR/MSE;
F_cr = finv(0.95,dof_r,dof_e);
% F0 is the F-statistic of the Null Hypothesis, meaning that y is not a
% function of X at all. If F>F0, then we can reject the null hypothesis,
% and the data suggests that y is a function of X. in the function finv,
% the first input 0.05, is to calculate those tolerance with 95%
% confidence.
if F > F_cr
    disp('We can reject the null hypothesis');
else
    disp('We cannot reject the null hypothesis');
end

%% T test for individual coefficients b_i
% The t statistic is the ration of a coefficient to its standard error
% The input to tinv is 0.975 is because it is a two-sided test. so .025 is
% subtracted from either end, meaning subtracting 0.05 total from both ends
% -> 95% confidence.
T = beta./se;
T_cr = tinv(0.975,dof_e);
sig = ~((T<T_cr) & (T>-T_cr));
% If sig == 1, then that coefficient rejects the null hypothesis that b =
% 0, so the data suggests that the terms matters to this regression.
% NOTE: This is a test of statistical significant of a term IN A
% MODEL THAT INCLUDES ALL THE TERMS. In other words, the t-test might have
% different results depending on which all terms were included in the
% regression. Obviously multicollinearity will greatly affect how
% statistically significant or insignificant a term appears to be, given
% differen data sets.