function ANOVA = anova(X,y)

fprintf('\nCalculating The Statistics ...\n')
% Statistical function collection.
% Reference: http://reliawiki.org/index.php/Multiple_Linear_Regression_Analysis
% These equations are not yet directly used by the code, and they are not
% really even functions, but scripts that hold different statistical
% analysis equations.

% % A breakpoint can be included in calc_xcalib.m at approximately line 34,
% % (right before xcalib_k is calculated). The scripts in this file will work
% % once the workspace has been populated to that point.
% % (For notation convenience, we'll store the loop iteration variable in
% % iter.}
% iter = k; clear k;
% 
% % Before beginning, the predictors and the responses will be re-named with
% % new variables for clarity.
% X = comIN_k;
% y = targetMatrix(:,iter);

%% Multiple Linear Regression
% Instead of calculating the coefficient matrix using pinv or backslash,
% the coeffciients will be solved for the normal way, that is:
% b = (X' X)^-1 X y
% NOTE: IF THIS METHOD THROWS UP WARNINGS ABOUT ILL-CONDITIONED, DO NOT
% SUPPRESS THEM, IT IS IMPORTANT INFORMATION TO HAVE.
beta = inv(X'*X)*X'*y;

%% TEMPORARY
% The current data is extremely ill-conditioned, so warnings are semi-suppressed
% for now if the data is bad, since it can be very repetitive.
% This will be removed later, since those
% warnings can be useful.
warnMsg = lastwarn;
if ~isempty(warnMsg)
    warning('off','MATLAB:nearlySingularMatrix');
    lastwarn('');
end
%%

% The approximation can then be made with
% y_hat = X b = X (X' X)^-1 X' y = H y,
% where H is called the "hat" matrix, H = X(X'X)^-1X'
H = X*inv(X'*X)*X';
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
p_F = 1 - fcdf(F,dof_r,dof_e);
F_cr = finv(0.95,dof_r,dof_e);
% F0 is the F-statistic of the Null Hypothesis, meaning that y is not a
% function of X at all. If F>F0, then we can reject the null hypothesis,
% and the data suggests that y is a function of X. in the function finv,
% the first input 0.05, is to calculate those tolerance with 95%
% confidence.
if F > F_cr
    disp('We can reject the null hypothesis: The residual is independent from the voltage');
else
    disp('We cannot reject the null hypothesis: The residual is related to the voltage');
end

%% T test for individual coefficients b_i
% The t statistic is the ration of a coefficient to its standard error
% The input to tinv is 0.975 is because it is a two-sided test. so .025 is
% subtracted from either end, meaning subtracting 0.05 total from both ends
% -> 95% confidence.
T = beta./se;
p_T = 2*(1 - tcdf(abs(T),dof_e));
T_cr = tinv(0.975,dof_e);
sig = ~((T<T_cr) & (T>-T_cr));
% If sig == 1, then that coefficient rejects the null hypothesis that b =
% 0, so the data suggests that the terms matters to this regression.
% NOTE: This is a test of statistical significance of a term IN A
% MODEL THAT INCLUDES ALL THE TERMS. In other words, the t-test might have
% different results depending on which all terms were included in the
% regression. Obviously multicollinearity will greatly affect how
% statistically significant or insignificant a term appears to be, given
% different data sets.

%% Confidence Intervals

% The confidence intervals for the coefficients are calculated as
beta_CI = T_cr.*se;

% confidence interval for the predicted values, must be calculated one data
% point at a time
for j = 1:n
    y_hat_CI(j,1) = T_cr*sqrt(sigma_hat_sq*X(j,:)*inv(X'*X)*X(j,:)');
end

%% Prediction Intervals
% Prediction intervals are different from confidence intervals. Confidence
% intervals are about how well we believe we regressed the given
% (calibration) data, prediction intervals would be giving a confidence
% interval for a new prediction, given new data.

for j = 1:n
    y_hat_PI(j,1) = T_cr*sqrt(sigma_hat_sq*(1+(X(j,:)*inv(X'*X)*X(j,:)')));
end

%% Measures of Model Adequacy
% Including R^2, adjusted R^2, PRESS, and PRESS R^2

% R^2 is calculated as the ratio of SSR to SST, or 1 - (SSE/SST), in either
% case, we'll need SST, which is just SSE + SSR;
SST = SSE + SSR;

% This bit has been commented out because it's no longer necessary.
% % to verify, we'll also calcualte SST the long way.
% SST_v = sum((y - y_bar).^2);
% % In the data set I used, I got a difference of O(e-5), so good enough.

R_sq = 1 - SSE/SST;

% Since R_sq increases as more terms are added to the model, adjusted R_sq
% takes degrees of freedom into account. = 1 - MSE/MST
dof_t = n-1;
MST = SST/dof_t;
R_sq_adj = 1 - MSE/MST;

% The PRESS residuals can be obtained using the diagonal elements of the
% hat matrix.
e_p = e./(1-diag(H));
PRESS = sum(e_p.^2);

% PRESS R_sq
R_sq_p = 1 - (PRESS/SST);

%% Multicollinearity detection / VIF calcualtion

VIF = vif(X);

%% Save all the relevant variables to a structure to pass out of the function

ANOVA.beta     = beta;     % Coefficients
ANOVA.e        = e;        % Residuals

ANOVA.PRESS = PRESS;       % PRESS Statistic
ANOVA.F = F;               % F-value of regression
ANOVA.p_F = p_F;           % P-value of regression

ANOVA.R_sq = R_sq;         % R-square
ANOVA.R_sq_adj = R_sq_adj; % Adjusted R-Square
ANOVA.R_sq_p = R_sq_p;     % PRESS R-square

ANOVA.beta_CI  = beta_CI;  % Coefficient Confidence Intervals 
ANOVA.T = T;               % T-statistic of coefficients
ANOVA.p_T = p_T;           % P-value of coefficients
ANOVA.VIF = VIF;           % Variance Inflation Factors
ANOVA.sig = sig;            %If term is significant

ANOVA.y_hat_PI=y_hat_PI; %Prediction interval for new datapoints
%% EXPERIMENTAL
% Saving variables to calculate prediction intervals live in approximation
ANOVA.PI.T_cr = T_cr;
ANOVA.PI.sigma_hat_sq = sigma_hat_sq;
ANOVA.PI.X = X;
ANOVA.PI.calc_pi = "T_cr*sqrt(sigma_hat_sq*(1+(x*inv(X'*X)*x')))";

%% Coded values for polynomial regressions
% Values of the variables are coded by centering or expressing the levels
% of the variable as deviations from the mean value of the variable and
% then scaling or dividing the deviations obtained by half the range of the
% variable. The reason for using coded predictor variables is that many
% times x and x^2 are highly correlated and, if uncoded values are used,
% there may be computational difficulties while calculating the (X'X)^-1
% matrix to obtain the estimates, beta, of the regression coefficients
% using the equation for the F distribution.

function VIF = vif(X)

[~,k] = size(X);
VIF = zeros(k,1);
for j = 1:k
    y = X(:,j);
    X_j = X; X_j(:,j) = []; %x_vif = [ones(n_data,1),x_vif];
    H = X_j*inv(X_j'*X_j)*X_j';
    
    y_bar = mean(y);
    y_hat = H*y;
    
    SSE = sum((y-y_hat).^2);
    SST = sum((y-y_bar).^2);
    
    R_sq = 1 - SSE/SST;
    VIF(j,1) = 1/(1-R_sq);
end

%% Turn warnings back on
warning('on','all');