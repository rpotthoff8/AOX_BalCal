function [VIF_met, VIF_max, sig_all, P_max ,search_metric,VIF,sig,p_T]=test_combo(X,y,pct, VIFthresh, nseries, FLAGS)
% Statistical function collection.
% Reference: http://reliawiki.org/index.php/Multiple_Linear_Regression_Analysis
% These equations are not yet directly used by the code, and they are not
% really even functions, but scripts that hold different statistical
% analysis equations.

%INPUTS:
%  X = Matrix of predictor variables. Each row is different observation, each column is different predictor
%  y  =  Target variable for regression.  Each row is observation
%  pct  =  Percent confidence level for ANOVA calculations
%  VIFthresh = Threshold for max allowed VIF
%  nseries  =  Number of series in dataset
%  search_metric_flag = Flag for search metric to return: 1=PRESS (minimize),
%       2=Square Root of Residual Mean Square (minimize), 3=F-Value (maximize)
%  VIF_stop_flag = Flag for if calculations should be stopped if VIF is too high


%OUTPUTS:
% VIF_met = Boolean: if all terms have VIF <= VIFthresh
% VIF_max = Value of maximum VIF
% sig_all = Boolean: if all terms (not including tare intercepts) are significant
% P_max = Value of max coefficient P value
% search_metric = Return value for search metric

search_metric_flag=FLAGS.search_metric;
VIF_stop_flag=FLAGS.VIF_stop;

%% Multicollinearity detection / VIF calcualtion
if FLAGS.glob_intercept==0 %If no global intercept term
VIF=vif_dl(X);
else
    VIF=vif_dl(X(:,2:end)); %Calculate VIF for all terms except global intercept
    VIF=[1,VIF]; %Set VIF for intercept=1
end

VIF_max=max(VIF,[],'all');
if VIF_max<=VIFthresh
    VIF_met=true;
else
    VIF_met=false;
    if VIF_stop_flag==1 %If haulting calculations for VIF> VIFthresh
        sig_all = false;
        P_max = 0;
        search_metric= Inf;
        sig=0;
        p_T=0;
        return
    end
end



%% Multiple Linear Regression
% Instead of calculating the coefficient matrix using pinv or backslash,
% the coeffciients will be solved for the normal way, that is:
% b = (X' X)^-1 X y
% NOTE: IF THIS METHOD THROWS UP WARNINGS ABOUT ILL-CONDITIONED, DO NOT
% SUPPRESS THEM, IT IS IMPORTANT INFORMATION TO HAVE.
invXtX = pinv(X'*X);
beta = invXtX*X'*y;

%% TEMPORARY
% The current data is extremely ill-conditioned, so warnings are semi-suppressed
% for now if the data is bad, since it can be very repetitive.
% This will be removed later, since those
% warnings can be useful.
% warnMsg = lastwarn;
% if ~isempty(warnMsg)
%     warning('off','MATLAB:nearlySingularMatrix');
%     lastwarn('');
% end
%%

%% Multicollinearity detection / VIF calcualtion
% Variance inflation factor is calculated first so that when we're
% iterating through test terms in the model, we can eliminate collinear
% data right away, before wasting computer time on unnecessary statistical
% calculations.
%VIF calculations performed if test_FLAG=0: Do not perform during iterated
%recommended equation solving for time saving
% if noVIF_FLAG==0
%     VIF = vif(X);


%
%     if any(VIF>=10)
%         warning('VIF calculation indicates strong multicollinearity. Analysis of Variance results cannot be trusted.')
%     elseif any(VIF>=4)
%         warning('VIF calculation indicates some multicollinearity. Analysis of Variance results may be inaccurate.')
%     end
% else
%     VIF={'VIF NOT CALCULATED'};
%     VIF=-1;
% end
% if test_FLAG == 1 && max(VIF) >= 4
%     ANOVA.test = -1;
%     return
% end

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

if search_metric_flag==2
    search_metric=sqrt(MSE);
end

C = sigma_hat_sq*invXtX;

% The standard error of each variable is then the square root of it's
% variance (the diagonal terms in the covariance matrix C.
se = sqrt(diag(C));

%% Hypothesis testing (t-testing and F-testing)
if search_metric_flag==3
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
    search_metric=-F;
end
% p_F = 1 - fcdf(F,dof_r,dof_e);
% F_cr = finv(pct/100,dof_r,dof_e);
% % F0 is the F-statistic of the Null Hypothesis, meaning that y is not a
% % function of X at all. If F>F0, then we can reject the null hypothesis,
% % and the data suggests that y is a function of X. in the function finv,
% % the first input 0.05, is to calculate those tolerance with 95%
% % confidence.
% if F > F_cr
%     %disp('We can reject the null hypothesis');
% else
%     disp('We cannot reject the null hypothesis');
% end

%% T test for individual coefficients b_i
% The t statistic is the ration of a coefficient to its standard error
% The input to tinv is 0.975 is because it is a two-sided test. so .025 is
% subtracted from either end, meaning subtracting 0.05 total from both ends
% -> 95% confidence.
T = beta./se;
p_T = 2*(1 - tcdf(abs(T),dof_e));
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

%% Confidence Intervals

% The confidence intervals for the coefficients are calculated as
% beta_CI = T_cr.*se;

% confidence interval for the predicted values, must be calculated one data
% point at a time
% for j = 1:n
%     y_hat_CI(j,1) = T_cr*sqrt(sigma_hat_sq*X(j,:)*invXtX*X(j,:)');
% end

%% Prediction Intervals
% Prediction intervals are different from confidence intervals. Confidence
% intervals are about how well we believe we regressed the given
% (calibration) data, prediction intervals would be giving a confidence
% interval for a new prediction, given new data.
%
% for j = 1:n
%     y_hat_PI(j,1) = T_cr*sqrt(sigma_hat_sq*(1+(X(j,:)*invXtX*X(j,:)')));
% end

%% Measures of Model Adequacy
% Including R^2, adjusted R^2, PRESS, and PRESS R^2

% R^2 is calculated as the ratio of SSR to SST, or 1 - (SSE/SST), in either
% case, we'll need SST, which is just SSE + SSR;
% SST = SSE + SSR;

% This bit has been commented out because it's no longer necessary.
% % to verify, we'll also calcualte SST the long way.
% SST_v = sum((y - y_bar).^2);
% % In the data set I used, I got a difference of O(e-5), so good enough.

% R_sq = 1 - SSE/SST;

% Since R_sq increases as more terms are added to the model, adjusted R_sq
% takes degrees of freedom into account. = 1 - MSE/MST
% dof_t = n-1;
% MST = SST/dof_t;
% R_sq_adj = 1 - MSE/MST;

% The PRESS residuals can be obtained using the diagonal elements of the
% hat matrix.
if search_metric_flag==1
    e_p = e./(1-diag(H));
    PRESS = sum(e_p.^2);
    search_metric=PRESS;
end

% % PRESS R_sq
% R_sq_p = 1 - (PRESS/SST);

%% OUTPUTS


P_max=max(p_T);
if FLAGS.tare_intercept==1
    sig_all=all(sig(1:size(sig,1)-nseries)); %Check if all terms are significant, excluding tare intercepts
else
    sig_all=all(sig); %Check if all terms are significant
end



end