function PRESS = press_check(X,y)
% Function returns PRESS for model
% Reference: http://reliawiki.org/index.php/Multiple_Linear_Regression_Analysis
% These equations are not yet directly used by the code, and they are not
% really even functions, but scripts that hold different statistical
% analysis equations.

%INPUTS:
%  X = Matrix of predictor variables. Each row is different observation, each column is different predictor
%  y  =  Target variable for regression.  Each row is observation

%OUTPUTS:
% PRESS = PRESS residuals


%% Multiple Linear Regression
% Instead of calculating the coefficient matrix using pinv or backslash,
% the coeffciients will be solved for the normal way, that is:
% b = (X' X)^-1 X y
% NOTE: IF THIS METHOD THROWS UP WARNINGS ABOUT ILL-CONDITIONED, DO NOT
% SUPPRESS THEM, IT IS IMPORTANT INFORMATION TO HAVE.
invXtX = pinv(X'*X);

%%


% The approximation can then be made with
% y_hat = X b = X (X' X)^-1 X' y = H y,
% where H is called the "hat" matrix, H = X(X'X)^-1X'
H = X*invXtX*X';
y_hat = H*y;

% The error e = y - y_hat
e = y - y_hat;

%% Measures of Model Adequacy
% Including R^2, adjusted R^2, PRESS, and PRESS R^2

% R^2 is calculated as the ratio of SSR to SST, or 1 - (SSE/SST), in either
% case, we'll need SST, which is just SSE + SSR;
% SST = SSE + SSR;


% The PRESS residuals can be obtained using the diagonal elements of the
% hat matrix.
e_p = e./(1-diag(H));
PRESS = sum(e_p.^2);

% PRESS R_sq
% R_sq_p = 1 - (PRESS/SST);


