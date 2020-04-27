%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
function y=balCal_meritFunction2(eps,residvec,R_square,h,dimFlag)
%Function is used in optimizing properties for RBFs to be added. Object of
%optimization is to minimize y

%INPUTS:
%  eps = Epsilon value: determines RBF width
%  residvec  =  Residuals attempting to fit
%  R_square  =  Matrix of distance squared between center and datapoints
%  h  =  RBF constant based on point spacing
%  dimFlag  =  Number of channels for data

%OUTPUTS:
%  y = Merit Value for success in RBF fitting residuals.  Object of optimization is to minimize y

s=dimFlag;
%    y=dot(residvec,residvec);

phi=exp(-((eps^2)*(R_square))/h^2); %From 'Iterated Approximate Moving Least Squares Approximation', Fasshauer and Zhang, Equation 22
% phi=((eps^s)/(sqrt(pi^s)))*exp(-((eps^2)*(R_square))/h^2); %From 'Iterated Approximate Moving Least Squares Approximation', Fasshauer and Zhang, Equation 22
% phi=phi-mean(phi); %Phi bias is mean value for phi

% %%OLD
% b=dot(phi,residvec); %Projection of phi onto residual
% z=dot(phi,phi); %Square of the magnitude of phi vector
% 
% p = residvec - (b/z)*phi; %Measure for remaining residual after adding phi
% 
% y = dot(p,p); % Scalar value to be minimized

%%NEW: 11 Feb 20
comIN=[ones(length(phi) ,1) ,phi];

c=comIN\residvec; %Solve for coefficient to best fit RBF to residuals

resid_remain=residvec-comIN*c; %Remainder of residual after adding new basis
y=dot(resid_remain,resid_remain); %RMS of remaining residual
end

%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
%