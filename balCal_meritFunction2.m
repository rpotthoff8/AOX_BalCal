%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
%
function y=balCal_meritFunction2(eps,residvec,R_square,h,dimFlag)
s=dimFlag;
%    y=dot(residvec,residvec);

phi=((eps^s)/(sqrt(pi^s)))*exp(-((eps^2)*(R_square))/h^2); %From 'Iterated Approximate Moving Least Squares Approximation', Fasshauer and Zhang, Equation 22
phi_b=phi-mean(phi,1); %Phi with bias term subtracted
b=dot(phi_b,residvec);
z=dot(phi_b,phi_b);

p = residvec - (b/z)*phi_b;

y = dot(p,p);

end

%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
%