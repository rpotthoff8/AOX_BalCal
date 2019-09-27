%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
%
function y=balCal_meritFunction2(eps,centerI,R_square_whole,h,dimFlag,residvec)
s=dimFlag;
R_square=R_square_whole(:,centerI);
%    y=dot(residvec,residvec);

phi=((eps^s)/(sqrt(pi^s)))*exp(-((eps^2)*(R_square))/h^2); %From 'Iterated Approximate Moving Least Squares Approximation', Fasshauer and Zhang, Equation 22

b=dot(phi,residvec);
z=dot(phi,phi);

p = residvec - (b/z)*phi;

y = dot(p,p);

end

%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
%