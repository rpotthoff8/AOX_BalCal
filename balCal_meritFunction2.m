%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
%
function y=balCal_meritFunction2(eps,centerLoc,dainputscalib,h,dimFlag,residvec)
s=dimFlag;

dist=centerLoc'-dainputscalib;
R_square=sum(dist.^2,2); %Eqn 17 from Javier's notes: squared distance between each point
    
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