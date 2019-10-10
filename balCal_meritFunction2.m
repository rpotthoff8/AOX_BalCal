%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
%
function y=balCal_meritFunction2(eps,centerI,dainputscalib,h,dimFlag,residvec)
s=dimFlag;

centerI_sub=sub2ind(size(dainputscalib),centerI,1:dimFlag);

dist=dainputscalib(centerI_sub)-dainputscalib;
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