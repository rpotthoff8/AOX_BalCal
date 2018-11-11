%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
%
function y=balCal_meritFunction2(w,residvec,etak)

   y=dot(residvec,residvec);

   u=exp(etak*log(abs(w)));
   b=dot(u,residvec);
   z=dot(u,u);
   
   p = residvec - (b/z)*u;

   y = dot(p,p);
   
end   

%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
%