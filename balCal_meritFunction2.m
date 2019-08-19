%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
%
function y=balCal_meritFunction2(xc_w,residvec,dainputs)
    w=xc_w(numel(xc_w));
    xc=xc_w(1:numel(xc_w)-1);
    
    adiffer=xc-dainputs;
    %     adiffervalid = dainputscalib(centerIndexHist(u,s),:)-dainputs;
    eta = dot(adiffer,adiffer,2);
    
%    y=dot(residvec,residvec);

   u=exp(eta*w);
   b=dot(u,residvec);
   z=dot(u,u);
   
   p = residvec - (b/z)*u;

   y = dot(p,p);
   
end   

%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
%