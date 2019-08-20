%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
%
function resSquare2=balCal_meritFunction2(xc_w,residvec,dainputs,aprxINminGZ2,targetMatrix0,series0)
    w=xc_w(numel(xc_w));
    xci=round(xc_w(1:numel(xc_w)-1));
    
    adiffer=dainputs(xci,:)-dainputs;
    %     adiffervalid = dainputscalib(centerIndexHist(u,s),:)-dainputs;
    eta = dot(adiffer,adiffer,2);
    
%    y=dot(residvec,residvec);

   rbfINminGZ=exp(eta*w);
   c = (rbfINminGZ\residvec);
   
    rbfc_INminGZ = c*rbfINminGZ;
   
    aprxINminGZ2 = aprxINminGZ2+rbfc_INminGZ;
    taresAllPointsGRBF = meantare(series0,aprxINminGZ2-targetMatrix0);
    targetRes2 = targetMatrix0-aprxINminGZ2+taresAllPointsGRBF;      %0=b-Ax
    newRes2 = targetRes2'*targetRes2;
    resSquare2 = diag(newRes2);
%    b=dot(u,residvec);
%    z=dot(u,u);
%    
%    p = residvec - (b/z)*u;
% 
%    y = dot(p,p);
   
end   
