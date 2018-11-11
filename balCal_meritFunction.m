%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
%
function y=balCal_meritFunction(w,residvec,etak,etaLZ)
%function y=restest_RBF(w,A,B,aeq,beq,lb,ub,nonlncon,opt,RR0,residvec,etak)
   y=0;

   u=exp(etak*log(abs(w))) - exp(etaLZ*log(abs(w)));
   %u=exp((etak-etaLZ)*log(abs(w)));
   %u=exp(etak*log(abs(w)));
   b=dot(u,residvec);
   z=dot(u,u);
   
   y = residvec - (b/z)*u;
   %y = y'*y;
   %y = diag(y);
   y = dot(y,y);
   
end   

%
% Copyright ©2016 Ali Arya Mokhtarzadeh.  All Rights Reserved.
%