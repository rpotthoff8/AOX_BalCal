
%
% Copyright ©2017 Andrew Meade.  All Rights Reserved.
%
function [comINvec,comLZvec,comGZvec]=balCal_algEquations3(themodel_FLAG,numeqns, ndim,numdatapts,numseries,thelasttare,inputmatrix,lzmatrix,gzmatrix)
%
   y=0;

%'made it'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if themodel_FLAG == 1;
%start of full set

%'made it this far'


%make the 160 combinations

looper = 1;
for m=1:ndim
comINvec(looper,:) = inputmatrix(:,m);
comLZvec(looper,:) = lzmatrix(:,m);
comGZvec(looper,:) = gzmatrix(:,m);
looper = looper+1;
end


for m=1:ndim
comINvec(looper,:) = abs(inputmatrix(:,m));
comLZvec(looper,:) = abs(lzmatrix(:,m));
comGZvec(looper,:) = abs(gzmatrix(:,m));
looper = looper+1;
end


for m=1:ndim
comINvec(looper,:) = inputmatrix(:,m).^2;
comLZvec(looper,:) = lzmatrix(:,m).^2;
comGZvec(looper,:) = gzmatrix(:,m).^2;
looper = looper+1;
end


for m=1:ndim
comINvec(looper,:) = inputmatrix(:,m).*abs(inputmatrix(:,m));
comLZvec(looper,:) = lzmatrix(:,m).*abs(lzmatrix(:,m));
comGZvec(looper,:) = gzmatrix(:,m).*abs(gzmatrix(:,m));
looper = looper+1;
end


for k=1:ndim-1
  for m=k+1:ndim
comINvec(looper,:) = inputmatrix(:,k).*inputmatrix(:,m);
comLZvec(looper,:) = lzmatrix(:,k).*lzmatrix(:,m);
comGZvec(looper,:) = gzmatrix(:,k).*gzmatrix(:,m);
looper = looper+1;
  end
end


for k=1:ndim-1
  for m=k+1:ndim
comINvec(looper,:) = abs((inputmatrix(:,k)).*(inputmatrix(:,m)));
comLZvec(looper,:) = abs((lzmatrix(:,k)).*(lzmatrix(:,m)));
comGZvec(looper,:) = abs((gzmatrix(:,k)).*(gzmatrix(:,m)));
looper = looper+1;
  end
end


for k=1:ndim-1
  for m=k+1:ndim
comINvec(looper,:) = inputmatrix(:,k).*abs(inputmatrix(:,m));
comLZvec(looper,:) = lzmatrix(:,k).*abs(lzmatrix(:,m));
comGZvec(looper,:) = gzmatrix(:,k).*abs(gzmatrix(:,m));
looper = looper+1;
  end
end


for k=1:ndim-1
  for m=k+1:ndim
comINvec(looper,:) = abs(inputmatrix(:,k)).*(inputmatrix(:,m));
comLZvec(looper,:) = abs(lzmatrix(:,k)).*(lzmatrix(:,m));
comGZvec(looper,:) = abs(gzmatrix(:,k)).*(gzmatrix(:,m));
looper = looper+1;
  end
end


  for m=1:ndim
comINvec(looper,:) = inputmatrix(:,m).^3;
comLZvec(looper,:) = lzmatrix(:,m).^3;
comGZvec(looper,:) = gzmatrix(:,m).^3;
looper = looper+1;
  end


  for m=1:ndim
comINvec(looper,:) = abs(inputmatrix(:,m).^3);
comLZvec(looper,:) = abs(lzmatrix(:,m).^3);
comGZvec(looper,:) = abs(gzmatrix(:,m).^3);
looper = looper+1;
  end

numeqns = looper-1; 


%end of full term set
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if themodel_FLAG == 2;

%start of truncated set


looper = 1;
for m=1:ndim
comINvec(looper,:) = inputmatrix(:,m);
comLZvec(looper,:) = lzmatrix(:,m);
comGZvec(looper,:) = gzmatrix(:,m);
looper = looper+1;
end


for m=1:ndim
comINvec(looper,:) = inputmatrix(:,m).^2;
comLZvec(looper,:) = lzmatrix(:,m).^2;
comGZvec(looper,:) = gzmatrix(:,m).^2;
looper = looper+1;
end


for k=1:ndim-1
  for m=k+1:ndim
comINvec(looper,:) = inputmatrix(:,k).*inputmatrix(:,m);
comLZvec(looper,:) = lzmatrix(:,k).*lzmatrix(:,m);
comGZvec(looper,:) = gzmatrix(:,k).*gzmatrix(:,m);
looper = looper+1;
  end
end


numeqns = looper-1; 

%end of truncated term set
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if themodel_FLAG == 3;

%start of linear set


looper = 1;
for m=1:ndim
comINvec(looper,:) = inputmatrix(:,m);
comLZvec(looper,:) = lzmatrix(:,m);
comGZvec(looper,:) = gzmatrix(:,m);
looper = looper+1;
end

numeqns = looper-1; 


%end of linear term set
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if thelasttare > 0

%% Effectively removes the original tare values so we can calculate the averages
for i=1:thelasttare
comINvec(numeqns+i,:) = 0;
comLZvec(numeqns+i,:) = 0;
comGZvec(numeqns+i,:) = 0;
end

%%

for loopk=1:numdatapts
           comLZvec(numeqns+1,loopk) = 1.0;
end

%%
%%

    end


% end of the subroutine
end   

%
% Copyright ©2017 Andrew Meade.  All Rights Reserved.
%


