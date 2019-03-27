% Copyright ©2017 Andrew Meade.  All Rights Reserved.
%Function returns values needed for calculating partial derivatives for
%uncertainty propagation
function uncert_comIN=balCal_algEquations_partialdiff(themodel_FLAG, ndim, inputmatrix)

if themodel_FLAG == 1;
%start of full set
%make the 160 combinations
looper = 1;
for m=1:ndim
comINvec(looper,:) = inputmatrix(:,m);
for n=1:ndim %ADD
    if m==n
    uncert_comIN(looper,:,n)=ones(1,size(comINvec,2));
    else
    uncert_comIN(looper,:,n)=zeros(1,size(comINvec,2));    
    end
end %END
looper = looper+1;
end

for m=1:ndim
comINvec(looper,:) = abs(inputmatrix(:,m));
for n=1:ndim %ADD
    if m==n
    uncert_comIN(looper,:,n)=sign(inputmatrix(:,m));
    else
    uncert_comIN(looper,:,n)=zeros(1,size(comINvec,2));    
    end
end %END
looper = looper+1;
end

for m=1:ndim
comINvec(looper,:) = inputmatrix(:,m).^2;
for n=1:ndim %ADD
    if m==n
    uncert_comIN(looper,:,n)=2*inputmatrix(:,m);
    else
    uncert_comIN(looper,:,n)=zeros(1,size(comINvec,2));    
    end
end %END
looper = looper+1;
end

for m=1:ndim
comINvec(looper,:) = inputmatrix(:,m).*abs(inputmatrix(:,m));
for n=1:ndim %ADD
    if m==n
    uncert_comIN(looper,:,n)=abs(inputmatrix(:,m))+inputmatrix(:,m).*sign(inputmatrix(:,m));
    else
    uncert_comIN(looper,:,n)=zeros(1,size(comINvec,2));    
    end
end %END
looper = looper+1;
end

for k=1:ndim-1
  for m=k+1:ndim
comINvec(looper,:) = inputmatrix(:,k).*inputmatrix(:,m);
for n=1:ndim %ADD
    if k==n
    uncert_comIN(looper,:,n)=inputmatrix(:,m);
    else
    uncert_comIN(looper,:,n)=zeros(1,size(comINvec,2));    
    end
end %END
looper = looper+1;
  end
end

for k=1:ndim-1
  for m=k+1:ndim
comINvec(looper,:) = abs((inputmatrix(:,k)).*(inputmatrix(:,m)));
for n=1:ndim %ADD
    if k==n
    uncert_comIN(looper,:,n)=abs(inputmatrix(:,m)).*sign(inputmatrix(:,k));
    else
    uncert_comIN(looper,:,n)=zeros(1,size(comINvec,2));    
    end
end %END
looper = looper+1;
  end
end

for k=1:ndim-1
  for m=k+1:ndim
comINvec(looper,:) = inputmatrix(:,k).*abs(inputmatrix(:,m));
for n=1:ndim %ADD
    if k==n
    uncert_comIN(looper,:,n)=abs(inputmatrix(:,m));
    else
    uncert_comIN(looper,:,n)=zeros(1,size(comINvec,2));    
    end
end %END
looper = looper+1;
  end
end

for k=1:ndim-1
  for m=k+1:ndim
comINvec(looper,:) = abs(inputmatrix(:,k)).*(inputmatrix(:,m));
for n=1:ndim %ADD
    if k==n
    uncert_comIN(looper,:,n)=inputmatrix(:,m).*sign(inputmatrix(:,k));
    else
    uncert_comIN(looper,:,n)=zeros(1,size(comINvec,2));    
    end
end %END
looper = looper+1;
  end
end

  for m=1:ndim
comINvec(looper,:) = inputmatrix(:,m).^3;
for n=1:ndim %ADD
    if m==n
    uncert_comIN(looper,:,n)=3*inputmatrix(:,m).^2;
    else
    uncert_comIN(looper,:,n)=zeros(1,size(comINvec,2));    
    end
end %END
looper = looper+1;
  end

  for m=1:ndim
comINvec(looper,:) = abs(inputmatrix(:,m).^3);
for n=1:ndim %ADD
    if m==n
    uncert_comIN(looper,:,n)=sign(inputmatrix(:,m)).*(3*inputmatrix(:,m).^2);
    else
    uncert_comIN(looper,:,n)=zeros(1,size(comINvec,2));    
    end
end %END
looper = looper+1;
  end
 
%end of full term set
end

if themodel_FLAG == 2;
%start of truncated set

looper = 1;
for m=1:ndim
comINvec(looper,:) = inputmatrix(:,m);
for n=1:ndim %ADD
    if m==n
    uncert_comIN(looper,:,n)=ones(1,size(comINvec,2));
    else
    uncert_comIN(looper,:,n)=zeros(1,size(comINvec,2));    
    end
end %END
looper = looper+1;
end

for m=1:ndim
comINvec(looper,:) = inputmatrix(:,m).^2;
for n=1:ndim %ADD
    if m==n
    uncert_comIN(looper,:,n)=2*inputmatrix(:,m);
    else
    uncert_comIN(looper,:,n)=zeros(1,size(comINvec,2));    
    end
end %END
looper = looper+1;
end

for k=1:ndim-1
  for m=k+1:ndim
comINvec(looper,:) = inputmatrix(:,k).*inputmatrix(:,m);
for n=1:ndim %ADD
    if k==n
    uncert_comIN(looper,:,n)=inputmatrix(:,m);
    else
    uncert_comIN(looper,:,n)=zeros(1,size(comINvec,2));    
    end
end %END
looper = looper+1;
  end
end

%end of truncated term set
end

if themodel_FLAG == 3;
%start of linear set

looper = 1;
for m=1:ndim
comINvec(looper,:) = inputmatrix(:,m);
for n=1:ndim %ADD
    if m==n
    uncert_comIN(looper,:,n)=ones(1,size(comINvec,2));
    else
    uncert_comIN(looper,:,n)=zeros(1,size(comINvec,2));    
    end
end %END
looper = looper+1;
end

%end of linear term set
end

% end of the subroutine
end   

%
% Copyright ©2017 Andrew Meade.  All Rights Reserved.
%


