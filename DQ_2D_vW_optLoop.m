function [xc, w, zstar,data] = DQ_2D_vW_optLoop(data,Res,DQcoeff,options)
method='wahba 5 Differential Quadrature, all derivative directions, different width, bounded search';
x=data.x;

if isfield(data,'dist')
    dist=data.dist;
else
    dist=zeros(size(x,1),size(x,1),size(x,2));
    for i=1:size(x,2)
        dist(:,:,i)=x(:,i)'-x(:,i); %solve distance in each dimension, Eqn 16 from Javier's notes
    end
    data.dist=dist;
end

if isfield(data,'dist_square')
    dist_square=data.dist_square;
else
    dist_square=dist.^2; %squared distance in each dimension;
    dist_square(dist_square==0)=NaN; %eliminate zero values (on diagonal)
    data.dist_square=dist_square;
end

if isfield(data,'min_dist_square')
    min_dist_square=data.min_dist_square;
else
    min_dist_square=min(dist_square,[],1);
    data.min_dist_square=min_dist_square;
end

% 
% spacing=x(2:size(x,1),:)-x(1:size(x,1)-1,:);
% spacing(spacing<=0)=NaN;
% minSpace_dim=min(spacing);
% minSpace=(minSpace_dim);

n = size(x,1);
err = 100.0;
[rmax,zstar] = max(abs(Res(:))); %Find and store location and value of maximum absolute value residual

%Start Differential Quadrature for 2nd Derivative
dim=size(x,2); %Dimension to take the 2nd derivative in
approx2D=zeros(size(x,2),size(x,1));
for i=1:dim
    approx2D(i,:)=sum(Res.*DQcoeff(:,:,i)); %2nd Derivative wrt R Approx: Eqn 9 from 'Development of RBF-DQ method... Y.L Wu
end
%End differential quadrature
dubya=(-1/2)*(1/rmax)*approx2D;

% dubya(dubya>0)=0; %Widths must be negative so set positive widths to zero
dubya=-abs(dubya); %Widths must be negative
%Test each xc and solved width
err_i=zeros(size(x,1),1);
for i = 1:n
    xc_i = x(i,:)';
    w_i = dubya(:,i);
    err_i(i) = meritFunction(x,Res,xc_i,w_i);
end

%Select location that gives lowest error
[err,i_best]=min(err_i);
w_guess=dubya(:,i_best);
xc = x(i_best,:)';

if isfield(options,'h')
    h=options.h;
else
    h = 0.25; 
end
wmin = (log(h)./(min_dist_square(:,i_best,:)));
wmin=reshape(wmin,size(w_guess));

 [w, err_k] = fminsearchbnd(@(w_best) meritFunction(x,Res,xc,w_best),w_guess,wmin,zeros(size(wmin))); %Multiple optimization variable: bounded but with different width in every dimension

end