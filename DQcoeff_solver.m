function [DQcoeff]=DQcoeff_solver(x)

%Start Differential Quadrature for 2nd Derivative
dist=zeros(size(x,1),size(x,1),size(x,2));
for i=1:size(x,2)
    dist(:,:,i)=x(:,i)'-x(:,i); %solve distance in each dimension, Eqn 16 from Javier's notes
end
R_square=sum(dist.^2,3); %Eqn 17 from Javier's notes
R=sqrt(R_square);
%Solving for derivative approx coefficient values
c=25; %Shape parameter: Figure 3 from 'Development of RBF-DQ method... Y.L Wu
B=exp(-c*R_square'); %Eqn 10 from 'Development of RBF-DQ method... Y.L Wu
%Solving 2nd Derivative wrt 'derdim'
DQcoeff=zeros(size(x,1),size(x,1),size(x,2));
for i=1:size(x,2)
phiR2D=exp(-c*R_square).*((-2*c*dist(:,:,i)).^2)-exp(-c*R_square).*2*c;
DQcoeff(:,:,i)=B\phiR2D;
% DQcoeff(:,:,i)=lsqminnorm(B,phiR2D);
end
%End differential quadrature
end