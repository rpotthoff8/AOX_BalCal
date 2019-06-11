function [U_Amod,PI_test]=uncertainty_test(X,y,MSE,se,beta)

SEE=sqrt(MSE);
ub=se;
x_ave=mean(X);
y_ave=mean(y);
[n,k] = size(X);

y_est=y_ave+(X-x_ave)*beta; %Eqn 39, matches y_hat, confirmed for linear

umod_square=(SEE/sqrt(n))^2+((X-x_ave).^2)*(ub.^2); %Eqn 42

U_Amod=T_cr*sqrt(umod_square); %Eqn 44, basically matches y_hat_CI, confirmed for linear


for j = 1:n
    U_Amod_2(j,1)=T_cr*sqrt((sqrt(SSE/dof_e)/sqrt(n))^2+((X(j,:)-mean(X)).^2)*(sqrt(diag(SSE/dof_e*inv(X'*X))).^2)); %Eqn 44, basically matches y_hat_CI, confirmed for linear, truncated
    y_hat_CI_2(j,1) = T_cr*sqrt((SSE/dof_e)*X(j,:)*inv(X'*X)*X(j,:)');
end

% U_Amod(j,1)=T_cr*sqrt((SEE/sqrt(n))^2+((X-mean(X)).^2)*(ub.^2)); %Eqn 44, basically matches y_hat_CI, confirmed for linear, truncated

max_CI_dif=max(abs(U_Amod-y_hat_CI))
mean_CI_dif=mean(abs(U_Amod-y_hat_CI))

U_Adata=T_cr*SEE; %Eqn 46

PI_test=sqrt(U_Adata.^2+U_Amod.^2); %Basically matches y_hat_PI, confirmed for linear, truncated

PI_dif=abs(PI_test-y_hat_PI);
max_PI_dif=max(abs(PI_test-y_hat_PI))
mean_PI_dif=mean(abs(PI_test-y_hat_PI))

end