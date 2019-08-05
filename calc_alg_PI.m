function [loadPI]=calc_alg_PI(ANOVA,pct,comIN,aprxIN)
%Function computes the load prediction interval from the ANOVA results and
%user specified % confidence level

%INPUTS:
%  ANOVA = Structure of ANOVA Results
%  pct  =  Percent Confidence level for PI width
%  comIN  =  Input voltages/combined voltages
%  aprxIN  =  Approximation of loads

%OUTPUTS:
%  loadPI = Load Prediction Interval (+/- in load units)

loadPI=zeros(size(aprxIN,1),size(aprxIN,2));
for i=1:size(aprxIN,2)
    T_cr = tinv(1-((1-(pct/100))/2),ANOVA(i).PI.dof_e); %T stat based on confidence level set
    for j = 1:size(comIN,1)
        loadPI(j,i)=T_cr*sqrt(ANOVA(i).PI.sigma_hat_sq*(1+(comIN(j,:)*ANOVA(i).PI.invXtX*comIN(j,:)')));
    end
end
end