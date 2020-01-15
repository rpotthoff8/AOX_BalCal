function [loadPI]=calc_PI(ANOVA,pct,comIN,aprxIN,calc_channel)
%Function computes the load prediction interval from the ANOVA results and
%user specified % confidence level

%INPUTS:
%  ANOVA = Structure of ANOVA Results
%  pct  =  Percent Confidence level for PI width
%  comIN  =  Input voltages/combined voltages
%  aprxIN  =  Approximation of loads
%  calc_channel  =  Logical array for if PI should be calculated in channel

%OUTPUTS:
%  loadPI = Load Prediction Interval (+/- in load units)

dimFlag=size(aprxIN,2); %Determine dimensions of data

if exist('calc_channel','var')==0 %If no variable is provided for which channels should be calculated
    calc_channel=ones(1,dimFlag); %Calculate PI in all channels
end

%Initialize variables
loadPI=zeros(size(aprxIN,1),size(aprxIN,2)); 
T_cr=zeros(1,dimFlag);
sigma_hat_sq=zeros(1,dimFlag);

for i=1:dimFlag %Calculate PI in each channel
    if calc_channel(i)==1
        T_cr(i) = tinv(1-((1-(pct/100))/2),ANOVA(i).PI.dof_e); %T stat based on confidence level set
        sigma_hat_sq(i)=ANOVA(i).PI.sigma_hat_sq; %Load sigma_hat_sq from ANOVA
%         for j = 1:size(comIN,1)
%             
%             loadPI(j,i)=T_cr(i)*sqrt(ANOVA(i).PI.sigma_hat_sq*(1+(comIN(j,:)*ANOVA(i).PI.invXtX*comIN(j,:)')));
%         end
        loadPI(:,i)=diag(T_cr(i)*sqrt(sigma_hat_sq(i)*(1+(comIN(:,:)*ANOVA(i).PI.invXtX*comIN(:,:)')))); %Calculate PI for each datapoint in channel 'i'
        
    end
end
end