function [combined_uncert,tare_uncert, FL_uncert,coeff_uncert]=uncert_prop_anova(xcalib,beta_CI,comIN,dimFlag,uncert_comIN,s_1st0,nterms,targetMatrix,series,voltTrust,Anova_Flag,Volt_Flag)
%Function calculates uncertainty in load from uncertainty in coefficients and input voltages

%Inputs:
%xcalib: Coefficients from calibration
%beta_ci: Confidence interval for coefficients
%comIN: Produced by balCal_algEquations3: Voltages provided put into matrix
%to multiply by coefficients for load output
%dimFlag: Number of dimensions (channels)
%uncert_comIN: Produced by balCal_algEquations3: Voltages provided put into
%matrix to multiply by coefficients for partial derivatives
%s_1st0: Index of where which datapoints are tare loads
%nterms: number of coefficients used for model
%targetMatrix: desired output loads
%voltTrust: 95% confidence (+/-) for voltage readings (in mV)
%Boot_Flag: 1 or 0: flag for if anova method was used to find coeff CI
%Volt_Flag: 1 or 0: If voltage uncertainty prop is to be performed

%Outputs:
%combined_uncert: Uncertainty (in pounds) due to coefficient and voltage uncertainty
%for each datapoint
%tare_uncert: Uncertainty (in pounds) due to coefficient and voltage uncertainty
%for tare loads only
%FL_uncert: Uncertainty for loads of (calculated-tare): combines
%uncertainty in calculated loads and tare loads

%START: coeff uncertainty propagation JRP 16 Jan 19
%Take Confidence interval found for coefficients using bootstrap, and store
%the absolute value of the larger error between +/- in a matrix
if Anova_Flag==1
        xcalib_error=beta_CI;
else
    xcalib_error=zeros(size(xcalib,1),size(xcalib,2));
end
comIN_square=comIN.^2; %Square input matrix
xcalib_error_square=xcalib_error.^2; %square error matrix
coeff_uncert_square=(comIN_square*xcalib_error_square); %Matrix for error in each calibration point from uncertainty in coefficients(+/-)
coeff_uncert=coeff_uncert_square.^.5;
%END:  coeff uncertainty propagation

%ADDED 23 Jan 19 JRP: Analytical calc of uncert in loads due to uncertainty in read
%voltages
if Volt_Flag==1
for i=1:max(series)
    for j=1:dimFlag
        uncert_comIN(nterms+i,:,j) = 0;
    end
end
uncert_comIN_use=uncert_comIN(1:size(xcalib,1),:,:);
volt_uncert_square=zeros(size(targetMatrix));
for i=1:dimFlag
    partial(:,:,i)=uncert_comIN_use(:,:,i)'*xcalib;
    uncert_channel_square(:,:,i)=((partial(:,:,i)).^2).*voltTrust^2;
    volt_uncert_square=volt_uncert_square+(uncert_channel_square(:,:,i));
end
else
    volt_uncert_square=zeros(size(targetMatrix));
end
%Combine uncertainty from coeff and input voltages for 1 total uncertainty
%value for every datapoint:
combined_uncert=(coeff_uncert_square+volt_uncert_square).^(.5);

tare_uncert=beta_CI(nterms+1:size(beta_CI,1),:); %Uncert for tare loads

%FL Uncert:
for i=1:numel(series)
    if i==s_1st0(series(i))
        FL_uncert(i,:)=combined_uncert(i,:);
    else
        FL_uncert(i,:)=(combined_uncert(i,:).^2+tare_uncert(series(i)).^2).^.5;
    end

end
%END added for uncert
end