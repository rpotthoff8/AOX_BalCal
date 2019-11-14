function [rbfc_INminGZ]=place_GRBF(u,dainputs,wHist,cHist,center_daHist,h)
%Function places a single GRBF for the validation or approximation section.

%INPUTS:
%  u  =  Index of current GRBF
%  dainputscalib  =  Calibration (voltage - (global zeros) )
%  dainputs  = Current Section (voltage - (global zeros) )
%  centerIndexHist  =  History of GRBF centers, as determined from calibration section 
%  wHist  =  History of GRBF widths, as determined from calibration section 
%  cHist  =  History of GRBF coefficients, as determined from calibration section 
%  center_daHist = Voltages of GRBF centers. Dim 1= RBF #, Dim 2= Channel for voltage, Dim 3= Dimension center is placed in ( what load channel it is helping approximate)

%OUTPUTS:
%  rbfc_INminGZvalid  =  GRBF contribution for load approximation (aprxINminGZ2 = aprxINminGZ2+rbfc_INminGZ)

%Initialize variables
etavalid=zeros(size(dainputs));
rbfINminGZvalid=zeros(size(dainputs));
rbfc_INminGZ=zeros(size(dainputs));
dimFlag=size(center_daHist,2);

for s=1:size(dainputs,2) % loops through the components
    
    adiffervalid=center_daHist(u,:,s)-dainputs;
%     adiffervalid = dainputscalib(centerIndexHist(u,s),:)-dainputs;
    etavalid(:,s) = dot(adiffervalid,adiffervalid,2);
    
    rbfINminGZvalid(:,s)=((wHist(u,s)^dimFlag)/(sqrt(pi^dimFlag)))*exp(-((wHist(u,s)^2)*(etavalid(:,s)))/h^2); %From 'Iterated Approximate Moving Least Squares Approximation', Fasshauer and Zhang, Equation 22
    rbfINminGZvalid(:,s)=rbfINminGZvalid(:,s)-mean(rbfINminGZvalid(:,s),1);
    rbfc_INminGZ(:,s) = cHist(u,s)*rbfINminGZvalid(:,s);
end

end