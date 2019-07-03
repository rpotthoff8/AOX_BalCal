function [rbfc_INminGZ]=place_GRBF(u,dainputscalib,dainputs,centerIndexHist,wHist,cHist)
%Function places a single GRBF for the validation or approximation section.

%INPUTS:
%  u  =  Index of current GRBF
%  dainputscalib  =  Calibration (voltage - (global zeros) )
%  dainputs  = Current Section (voltage - (global zeros) )
%  centerIndexHist  =  History of GRBF centers, as determined from calibration section 
%  wHist  =  History of GRBF widths, as determined from calibration section 
%  cHist  =  History of GRBF coefficients, as determined from calibration section 

%OUTPUTS:
%  rbfc_INminGZvalid  =  GRBF contribution for load approximation (aprxINminGZ2 = aprxINminGZ2+rbfc_INminGZ)

%Initialize variables
etavalid=zeros(size(dainputs));
rbfINminGZvalid=zeros(size(dainputs));
rbfc_INminGZ=zeros(size(dainputs));

for s=1:size(dainputs,2) % loops through the components
    
    adiffervalid = dainputscalib(centerIndexHist(u,s),:)-dainputs;
    etavalid(:,s) = dot(adiffervalid,adiffervalid,2);
    
    rbfINminGZvalid(:,s)=exp(etavalid(:,s)*log(abs(wHist(u,s))));
    
    rbfc_INminGZ(:,s) = cHist(u,s)*rbfINminGZvalid(:,s);
end

end