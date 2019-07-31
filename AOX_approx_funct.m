function []=AOX_approx_funct(coeff,natzerosapprox,excessVecapprox,FLAGS,seriesapprox,series2approx,pointIDapprox,loadlist,output_location,GRBF,ANOVA)
%Function performs all calculations and outputs for approximation.  This
%function can be called by both the main AOX_balcal code or the standalone
%approximation code

%INPUTS:
%  coeff = Calibration Coefficients
%  natzerosapprox  =  Natural zero voltages
%  excessVecapprox  =  Approximation voltages
%  FLAGS  =  User option flags
%  seriesapprox  =  Series labels
%  series2approx  =  Series2 labels
%  pointIDapprox  =  data point IDs 
%  loadlist  =  channel load labels
%  output_location  =  save location for output files
%  GRBF  =  Structure containing GRBF centers, widths, and coefficients if RBFs were placed in calibration
%  ANOVA  =  ANOVA results needed for calculating PI

%OUTPUTS:


[~, dimFlag] = size(excessVecapprox);

%natural zeros (also called global zeros)
globalZerosapprox = mean(natzerosapprox,1);

% Subtract the Global Zeros from the Inputs
dainputsapprox = excessVecapprox-globalZerosapprox;

% Call the Algebraic Subroutine
comINapprox = balCal_algEqns(FLAGS.model,dainputsapprox,seriesapprox,0);

%LOAD APPROXIMATION
%define the approximation for inputs minus global zeros
aprxINapprox = comINapprox*coeff;        %to find approximation AJM111516
aprxINminGZapprox = aprxINapprox;

if FLAGS.loadPI==1
    loadPI_approx=zeros(size(aprxINapprox,1),size(aprxINapprox,2));
    for i=1:dimFlag
        for j = 1:size(aprxINapprox,1)
            loadPI_approx(j,i)=ANOVA(i).PI.T_cr*sqrt(ANOVA(i).PI.sigma_hat_sq*(1+(comINapprox(j,:)*ANOVA(i).PI.invXtX*comINapprox(j,:)')));
        end
    end
end

%OUTPUT
fprintf('\n ********************************************************************* \n');
if FLAGS.excel == 1
    %Output approximation load approximation
    filename = 'GLOBAL_ALG_APPROX.csv';
    approxinput=aprxINminGZapprox;
    description='APPROXIMATION ALGEBRAIC MODEL LOAD APPROXIMATION';
    print_approxcsv(filename,approxinput,description,pointIDapprox,seriesapprox,series2approx,loadlist,output_location);
else
    fprintf('\nAPPROXIMATION ALGEBRAIC MODEL LOAD APPROXIMATION RESULTS: Check aprxINminGZapprox in Workspace \n');
end

%OUTPUTING APPROXIMATION WITH PI
if FLAGS.approx_and_PI_print==1
    approxinput=cellstr(string(aprxINminGZapprox)+' +/- '+string(loadPI_approx));
    filename = 'APPROX_AOX_GLOBAL_ALG_RESULT_w_PI.csv';
    description='ALG APPROXIMATION LOAD APPROX WITH PREDICTION INTERVALS';
    print_approxcsv(filename,approxinput,description,pointIDapprox,seriesapprox,series2approx,loadlist,output_location);
end

%OUTPUTING PI VALUE
if FLAGS.PI_print==1
    filename = 'APPROX_ALG_PREDICTION_INTERVAL.csv';
    approxinput=loadPI_approx;
    description='APPROXIMATION ALGEBRAIC MODEL APPROXIMATION PREDICTION INTERVAL';
    print_approxcsv(filename,approxinput,description,pointIDapprox,seriesapprox,series2approx,loadlist,output_location);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    RBF SECTION FOR APPROXIMATION     AJM 6/29/17                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%goal to use centers, width and coefficients to approxate parameters against
%independent data

if FLAGS.balCal == 2
    
    numBasis=size(GRBF.wHist,1);
    
    aprxINminGZ2approx = aprxINminGZapprox;
    aprxINminGZ_Histapprox = cell(numBasis,1);
    
    for u=1:numBasis
        
        %Call function to place single GRBF
        [rbfc_INminGZapprox]=place_GRBF(u,dainputsapprox,GRBF.wHist,GRBF.cHist,GRBF.center_daHist);
        
        %update the approximation
        aprxINminGZ2approx = aprxINminGZ2approx+rbfc_INminGZapprox;
        aprxINminGZ_Histapprox{u} = aprxINminGZ2approx;
        
    end
    
    %OUTPUT
    fprintf('\n ********************************************************************* \n');
    if FLAGS.excel == 1
        %Output approximation load approximation
        filename = 'GLOBAL_ALG+GRBF_APPROX.csv';
        approxinput=aprxINminGZ2approx;
        description='APPROXIMATION ALGEBRAIC+GRBF MODEL LOAD APPROXIMATION';
        print_approxcsv(filename,approxinput,description,pointIDapprox,seriesapprox,series2approx,loadlist,output_location);
    else
        fprintf('\nAPPROXIMATION ALGEBRAIC+GRBF MODEL LOAD APPROXIMATION RESULTS: Check aprxINminGZapprox in Workspace \n');
    end
    
end
% END APPROXIMATION GRBF SECTION