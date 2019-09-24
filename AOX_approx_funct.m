function [aprxINminGZapprox,loadPI_approx]=AOX_approx_funct(coeff,natzerosapprox,excessVecapprox,FLAGS,seriesapprox,series2approx,pointIDapprox,loadlist,output_location,GRBF,ANOVA,pct)
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
%  pct  =  Percent confidence level for PI

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
aprxINminGZapprox.ALG = aprxINapprox;

if FLAGS.loadPI==1
    loadPI_approx.ALG=calc_alg_PI(ANOVA,pct,comINapprox,aprxINapprox); %Calculate load PI
else
    loadPI_approx='PI NOT COMPUTED';
end

%OUTPUT
fprintf('\n ********************************************************************* \n');
if FLAGS.excel == 1
    %Output approximation load approximation
    filename = 'GLOBAL_ALG_APPROX.csv';
    approxinput=aprxINminGZapprox.ALG;
    description='APPROXIMATION ALGEBRAIC MODEL LOAD APPROXIMATION';
    print_approxcsv(filename,approxinput,description,pointIDapprox,seriesapprox,series2approx,loadlist,output_location);
else
    fprintf('\nAPPROXIMATION ALGEBRAIC MODEL LOAD APPROXIMATION RESULTS: Check aprxINminGZapprox in Workspace \n');
end

%OUTPUTING APPROXIMATION WITH PI
if FLAGS.approx_and_PI_print==1
    approxinput=cellstr(string(aprxINminGZapprox.ALG)+' +/- '+string(loadPI_approx.ALG));
    filename = 'APPROX_AOX_GLOBAL_ALG_RESULT_w_PI.csv';
    description='ALG APPROXIMATION LOAD APPROX WITH PREDICTION INTERVALS';
    print_approxcsv(filename,approxinput,description,pointIDapprox,seriesapprox,series2approx,loadlist,output_location);
end

%OUTPUTING PI VALUE
if FLAGS.PI_print==1
    filename = 'APPROX_ALG_PREDICTION_INTERVAL.csv';
    approxinput=loadPI_approx.ALG;
    description='APPROXIMATION ALGEBRAIC MODEL APPROXIMATION PREDICTION INTERVAL';
    print_approxcsv(filename,approxinput,description,pointIDapprox,seriesapprox,series2approx,loadlist,output_location);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    RBF SECTION FOR APPROXIMATION     AJM 6/29/17                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%goal to use centers, width and coefficients to approxate parameters against
%independent data

if FLAGS.balCal == 2
    
    numBasis=size(GRBF.epsHist,1);
    
    aprxINminGZ2approx = aprxINminGZapprox.ALG;
    aprxINminGZ_Histapprox = cell(numBasis,1);
    
    for u=1:numBasis
        
        %Call function to place single GRBF
        [rbfc_INminGZapprox]=place_GRBF(u,dainputsapprox,GRBF.epsHist,GRBF.cHist,GRBF.center_daHist,GRBF.h_GRBF);
        
        %update the approximation
        aprxINminGZ2approx = aprxINminGZ2approx+rbfc_INminGZapprox;
        aprxINminGZ_Histapprox{u} = aprxINminGZ2approx;
        
    end
    
    %OUTPUT
    aprxINminGZapprox.GRBF=aprxINminGZ2approx;
    
    
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