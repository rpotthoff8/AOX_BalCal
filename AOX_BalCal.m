%balanceCalibration_with_RBF_8D.m
%requires "balCal_meritFunction.m" to run
%input file: "BuffetBalance-CalDataOfOct2015-3F3M.csv"
%output file: "balCal_output_ALGB.xls"
%output file: "balCal_output_GRBF.xls"
%%
%initialize the workspace
clc;
clearvars;
close all;
% workspace;
fprintf('Copyright 2019 Andrew Meade, Ali Arya Mokhtarzadeh, Javier Villarreal, and John Potthoff.  All Rights Reserved.\n')
% Because of measurement noise in the voltage the APPROXIMATION tare is computed
% by post-processing. The average and stddev is taken of all channels per section.
% If the stddev is less than 0.25% of the capacity for any station the tare
% is equal to the average for that channel. If the stddev is greater than 0.25%
% then the approximation at the local zero is taken as the tare for that channel
% and section. The global approximation is left alone but the tare is subtracted
% from the values to make the local loads.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       USER INPUT SECTION
out = AOX_GUI;
if out.cancel == 1
    return
end
%TO SELECT Algebraic Model                                  set FLAGS.balCal = 1;
%TO SELECT Algebraic and GRBF Model                         set FLAGS.balCal = 2;
FLAGS.balCal = out.grbf;
%DEFINE THE NUMBER OF BASIS FUNCTIONS
numBasis = out.basis;
%
%SELECT ALGEBRAIC MODE                                      set FLAGS.model = 1 (full)
%                                                                             2 (trunc)
%                                                                             3 (linear)
%                                                                             4 (custom)
FLAGS.model = out.model;
%
%TO PRINT LOAD PERFORMANCE PARAMETERS TO CSV                set FLAGS.print = 1;
FLAGS.print = out.print;
%
%TO DISPLAY LOAD PERFORMANCE PARAMETERS IN COMMAND WINDOW   set FLAGS.disp = 1;
FLAGS.disp= out.disp;
%
%TO SAVE DATA TO CSV                                        set FLAGS.excel = 1;
FLAGS.excel = out.excel;
%
%TO PRINT INPUT/OUTPUT CORRELATION PLOTS                    set FLAGS.corr = 1;
FLAGS.corr = out.corr;
%
%TO PRINT INPUT/RESIDUALS CORRELATION PLOTS                 set FLAGS.rescorr = 1;
FLAGS.rescorr = out.rescorr;
%
%TO PRINT ORDER/RESIDUALS PLOTS                             set rest_FLAG = 1;
FLAGS.res = out.res;
%
%TO PRINT RESIDUAL HISTOGRAMS                               set FLAGS.hist = 1;
FLAGS.hist = out.hist;
%
%TO SELECT Validation of the Model                          set FLAGS.balVal = 1;
FLAGS.balVal = out.valid;
%
%TO SELECT Approximation from Cal Data                      set FLAGS.balApprox = 1;
FLAGS.balApprox = out.approx;
%
%TO FLAG POTENTIAL OUTLIERS                                 set FLAGS.balOut = 1;
FLAGS.balOut = out.outlier;
numSTD = out.numSTD;  %Number of St.D. for outlier threshold.
%
%TO REMOVE POTENTIAL OUTLIERS                               set FLAGS.zeroed = 1;
FLAGS.zeroed = out.zeroed;
%
%Uncertainty button outputs
FLAGS.volt=out.voltFlag;
voltTrust=out.voltTrust;

FLAGS.anova = out.anova;
FLAGS.loadPI = out.loadPI;
FLAGS.BALFIT_Matrix=out.BALFIT_Matrix;
FLAGS.BALFIT_ANOVA=out.BALFIT_ANOVA;
FLAGS.Rec_Model=out.Rec_Model;
anova_pct=out.anova_pct;
FLAGS.approx_and_PI_print=out.approx_and_PI_print;
FLAGS.custom_eqn_iter=out.stableRec_FLAGcheck;

REPORT_NO=out.REPORT_NO;
output_location=out.output_location;

%TO SAVE .MAT FILE OF CALIBRATION MODEL
FLAGS.calib_model_save=out.calib_model_save_FLAG;
%TO SAVE INPUT .CAL, .VAL, .APP FILE IN OUTPUT LOCATION
FLAGS.input_save=out.input_save_FLAG;
%                       END USER INPUT SECTION
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       INITIALIZATION SECTION
fprintf('\nWorking ...\n')
% Load data and characterize series
load(out.savePathcal,'-mat');
series0 = series;
series20=series2;
pointID0=pointID;
[~,s_1st0,~] = unique(series0);
nseries0 = length(s_1st0);
[numpts0, dimFlag] = size(excessVec0);
% Loads:
% loadlabes, voltlabels (if they exist)
% loadCapacities, natzeros, targetMatrix0, excessVec0, series0

% Load the custom equation matrix if using a custom algebraic model
% SEE: CustomEquationMatrixTemplate.csv
if FLAGS.model == 4
    customMatrix = out.customMatrix;
    customMatrix = [customMatrix; ones(nseries0,dimFlag)];
else
    customMatrix = 1;
end

% Load data labels if present, otherwise use default values.
if exist('loadlabels','var')
    loadlist = loadlabels;
    voltagelist = voltlabels;
    reslist = strcat('res',loadlist);
else
    loadlist = {'NF','BM','S1','S2','RM','AF','PLM', 'PCM', 'MLM', 'MCM'};
    voltagelist = {'rNF','rBM','rS1','rS2','rRM','rAF','rPLM','rPCM','rMLM','rMCM'};
    reslist = strcat('res',loadlist);
end

if exist('loadunits','var')==0
    loadunits = {'lbs','in-lbs','lbs','lbs','in-lbs','lbs','in-lbs', 'in-lbs', 'in-lbs', 'in-lbs'};
    voltunits = {'microV/V','microV/V','microV/V','microV/V','microV/V','microV/V','microV/V','microV/V','microV/V','microV/V'};
end

% Prints output vs. input and calculates correlations
if FLAGS.corr == 1
    figure('Name','Correlation plot','NumberTitle','off','WindowState','maximized');
    correlationPlot(targetMatrix0, excessVec0, loadlist, voltagelist);
end

%                       END INITIALIZATION SECTION
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  CALIBRATION - ALGEBRAIC SECTION                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize structure for unique outputs for section
uniqueOut=struct();

% Finds the average  of the natural zeros (called global zeros)
globalZeros = mean(natzeros,1);

% Subtracts global zeros from signal.
dainputs0 = excessVec0 - globalZeros;

% Determines how many terms are in the algebraic model; this will help
% determine the size of the calibration matrix
switch FLAGS.model
    case {1,4}
        % Full Algebraic Model or Custom Algebraic Model
        % The Custom model calculates all terms, and then excludes them in
        % the calibration process as determined by the customMatrix.
        nterms = 2*dimFlag*(dimFlag+2);
    case 2
        % Truncated Algebraic Model
        nterms = dimFlag*(dimFlag+3)/2;
    case 3
        % Linear Algebraic Model
        nterms = dimFlag;
end

% Creates the algebraic combination terms of the inputs.
% Also creates intercept terms; a different intercept for each series.
comIN0 = balCal_algEqns(FLAGS.model,dainputs0,series0,1);

%%% Balfit Stats and Regression Coeff Matrix
balfitdainputs0 = targetMatrix0;
balfittargetMatrix0 = balCal_algEqns(3,dainputs0,series0,0);
balfitcomIN0 = balCal_algEqns(FLAGS.model,balfitdainputs0,series0,1);
%%% Balfit Stats and Regression Coeff Matrix

fprintf('\nStarting Calibration Calculations\n')

%Creates vectors that will not have outliers removed for balfit
series = series0;
targetMatrix = targetMatrix0;
comIN = comIN0;

%Calculate xcalib (coefficients)
[xcalib, ANOVA] = calc_xcalib(comIN       ,targetMatrix       ,series,...
    nterms,nseries0,dimFlag,FLAGS,customMatrix,anova_pct,loadlist,'Direct');

[balfitxcalib, balfitANOVA] = calc_xcalib(balfitcomIN0,balfittargetMatrix0,series,...
    nterms,nseries0,dimFlag,FLAGS,customMatrix,anova_pct,voltagelist,'BALFIT');

% APPROXIMATION
% define the approximation for inputs minus global zeros (includes
% intercept terms)
aprxIN = comIN0*xcalib;

% RESIDUAL
targetRes = targetMatrix0-aprxIN;

% Identify Outliers After Filtering
% (Threshold approach)
if FLAGS.balOut == 1

    %Identify outliers based on residuals
    fprintf('\n Identifying Outliers....')

    [OUTLIER_ROWS,num_outliers,prcnt_outliers,rowOut,colOut] = ID_outliers(targetRes,loadCapacities,numpts0,dimFlag,numSTD,FLAGS);

    newStruct = struct('num_outliers',num_outliers,...
        'prcnt_outliers',prcnt_outliers,...
        'rowOut',rowOut,...
        'colOut',colOut,...
        'numSTD',numSTD);
    uniqueOut = cell2struct([struct2cell(uniqueOut);struct2cell(newStruct)],...
        [fieldnames(uniqueOut); fieldnames(newStruct)],1);

    fprintf('Complete\n')

    % Use the reduced input and target files
    if FLAGS.zeroed == 1
        fprintf('\n Removing Outliers....')
        % Remove outlier rows for recalculation and all future calculations:
        numpts0 =  numpts0 - num_outliers;
        targetMatrix0(OUTLIER_ROWS,:) = [];
        excessVec0(OUTLIER_ROWS,:) = [];
        series0(OUTLIER_ROWS) = [];
        series20(OUTLIER_ROWS)=[];
        pointID0(OUTLIER_ROWS)=[];
        comIN0(OUTLIER_ROWS,:) = [];
        [~,s_1st0,~] = unique(series0);
        nseries0 = length(s_1st0);
        fprintf('Complete\n')

        %Calculate xcalib (coefficients)
        [xcalib,ANOVA] = calc_xcalib(comIN0,targetMatrix0,series0,...
            nterms,nseries0,dimFlag,FLAGS,customMatrix,anova_pct,loadlist,'Direct');

        %%% Balfit Stats and Regression Coeff Matrix
        [balfitxcalib, balfitANOVA] = calc_xcalib(balfitcomIN0,balfittargetMatrix0,series,...
            nterms,nseries0,dimFlag,FLAGS,customMatrix,anova_pct,voltagelist,'BALFIT');
        %%% Balfit Stats and Matrix

        % APPROXIMATION
        % define the approximation for inputs minus global zeros (includes
        % intercept terms)
        aprxIN = comIN0*xcalib;

        % RESIDUAL
        targetRes = targetMatrix0-aprxIN;
    end
end

%Iterate to find stable recommended equation
if FLAGS.custom_eqn_iter==1
    FLAGS_iter.anova=1;
    FLAGS_iter.model=4;
    FLAGS_iter.test_FLAG=1;
    customMatrix_iter=[zeros(nterms,dimFlag);ones(nseries0,dimFlag)];
    customMatrix_last=customMatrix_iter;
    ANOVA_iter=ANOVA;
    for i=1:dimFlag
        customMatrix_iter(1:nterms,i)=ANOVA_iter(i).sig(1:nterms);
    end

    samRec=0;
    iter_count=0;
    while samRec==0
        iter_count=iter_count+1;
        fprintf('\n Searching for stable recommended equation. Iteration '); fprintf(string(iter_count)); fprintf('\n');
        customMatrix_last=customMatrix_iter;

        [~,ANOVA_iter] = calc_xcalib(comIN0,targetMatrix0,series0,...
            nterms,nseries0,dimFlag,FLAGS_iter,customMatrix_iter,anova_pct,loadlist,'Direct');
        for i=1:dimFlag
            customMatrix_iter(1:nterms,i)=ANOVA_iter(i).sig(1:nterms);
        end

        samRec=isequal(customMatrix_last,customMatrix_iter);


    end
    RECOMM_ALG_EQN_STABLE=customMatrix_iter(1:nterms,:);
    newStruct = struct('RECOMM_ALG_EQN_STABLE',RECOMM_ALG_EQN_STABLE);
    uniqueOut = cell2struct([struct2cell(uniqueOut);struct2cell(newStruct)],...
        [fieldnames(uniqueOut); fieldnames(newStruct)],1);
end



% Splits xcalib into Coefficients and Intercepts (which are negative Tares)
coeff = xcalib(1:nterms,:);
tares = -xcalib(nterms+1:end,:);
intercepts=-tares;
taretal=tares(series0,:);
aprxINminGZ=aprxIN+taretal; %Approximation that does not include intercept terms

%    QUESTION: JRP; IS THIS NECESSARY/USEFUL?
[~,tares_STDDEV_all] = meantare(series0,aprxINminGZ-targetMatrix0);
tares_STDDEV = tares_STDDEV_all(s_1st0,:);

%%% Balfit Stats and Regression Coeff Matrix
balfit_C1INV = xcalib((1:dimFlag), :);
balfit_D1 = zeros(dimFlag,dimFlag);
balfit_INTERCEPT = globalZeros;
balfit_C1INVC2 = balfitxcalib((dimFlag+1:nterms), :)*balfit_C1INV;
balfit_regress_matrix = [globalZeros ; balfit_INTERCEPT ; balfit_C1INV ; balfit_D1 ; balfit_C1INVC2 ];
%%% Balfit Stats and Matrix

%Start uncertainty section
%ANOVA data for uncertainty
beta_CI_comb=zeros(size(xcalib,1),dimFlag);
y_hat_PI_comb=zeros(size(targetMatrix0));
if FLAGS.anova==1
    for j=1:dimFlag
        beta_CI_comb(:,j)=ANOVA(j).beta_CI;
        y_hat_PI_comb(:,j)=ANOVA(j).y_hat_PI;
    end
end
%END: ANOVA data for uncertainty

if FLAGS.volt==1
    %uncertainty due to uncertainty in volt readings
    uncert_comIN=balCal_algEquations_partialdiff(FLAGS.model, dimFlag, dainputs0);
else
    uncert_comIN=zeros(nterms,numpts0,dimFlag);
end

[combined_uncert_anova,tare_uncert_anova, FL_uncert_anova,coeff_uncert_anova]=...
    uncert_prop_anova(xcalib,beta_CI_comb,comIN0,dimFlag,uncert_comIN,s_1st0,nterms,targetMatrix0,series0,voltTrust,FLAGS.anova,FLAGS.volt);
%end uncertainty section

%OUTPUT FUNCTION
%Function creates all outputs for calibration, algebraic section
section={'Calibration Algebraic'};

newStruct=struct('aprxIN',aprxIN,...
    'coeff',coeff,...
    'nterms',nterms,...
    'ANOVA',ANOVA,...
    'balfitcomIN',balfitcomIN0,...
    'balfitxcalib',balfitxcalib,...
    'balfittargetMatrix',balfittargetMatrix0,...
    'balfitANOVA',balfitANOVA,...
    'balfit_regress_matrix',balfit_regress_matrix,...
    'targetMatrix0',targetMatrix0,...
    'loadunits',{loadunits(:)},...
    'voltunits',{voltunits(:)},...
    'balance_type',balance_type,...
    'description',description);
uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
    [fieldnames(uniqueOut); fieldnames(newStruct)],1);

output(section,FLAGS,targetRes,loadCapacities,fileName,numpts0,nseries0,...
    tares,tares_STDDEV,loadlist,series0,excessVec0,dimFlag,voltagelist,...
    reslist,numBasis,pointID0,series20,output_location,REPORT_NO,uniqueOut)

%END CALIBRATION ALGEBRAIC SECTION
%%
if FLAGS.balCal == 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       CALIBRATION - RBF SECTION                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %goal to minimize: minimize the sum of the squares (dot product) of each of the 8
    %residual vectors 'targetRes' 'target1' ... 'target8'
    %dt1 = dot(target1,target1);
    %find centers by finding the index of max residual, using that index to
    %subtract excess(counter)-excess(indexMaxResid) and then taking the dot
    %product of the resulting column vector

    %Initialize structure for unique outputs for section
    uniqueOut=struct();

    targetRes2=targetRes;
    aprxINminGZ2 = aprxINminGZ;
    dainputscalib = excessVec0-globalZeros;

    %Initialize Variables
    aprxINminGZ_Hist = cell(numBasis,1);
    tareGRBFHist = cell(numBasis,1);
    centerLoop=zeros(dimFlag,dimFlag);
    centerIndexLoop_sub=zeros(dimFlag,dimFlag);    
    eta=zeros(length(excessVec0(:,1)),dimFlag);
    eps=zeros(1,dimFlag);
    rbfINminGZ=zeros(length(excessVec0(:,1)),dimFlag);
    coeffRBF=zeros(1,dimFlag);
    rbfc_INminGZ=zeros(length(excessVec0(:,1)),dimFlag);
    epsHist=zeros(numBasis,dimFlag);
    cHist=zeros(numBasis,dimFlag);
    centerIndexHist=zeros(numBasis,dimFlag,dimFlag);
    center_daHist=zeros(numBasis,dimFlag,dimFlag);
    resSquareHist=zeros(numBasis,dimFlag);
    resStdHist=zeros(numBasis,dimFlag);

    dist=zeros(size(dainputscalib,1),size(dainputscalib,1),size(dainputscalib,2));
    for i=1:size(dainputscalib,2)
        dist(:,:,i)=dainputscalib(:,i)'-dainputscalib(:,i); %solve distance in each dimension, Eqn 16 from Javier's notes
    end
    R_square=sum(dist.^2,3); %Eqn 17 from Javier's notes: squared distance between each point
    R_square_find=R_square;
    R_square_find(R_square_find==0)=NaN; %Eliminate zero values (on diagonal)
    min_R_square=min(R_square_find); %Find distance to closest point
    %Set limits on width (shape factor)
    h_GRBF=sqrt(max(min(R_square_find)));
    eps_min=0.1; %Fasshauer pg 234
    eps_max=1.0;
    center_min=min(dainputscalib);
    center_max=max(dainputscalib);
    
    
    for u=1:numBasis
        for s=1:dimFlag
            [~,maxI]=max(targetRes2(:,s));
            
            %find epsilon and center index via ga
            fmin_options=optimset('Display','off');
            ga_options=optimoptions('ga','Display','off');
            ga_merit=@(x) balCal_meritFunction2(x(1),x(2:dimFlag+1),dainputscalib,h_GRBF,dimFlag,targetRes2(:,s));
%             ga_out=ga(ga_merit,1+dimFlag,[],[],[],[],[eps_min; center_min'],[eps_max;center_max'],[],2,ga_options);
            ga_out=fminsearchbnd(ga_merit,[0.5;dainputscalib(maxI,:)'],[eps_min; center_min'],[eps_max;center_max'],fmin_options);
            eps(s)=ga_out(1);
            centerLoop(:,s)=ga_out(2:dimFlag+1);
            
            dist=centerLoop(:,s)'-dainputscalib;
            R_square=sum(dist.^2,2); %Eqn 17 from Javier's notes: squared distance between each point
            
            rbfINminGZ(:,s)=((eps(s)^dimFlag)/(sqrt(pi^dimFlag)))*exp(-((eps(s)^2)*(R_square))/h_GRBF^2); %From 'Iterated Approximate Moving Least Squares Approximation', Fasshauer and Zhang, Equation 22
            coeffRBF(s) = lsqminnorm(rbfINminGZ(:,s),targetRes2(:,s));
            
            rbfc_INminGZ(:,s) = coeffRBF(s)*rbfINminGZ(:,s);
        end

        %Store basis parameters in Hist variables
        epsHist(u,:) = eps;
        cHist(u,:) = coeffRBF;
        for s=1:dimFlag
            centerIndexHist(u,:,s) = centerLoop(:,s);
            center_daHist(u,:,s)=centerLoop(:,s); %Variable stores the voltages of the RBF centers.
            %Dim 1= RBF #
            %Dim 2= Channel for voltage
            %Dim 3= Dimension center is placed in ( what load channel it is helping approximate)
        end

        %update the approximation
        aprxINminGZ2 = aprxINminGZ2+rbfc_INminGZ;
        aprxINminGZ_Hist{u} = aprxINminGZ2;

        % SOLVE FOR TARES BY TAKING THE MEAN
        [taresAllPointsGRBF,taretalGRBFSTDDEV] = meantare(series0,aprxINminGZ2-targetMatrix0);
        taresGRBF = taresAllPointsGRBF(s_1st0,:);
        taresGRBFSTDEV = taretalGRBFSTDDEV(s_1st0,:);
        tareGRBFHist{u} = taresGRBF;

        %Calculate tare corrected load approximation
        aprxINminTARE2=aprxINminGZ2-taresAllPointsGRBF;

        %Calculate and store residuals
        targetRes2 = targetMatrix0-aprxINminTARE2;      %0=b-Ax
        newRes2 = targetRes2'*targetRes2;
        resSquare2 = diag(newRes2);
        resSquareHist(u,:) = resSquare2;
        resStdHist(u,:)=std(targetRes2);
    end

    %OUTPUT FUNCTION
    %Function creates all outputs for calibration, GRBF section
    section={'Calibration GRBF'};
    newStruct=struct('aprxINminTARE2',aprxINminTARE2,...
        'epsHist',epsHist,...
        'cHist',cHist,...
        'centerIndexHist',centerIndexHist,...
        'center_daHist',center_daHist,...
        'ANOVA',ANOVA,...
        'coeff',coeff, 'h_GRBF',h_GRBF);
    uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
        [fieldnames(uniqueOut); fieldnames(newStruct)],1);
    output(section,FLAGS,targetRes2,loadCapacities,fileName,numpts0,nseries0,...
        taresGRBF,taresGRBFSTDEV,loadlist,series0,excessVec0,dimFlag,voltagelist,...
        reslist,numBasis,pointID0,series20,output_location,REPORT_NO,uniqueOut)

end
%END CALIBRATION GRBF SECTION

%%
if FLAGS.balVal == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       VALIDATION SECTION                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Initialize structure for unique outputs for section
    uniqueOut=struct();

    load(out.savePathval,'-mat');
    [validSeries,s_1stV,~] = unique(seriesvalid);
    xvalid=coeff; %JUST USE COEFF FOR VALIDATION (NO ITERCEPTS)

    % num of data points
    numptsvalid = length(seriesvalid);
    dimFlagvalid = length(excessVecvalid(1,:));

    %find the average natural zeros (also called global zeros)
    globalZerosvalid = mean(natzerosvalid,1);

    %load capacities
    loadCapacitiesvalid(loadCapacitiesvalid == 0) = realmin;

    %find number of series0; this will tell us the number of tares
    nseriesvalid = max(seriesvalid);

    %find zero points of each series0 and number of points in a series0
    %localZerosAllPoints is the same as localZeroMatrix defined in the RBF
    %section
    globalZerosAllPointsvalid = ones(numptsvalid,1)*globalZerosvalid;

    % Subtract the Global Zeros from the Inputs and Local Zeros
    dainputsvalid = excessVecvalid-globalZerosvalid;

    %%% 5/16/18
    %Remember that  excessVec0 = excessVec0_complete - globalZerosAllPoints;
    excessVecvalidkeep = excessVecvalid  - globalZerosAllPointsvalid;
    %%%

    % Call the Algebraic Subroutine
    comINvalid = balCal_algEqns(FLAGS.model,dainputsvalid,seriesvalid,0);

    %VALIDATION APPROXIMATION
    %define the approximation for inputs minus global zeros
    aprxINvalid = comINvalid*xvalid;        %to find approximation

    %%%%% 3/23/17 Zap intercepts %%%
    aprxINminGZvalid = aprxINvalid;
    checkitvalid = aprxINminGZvalid-targetMatrixvalid;

    % SOLVE FOR TARES BY TAKING THE MEAN
    [taresAllPointsvalid,taretalstdvalid] = meantare(seriesvalid,checkitvalid);
    taresvalid     = taresAllPointsvalid(s_1stV,:);
    tares_STDEV_valid = taretalstdvalid(s_1stV,:);

    %Tare corrected approximation
    aprxINminTAREvalid=aprxINminGZvalid-taresAllPointsvalid;

    %RESIDUAL
    targetResvalid = targetMatrixvalid-aprxINminTAREvalid;

    %CALCULATE PREDICTION INTERVAL FOR POINTS
    if FLAGS.loadPI==1

        [loadPI_valid]=calc_alg_PI(ANOVA,anova_pct,comINvalid,aprxINvalid); %Calculate prediction interval for loads

        newStruct=struct('loadPI_valid',loadPI_valid);
        uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
            [fieldnames(uniqueOut); fieldnames(newStruct)],1);
    end

    %OUTPUT FUNCTION
    %Function creates all outputs for validation, algebraic section
    newStruct=struct('aprxINminTAREvalid',aprxINminTAREvalid);
    uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
        [fieldnames(uniqueOut); fieldnames(newStruct)],1);
    section={'Validation Algebraic'};
    output(section,FLAGS,targetResvalid,loadCapacitiesvalid,fileNamevalid,numptsvalid,nseriesvalid,...
        taresvalid,tares_STDEV_valid,loadlist, seriesvalid ,excessVecvalidkeep,dimFlag,voltagelist,...
        reslist,numBasis,pointIDvalid,series2valid,output_location,REPORT_NO,uniqueOut)

    %END VALIDATION ALGEBRAIC SECTION

    %%
    if FLAGS.balCal == 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                    RBF SECTION FOR VALIDATION                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %goal to use centers, width and coefficients to validate parameters against
        %independent data

        %Initialize structure for unique outputs for section
        uniqueOut=struct();

        targetRes2valid = targetResvalid;
        aprxINminGZ2valid = aprxINminGZvalid;

        %Initialize Variables
        aprxINminGZ_Histvalid = cell(numBasis,1);
        tareHistvalid = cell(numBasis,1);
        resSquareHistvalid=zeros(numBasis,dimFlagvalid);
        resStdHistvalid=zeros(numBasis,dimFlagvalid);
        for u=1:numBasis
            %Call function to place single GRBF
            [rbfc_INminGZvalid]=place_GRBF(u,dainputsvalid,epsHist,cHist,center_daHist,h_GRBF);

            %update the approximation
            aprxINminGZ2valid = aprxINminGZ2valid+rbfc_INminGZvalid;
            aprxINminGZ_Histvalid{u} = aprxINminGZ2valid;

            % SOLVE FOR TARES BY TAKING THE MEAN
            [~,s_1st,~] = unique(seriesvalid);
            [taresAllPointsvalid2,taretalstdvalid2] = meantare(seriesvalid,aprxINminGZ2valid-targetMatrixvalid);
            taresGRBFvalid = taresAllPointsvalid2(s_1st,:);
            taresGRBFSTDEVvalid = taretalstdvalid2(s_1st,:);
            tareHistvalid{u} = taresGRBFvalid;

            %Calculate tare corrected load approximation
            aprxINminTARE2valid=aprxINminGZ2valid-taresAllPointsvalid2;

            %Residuals
            targetRes2valid = targetMatrixvalid-aprxINminTARE2valid;      %0=b-Ax
            newRes2valid = targetRes2valid'*targetRes2valid;
            resSquare2valid = diag(newRes2valid);
            resSquareHistvalid(u,:) = resSquare2valid;
            resStdHistvalid(u,:)=std(targetRes2valid);
        end

        %OUTPUT FUNCTION
        %Function creates all outputs for validation, GRBF section
        section={'Validation GRBF'};
        newStruct=struct('aprxINminTARE2valid',aprxINminTARE2valid);
        uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
            [fieldnames(uniqueOut); fieldnames(newStruct)],1);
        output(section,FLAGS,targetRes2valid,loadCapacitiesvalid,fileNamevalid,numptsvalid,nseriesvalid,...
            taresGRBFvalid,taresGRBFSTDEVvalid,loadlist,seriesvalid,excessVecvalid,dimFlagvalid,voltagelist,...
            reslist,numBasis,pointIDvalid,series2valid,output_location,REPORT_NO,uniqueOut)
    end
    %END GRBF SECTION FOR VALIDATION
end
%END VALIDATION SECTION

%%
if FLAGS.balApprox == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                        APPROXIMATION SECTION                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %DEFINE THE PRODUCTION CSV INPUT FILE AND SELECT THE RANGE OF DATA VALUES TO READ
    load(out.savePathapp,'-mat');

    if FLAGS.balCal == 2 %If RBFs were placed, put parameters in structure
        GRBF.epsHist=epsHist;
        GRBF.cHist=cHist;
        GRBF.center_daHist=center_daHist;
        GRBF.h_GRBF=h_GRBF;
    else
        GRBF='GRBFS NOT PLACED';
    end

    %Function that performs all ANOVA calculations and outputs
    [aprxINminGZapprox,loadPI_approx]=AOX_approx_funct(coeff,natzerosapprox,excessVecapprox,FLAGS,seriesapprox,...
        series2approx,pointIDapprox,loadlist,output_location,GRBF,ANOVA,anova_pct);

end
%END APPROXIMATION SECTION

%File Cleanup
if FLAGS.input_save==0
    if isfield(out,'cal_create')==1
        try
            delete(out.savePathcal);
        end
    end
    if isfield(out,'val_create')==1
        try
            delete(out.savePathval);
        end
    end
    if isfield(out,'app_create')==1
        try
            delete(out.savePathapp);
        end
    end
end

%

fprintf('\n  ');
fprintf('\nCalculations Complete.\n');
fprintf('%s',strcat('Check '," ",output_location,' for output files.'))
fprintf('\n \n');

if isdeployed % Optional, use if you want the non-deployed version to exit immediately
    input('Press enter to finish and close');
end
