%requires "balCal_meritFunction.m" to run

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
tic;

%TO SELECT Algebraic Model                                  set FLAGS.balCal = 1;
%TO SELECT Algebraic and GRBF Model                         set FLAGS.balCal = 2;
FLAGS.balCal = out.grbf;
%DEFINE THE NUMBER OF BASIS FUNCTIONS
numBasis = out.basis;
%SET SELF TERMINATE OPTION FOR RBFS
pos_str={'No Self-Termination','Validation Error Self-Terminate','Prediction Interval Self-Terminate','VIF + Prediction Interval Self-Terminate'}; %Possible self-termination options
match=strcmp(pos_str,out.selfTerm_str);
FLAGS.valid_selfTerm=match(2);
FLAGS.PI_selfTerm=match(3);
FLAGS.VIF_selfTerm=match(4);
if out.valid==0
    FLAGS.valid_selfTerm=0;
end
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

FLAGS.anova = out.anova;
FLAGS.loadPI = out.loadPI;
FLAGS.BALFIT_Matrix=out.BALFIT_Matrix;
FLAGS.BALFIT_ANOVA=out.BALFIT_ANOVA;
FLAGS.Rec_Model=out.Rec_Model;
anova_pct=out.anova_pct;
FLAGS.approx_and_PI_print=out.approx_and_PI_print;
FLAGS.custom_eqn_iter=out.stableRec_FLAGcheck;

REPORT_NO=out.REPORT_NO;
file_output_location=out.output_location;

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
if FLAGS.model==6 %If user has selected a custom model
    termInclude=out.termInclude;
    %Assemble custom matrix
    customMatrix=customMatrix_builder(dimFlag,termInclude);
    customMatrix = [customMatrix; ones(nseries0,dimFlag)];
    
    %Proceed through code with custom equation
    FLAGS.model = 4;
    algebraic_model={'CUSTOM'};
elseif FLAGS.model==5
    %Build bustom equation matrix based on the balance type selected
    balanceType=out.balanceEqn;
    %Select the terms to be included
    %Terms are listed in following order:
    % F, |F|, F*F, F*|F|, F*G, |F*G|, F*|G|, |F|*G, F*F*F, |F*F*F|
    termInclude=zeros(10,1);
    if balanceType==1
        termInclude([1,3,5])=1;
        algebraic_model={'TRUNCATED (BALANCE TYPE 1-A)'};
    elseif balanceType==2
        termInclude([1,3,5,9])=1;
        algebraic_model={'BALANCE TYPE 1-B'};
    elseif balanceType==3
        termInclude([1,5])=1;
        algebraic_model={'BALANCE TYPE 1-C'};
    elseif balanceType==4
        termInclude([1,3])=1;
        algebraic_model={'BALANCE TYPE 1-D'};
    elseif balanceType==5
        termInclude([1,2,3,5])=1;
        algebraic_model={'BALANCE TYPE 2-A'};
    elseif balanceType==6
        termInclude([1,2,3,4,5])=1;
        algebraic_model={'BALANCE TYPE 2-B'};
    elseif balanceType==7
        termInclude(1:8)=1;
        algebraic_model={'BALANCE TYPE 2-C'};
    elseif balanceType==8
        termInclude(1:10)=1;
        algebraic_model={'FULL (BALANCE TYPE 2-D)'};
    elseif balanceType==9
        termInclude([1,2,4,5])=1;
        algebraic_model={'BALANCE TYPE 2-E'};
    elseif balanceType==10
        termInclude([1,2,5])=1;
        algebraic_model={'BALANCE TYPE 2-F'};
    end
    %Assemble custom matrix
    customMatrix=customMatrix_builder(dimFlag,termInclude);
    customMatrix = [customMatrix; ones(nseries0,dimFlag)];
    %Proceed through code with custom equation
    FLAGS.model = 4;
elseif FLAGS.model == 4
    % Load the custom equation matrix if using a custom algebraic model
    % SEE: CustomEquationMatrixTemplate.csv
    customMatrix = out.customMatrix;
    customMatrix = [customMatrix; ones(nseries0,dimFlag)];
    algebraic_model={'CUSTOM'};
else
    customMatrix = 1;
    if FLAGS.model == 1
        algebraic_model={'FULL (BALANCE TYPE 2-D)'};
    elseif FLAGS.model == 2
        algebraic_model={'TRUNCATED (BALANCE TYPE 1-A)'};
    elseif FLAGS.model == 3
        algebraic_model={'LINEAR'};
    end
end

% Load data labels if present, otherwise use default values.
if exist('loadlabels','var')==0 || isempty(loadlabels)==1
    loadlist = {'NF','BM','S1','S2','RM','AF','PLM', 'PCM', 'MLM', 'MCM'};
    voltagelist = {'rNF','rBM','rS1','rS2','rRM','rAF','rPLM','rPCM','rMLM','rMCM'};
    reslist = strcat('res',loadlist);
else
    loadlist = loadlabels;
    voltagelist = voltlabels;
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
fprintf('\n ********** Starting Calibration Algebraic Calculations **********\n')

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
[comIN0,high,high_CELL] = balCal_algEqns(FLAGS.model,dainputs0,series0,1,voltagelist);

%%% Balfit Stats and Regression Coeff Matrix
balfitdainputs0 = targetMatrix0;
balfittargetMatrix0 = balCal_algEqns(3,dainputs0,series0,0);
balfitcomIN0 = balCal_algEqns(FLAGS.model,balfitdainputs0,series0,1);
%%% Balfit Stats and Regression Coeff Matrix



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
        dainputs0(OUTLIER_ROWS,:)= [];
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
    reslist,numBasis,pointID0,series20,file_output_location,REPORT_NO,algebraic_model,uniqueOut)

%END CALIBRATION ALGEBRAIC SECTION

%%
if FLAGS.balVal == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       VALIDATION - ALGEBRAIC SECTION                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n ********** Starting Validation Algebraic Calculations **********\n')
    
    %Initialize structure for unique outputs for section
    uniqueOut=struct();
    
    load(out.savePathval,'-mat');
    [validSeries,s_1stV,~] = unique(seriesvalid);
    
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
    aprxINvalid = comINvalid*coeff;        %to find approximation, JUST USE COEFF FOR VALIDATION (NO ITERCEPTS)
    
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
    std_targetResvalid=std(targetResvalid);
    
    %CALCULATE PREDICTION INTERVAL FOR POINTS
    if FLAGS.loadPI==1
        
        [loadPI_valid]=calc_PI(ANOVA,anova_pct,comINvalid,aprxINvalid); %Calculate prediction interval for loads
        
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
        reslist,numBasis,pointIDvalid,series2valid,file_output_location,REPORT_NO,algebraic_model,uniqueOut)
    
end
%END VALIDATION ALGEBRAIC SECTION

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
    
    fprintf('\n ********** Starting Calibration GRBF Calculations **********\n')
    
    %Initialize structure for unique outputs for section
    uniqueOut=struct();
    add_RBFs=1;
    
    targetRes2=targetRes;
    aprxINminGZ2 = aprxINminGZ;
    
    %Initialize Variables
    aprxIN2_Hist = cell(numBasis,1);
    tareGRBFHist = cell(numBasis,1);
    cHist_tot=cell(numBasis,1);
    centerIndexLoop=zeros(1,dimFlag);
    eta=zeros(length(excessVec0(:,1)),dimFlag);
    eps=zeros(1,dimFlag);
    rbfINminGZ=zeros(length(excessVec0(:,1)),numBasis,dimFlag);
    rbfc_INminGZ=zeros(length(excessVec0(:,1)),dimFlag);
    epsHist=zeros(numBasis,dimFlag);
    centerIndexHist=zeros(numBasis,dimFlag);
    center_daHist=zeros(numBasis,dimFlag,dimFlag);
    resSquareHist=zeros(numBasis,dimFlag);
    resStdHist=zeros(numBasis,dimFlag);
    comINvalid_RBF=zeros(size(dainputsvalid,1),numBasis*dimFlag);
    dist=zeros(size(dainputs0,1),size(dainputs0,1),size(dainputs0,2));
    
    for i=1:size(dainputs0,2)
        dist(:,:,i)=dainputs0(:,i)'-dainputs0(:,i); %solve distance in each dimension, Eqn 16 from Javier's notes
    end
    R_square=sum(dist.^2,3); %Eqn 17 from Javier's notes: squared distance between each point
    R_square_find=R_square;
    R_square_find(R_square_find==0)=NaN; %Eliminate zero values (on diagonal)
    min_R_square=min(R_square_find); %Find distance to closest point
    %Set limits on width (shape factor)
    h_GRBF=sqrt(max(min(R_square_find)));
    eps_min=0.1; %Fasshauer pg 234, large epsilon= 'spiky'
    eps_max=1.0;
    
    max_mult=5;
    maxPer=ceil(max_mult*numBasis/size(dainputs0,1)); %Max number of RBFs that can be placed at any 1 location: max_mult* each point's true 'share' or RBFs
    %     maxPer=ceil(0.05*numBasis); %Max number of RBFs that can be placed at any 1 location
    maxPer=1;
    
    %Initialize self terminate variables:
    self_Terminate=false(1,dimFlag);
    RBFs_added=zeros(1,dimFlag);
    if FLAGS.valid_selfTerm==1
        resStdHistvalid=zeros(numBasis,dimFlag);
        period_change=zeros(numBasis,dimFlag);
        period_length=max([10,0.1*numBasis]); %CHANGE?
    end
    if (FLAGS.PI_selfTerm==1 || FLAGS.VIF_selfTerm==1)
        calib_PI_rms_Hist=zeros(numBasis,dimFlag);
        period_change=zeros(numBasis,dimFlag);
        period_length=max([10,0.1*numBasis]); %CHANGE?
        [loadPI_ALG]=calc_PI(ANOVA,anova_pct,comIN0(:,1:nterms),aprxIN); %Calculate prediction interval for loads
        calib_ALG_PI_rms=sqrt(sum((loadPI_ALG).^2,1)/numpts0); %RMS for calibration PI
    end
    
    if FLAGS.VIF_selfTerm==1 %Initialize variables for self terminating based on VIF
        max_VIF_hist=zeros(numBasis,dimFlag); %History variable of VIF as RBFs are added
        comIN0_RBF=comIN0; %Initialize 'X' matrix for RBF predictor variables
        VIF_lim=9.95; %Limit for acceptable VIF
        
        for i=1:dimFlag %Check Algebraic model max VIF in each channel
            if max(ANOVA(i).VIF)>VIF_lim %If Algebraic model already exceeds max VIF
                self_Terminate(i)=1; %Terminate channel initially
                fprintf(strcat('\n Channel'," ", string(i), ' Reached VIF termination criteria with Algebraic Model, no RBFs will be added in Channel'," ",string(i))); %Output message
            end
        end
        if all(self_Terminate)==1 %Check if all channels have self terminated
            fprintf(strcat('\n All Channels Reached VIF termination criteria with Algebraic Model, no RBFs will be added')); %Output message
            calib_PI_rms_Hist=calib_ALG_PI_rms;
        end
        
        %Initialize customMatrix for solving terms
        if FLAGS.model==4
            customMatrix_RBF=customMatrix;
        else
            customMatrix_RBF=[ones(nterms,dimFlag);ones(nseries0,dimFlag)];
        end
    end
    
    count=zeros(size(dainputs0)); %Initialize matrix to count how many RBFs have been placed at each location
    for u=1:numBasis
        RBFs_added(not(self_Terminate))=u; %Counter for how many RBFs have been placed in each channel
        if FLAGS.VIF_selfTerm==1 %If self terminating based on VIF
            comIN0_RBF_VIFtest=[comIN0_RBF,zeros(numpts0,1)]; %Initialize
        end
        for s=1:dimFlag %Loop places center and determines width for RBF in each channel
            if self_Terminate(s)==0 %If channel has not been self-terminated
                VIF_good=0; %Initialize Flag for if VIF is acceptable
                while VIF_good==0 %Repeat until VIF is acceptable
                    
                    %PLACE CENTER BASED ON LOCATION OF MAX RESIDUAL
                    targetRes2_find=targetRes2;
                    targetRes2_find(count(:,s)>=maxPer,s)=0; %Zero out residuals that have reached max number of RBFs
                    [~,centerIndexLoop(s)] = max(abs(targetRes2_find(:,s))); %Place RBF at max residual location
                    count(centerIndexLoop(s),s)=count(centerIndexLoop(s),s)+1; %Advance count for center location
                    
                    %DEFINE DISTANCE BETWEEN RBF CENTER AND OTHER DATAPOINTS
                    eta(:,s)=R_square(:,centerIndexLoop(s)); %Distance squared between RBF center and datapoints
                    
                    %find widths 'w' by optimization routine
                    eps(s) = fminbnd(@(eps) balCal_meritFunction2(eps,targetRes2(:,s),eta(:,s),h_GRBF,dimFlag),eps_min,eps_max );
                    
                    %DEFINE RBF W/O COEFFFICIENT FOR MATRIX ('X') OF PREDICTOR VARIABLES
                    rbfINminGZ_temp=((eps(s)^dimFlag)/(sqrt(pi^dimFlag)))*exp(-((eps(s)^2)*(eta(:,s)))/h_GRBF^2); %From 'Iterated Approximate Moving Least Squares Approximation', Fasshauer and Zhang, Equation 22
                    rbfINminGZ_temp=rbfINminGZ_temp-mean(rbfINminGZ_temp); %Bias is mean of RBF
                    
                    if FLAGS.VIF_selfTerm==1 %If self terminating based on VIF
                        comIN0_RBF_VIFtest(:,end)=rbfINminGZ_temp; %Define input matrix ('X') with new RBF, algebraic terms, and previous RBFs
                        customMatrix_RBF_VIFtest=[customMatrix_RBF(:,s);1]; %Define customMatrix for solving terms with new RBF
                        VIFtest=vif_dl(comIN0_RBF_VIFtest(:,logical(customMatrix_RBF_VIFtest))); %Calculate VIFs for predictor terms with new RBF
                        
                        max_VIF_hist(u,s)=max(VIFtest([1:end-nseries0-1,end])); %Store maximum VIF in history
                        if  max_VIF_hist(u,s)<=VIF_lim %If max VIF is <= VIF limit
                            VIF_good=1; %VIF criteria is satisfied, this will exit 'while' loop
                        else
                            if all(count(:,s)==1) %If all points have been tested and none can be added without exceeding VIF limit
                                self_Terminate(s)=1; %Self terminate channel
                                RBFs_added(s)=u-1; %RBFs added in that channel is back one from current iteration
                                rbfINminGZ_temp=0; %Zero out new RBF since VIF limit not met
                                fprintf(strcat('\n Channel'," ", string(s), ' Reached VIF termination criteria, # RBF=',string(u-1))); %Output message
                                if all(self_Terminate) %If all channels have now terminated
                                    calib_PI_rms_Hist(u:end,:)=[];
                                end
                                break %Exit while loop
                            end
                        end
                    else %If not self terminating based on VIF
                        VIF_good=1; %Exit while loop, VIF criteria not considered
                    end %END VIF_selfTerm section
                    
                end %END loop for iterating until good VIF
                rbfINminGZ(:,u,s)=rbfINminGZ_temp; %Store temp RBF
            end
        end
        
        %Make custom Matrix to solve for only RBF coefficinets in correct channel
        RBF_custom=repmat(eye(dimFlag,dimFlag),u,1);
        for i=1:dimFlag
            RBF_custom(dimFlag*RBFs_added(i)+1:end,i)=0;
        end
        if FLAGS.model==4
            customMatrix_RBF=[customMatrix(1:nterms,:);RBF_custom;customMatrix(nterms+1:end,:)];
        else
            customMatrix_RBF=[ones(nterms,dimFlag);RBF_custom;ones(nseries0,dimFlag)];
        end
        
        %Add RBFs to comIN0 variable to solve with alg coefficients
        comIN0_RBF=[comIN0(:,1:nterms),zeros(size(comIN0,1),u*dimFlag),comIN0(:,nterms+1:end)];
        for i=1:u
            comIN0_RBF(:,nterms+1+dimFlag*(i-1):nterms+dimFlag*(i))=rbfINminGZ(:,i,:);
        end
        
        %New flag structure for calc_xcalib
        FLAGS_RBF.model=4;
        if u==numBasis %If final RBF placed
            if any(self_Terminate) %If self terminated, stats will be recalculated below
                calc_channel=not(self_Terminate);
                if (FLAGS.PI_selfTerm==1 || FLAGS.VIF_selfTerm==1) %If self terminating based on Prediction Interval
                    FLAGS_RBF.anova=1; %perform ANOVA analysis
                    FLAGS_RBF.test_FLAG=1; %Do not calculate VIF for time savings
                else %If not self terminating with PI
                    FLAGS_RBF.anova=0; %do not perform ANOVA analysis
                end
            else %Otherwise, final calculation with RBFs
                FLAGS_RBF.anova=FLAGS.anova; %Calculate ANOVA based on user preference
                FLAGS_RBF.test_FLAG=0; %calculate VIF
                calc_channel=ones(1,dimFlag); %Calculate every channel
            end
            
        else %NOT final RBF Placed
            if (FLAGS.PI_selfTerm==1 || FLAGS.VIF_selfTerm==1) %If self terminating based on Prediction Interval
                FLAGS_RBF.anova=1; %perform ANOVA analysis
                FLAGS_RBF.test_FLAG=1; %Do not calculate VIF for time savings
            else
                FLAGS_RBF.anova=0; %Do not calculate ANOVA
            end
            calc_channel=not(self_Terminate); %Calculate channels that have not been terminated
        end
        nterms_RBF=nterms+u*dimFlag; %New number of terms to solve for
        
        %Calculate Algebraic and RBF coefficients with calc_xcalib function
        [xcalib_RBF, ANOVA_GRBF] = calc_xcalib(comIN0_RBF,targetMatrix0,series0,...
            nterms_RBF,nseries0,dimFlag,FLAGS_RBF,customMatrix_RBF,anova_pct,loadlist,'Direct w RBF',calc_channel);
        
        if u>1 && any(self_Terminate) %Coefficients for self terminated channels are retained from previous run
            xcalib_RBF(1:size(xcalib_RBF_last,1)-nseries0,self_Terminate)=xcalib_RBF_last(1:end-nseries0,self_Terminate);
            xcalib_RBF(end-nseries0:end,self_Terminate)=xcalib_RBF_last(end-nseries0:end,self_Terminate);
        end
        xcalib_RBF_last=xcalib_RBF;
        %Extract RBF coefficients
        coeff_algRBFmodel=xcalib_RBF(1:nterms_RBF,:); %Algebraic and RBF coefficient matrix
        coeff_algRBFmodel_alg=xcalib_RBF(1:nterms,:); %new algebraic coefficients
        coeff_algRBFmodel_RBF_diag=xcalib_RBF(nterms+1:nterms_RBF,:); %new RBF coefficients, spaced on diagonals
        %Extract only RBF coefficients in compact matrix
        coeff_algRBFmodel_RBF=zeros(u,dimFlag);
        for i=1:u
            coeff_algRBFmodel_RBF(i,:)=diag(coeff_algRBFmodel_RBF_diag(1+dimFlag*(i-1):dimFlag*i,:));
        end
        
        %Store basis parameters in Hist variables
        epsHist(u,not(self_Terminate)) = eps(not(self_Terminate));
        cHist_tot{u} = coeff_algRBFmodel;
        centerIndexHist(u,not(self_Terminate)) = centerIndexLoop(not(self_Terminate));
        for s=1:dimFlag
            if self_Terminate(s)==0
                center_daHist(u,:,s)=dainputs0(centerIndexLoop(s),:); %Variable stores the voltages of the RBF centers.
                %Dim 1= RBF #
                %Dim 2= Channel for voltage
                %Dim 3= Dimension center is placed in ( what load channel it is helping approximate)
            end
        end
        
        %update the approximation
        aprxIN2=comIN0_RBF*xcalib_RBF;
        aprxIN2_Hist{u} = aprxIN2;
        
        taresGRBF = -xcalib_RBF(nterms+u*dimFlag+1:end,:);
        taretalRBF=taresGRBF(series0,:);
        aprxINminGZ2=aprxIN2+taretalRBF; %Approximation that does not include intercept terms
        tareGRBFHist{u} = taresGRBF;
        
        %    QUESTION: JRP; IS THIS NECESSARY/USEFUL?
        [~,taresGRBF_STDDEV_all] = meantare(series0,aprxINminGZ2-targetMatrix0);
        taresGRBFSTDEV = taresGRBF_STDDEV_all(s_1st0,:);
        
        %Calculate tare corrected load approximation
        aprxINminTARE2=aprxINminGZ2-taretalRBF;
        
        %Calculate and store residuals
        targetRes2 = targetMatrix0-aprxINminTARE2;
        newRes2 = targetRes2'*targetRes2;
        resSquare2 = diag(newRes2);
        resSquareHist(u,:) = resSquare2;
        resStdHist(u,:)=std(targetRes2);
        
        %Validation Error Self-Termination Check
        if FLAGS.valid_selfTerm==1
            %Test on validation data
            comINvalid_RBF(:,(u-1)*dimFlag+1:u*dimFlag)=create_comIN_RBF(dainputsvalid,epsHist(u,:),center_daHist(u,:,:),h_GRBF); %Generate comIN for RBFs
            comINvalid_algRBF=[comINvalid, comINvalid_RBF(:,1:u*dimFlag)]; %Combine comIN from algebraic terms and RBF terms to multiply by coefficients
            
            aprxINminGZ2valid=comINvalid_algRBF*coeff_algRBFmodel; %find approximation with alg and RBF Coefficients
            
            % SOLVE FOR TARES BY TAKING THE MEAN
            [~,s_1st,~] = unique(seriesvalid);
            [taresAllPointsvalid2,taretalstdvalid2] = meantare(seriesvalid,aprxINminGZ2valid-targetMatrixvalid);
            
            %Calculate tare corrected load approximation
            aprxINminTARE2valid=aprxINminGZ2valid-taresAllPointsvalid2;
            
            %Residuals
            targetRes2valid = targetMatrixvalid-aprxINminTARE2valid;      %0=b-Ax
            newRes2valid = targetRes2valid'*targetRes2valid;
            resSquare2valid = diag(newRes2valid);
            resStdHistvalid(u,:)=std(targetRes2valid);
            
            %Self termination criteria
            %Calculate period_change, the difference between the minimum
            %error in the last n iterations and the error n+1 iterations ago
            if u>period_length
                period_change(u,:)=min(resStdHistvalid(u-(period_length-1):u,:))-resStdHistvalid(u-period_length,:);
            elseif u==period_length
                period_change(u,:)=min(resStdHistvalid(1:u,:))-std_targetResvalid;
            end
            
            %Self Terminate if validation error has only gotten worse over
            %the last n+1 iterations
            for i=1:dimFlag
                if period_change(u,i)>0 && self_Terminate(i)==0
                    fprintf(strcat('\n Channel'," ", string(i), ' Reached validation period change termination criteria, # RBF=',string(u)));
                    self_Terminate(i)=1;
                end
            end
        end
        
        %Prediction Error Self-Termination Check
        if (FLAGS.PI_selfTerm==1 || FLAGS.VIF_selfTerm==1) && any(~self_Terminate)
            [loadPI_GRBF_iter]=calc_PI(ANOVA_GRBF,anova_pct,comIN0_RBF(:,1:nterms_RBF),aprxIN2,calc_channel); %Calculate prediction interval for loads
            for i=1:dimFlag
                if self_Terminate(i)==1 %If channel is self terminated, use PI RMS from previous iteration
                    %                     calib_PI_rms_Hist(u,i)=calib_PI_rms_Hist(u-1,i);
                    calib_PI_rms_Hist(u,i)=calib_PI_rms_Hist(u-1,i);
                else
                    %                     calib_PI_rms_Hist(u,i)=sqrt(sum((ANOVA_GRBF(i).y_hat_PI).^2)/numpts0); %RMS for calibration PI: from ANOVA
                    calib_PI_rms_Hist(u,i)=sqrt(sum((loadPI_GRBF_iter(:,i)).^2)/numpts0); %RMS for calibration PI
                end
            end
            
            %Self termination criteria
            %Calculate period_change, the difference between the minimum
            %PI in the last n iterations and the PI n+1 iterations ago
            if u>period_length
                period_change(u,:)=min(calib_PI_rms_Hist(u-(period_length-1):u,:))-calib_PI_rms_Hist(u-period_length,:);
            elseif u==period_length
                period_change(u,:)=min(calib_PI_rms_Hist(1:u,:))-calib_ALG_PI_rms;
            end
            
            %Self Terminate if validation error has only gotten worse over
            %the last n+1 iterations
            for i=1:dimFlag
                if period_change(u,i)>0 && self_Terminate(i)==0
                    fprintf(strcat('\n Channel'," ", string(i), ' Reached Prediction Interval period change termination criteria, # RBF=',string(u)));
                    self_Terminate(i)=1;
                    if all(self_Terminate)
                        calib_PI_rms_Hist(u+1:end,:)=[];
                    end
                end
            end
        end
        
        if all(self_Terminate)==1 %Check if all channels have self terminated
            %Trim Variables
            aprxIN2_Hist(u+1:end)=[];
            tareGRBFHist(u+1:end)=[];
            cHist_tot(u+1:end)=[];
            rbfINminGZ(:,u+1:end,:)=[];
            epsHist(u+1:end,:)=[];
            centerIndexHist(u+1:end,:)=[];
            center_daHist(u+1:end,:,:)=[];
            resSquareHist(u+1:end,:)=[];
            resStdHist(u+1:end,:)=[];
            fprintf('\n');
            break %Exit loop placing RBFs
        end
    end
    final_RBFs_added=RBFs_added;
    
    %If validation self-termination selected, recalculate for RBF number of min
    %validation STD
    if FLAGS.valid_selfTerm==1
        fprintf(strcat('\n Trimming RBFs for minimum validation STD'));
        %Find RBF number for lowest Validation STD
        min_validSTD_num=zeros(1,dimFlag);
        for i=1:dimFlag
            %Find RBF number for lowest Validation STD
            [~,min_validSTD_num(i)]=min([std_targetResvalid(i);resStdHistvalid(1:u,i)],[],1);
            min_validSTD_num(i)=min_validSTD_num(i)-1;
            fprintf(strcat('\n Channel'," ", string(i), ' Final # RBF=',string(min_validSTD_num(i))));
        end
        final_RBFs_added=min_validSTD_num; %Final RBF model is model that results in lowest validation STD
        fprintf('\n');
    end
    
    %If PI self-termination selected, recalculate for RBF number of min
    %PI STD
    
    if (FLAGS.PI_selfTerm==1 || FLAGS.VIF_selfTerm==1)
        fprintf(strcat('\n Trimming RBFs for minimum calibration prediction interval RMS'));
        %Find RBF number for lowest calibration PI rms
        min_calibPI_num=zeros(1,dimFlag);
        final_calibPI_rms=zeros(1,dimFlag);
        for i=1:dimFlag
            %Find RBF number for lowest Validation STD
            [final_calibPI_rms(i),min_calibPI_num(i)]=min([calib_ALG_PI_rms(i);calib_PI_rms_Hist(1:end,i)],[],1);
            min_calibPI_num(i)=min_calibPI_num(i)-1;
            fprintf(strcat('\n Channel'," ", string(i), ' Final # RBF=',string(min_calibPI_num(i))));
        end
        final_RBFs_added=min_calibPI_num; %Final RBF model is model that results in lowest calib PI RMS
        fprintf('\n');
    end
    
    %If any channel self-terminated, recalculate with final ALG+GRBF model
    if any(final_RBFs_added<numBasis)
        %Make custom Matrix to solve for only RBF coefficinets in correct channel
        RBF_custom=repmat(eye(dimFlag,dimFlag),max(final_RBFs_added),1);
        for i=1:dimFlag
            RBF_custom(dimFlag*final_RBFs_added(i)+1:end,i)=0;
        end
        if FLAGS.model==4
            customMatrix_RBF=[customMatrix(1:nterms,:);RBF_custom;customMatrix(nterms+1:end,:)];
        else
            customMatrix_RBF=[ones(nterms,dimFlag);RBF_custom;ones(nseries0,dimFlag)];
        end
        
        %New flag structure for calc_xcalib
        FLAGS_RBF.model=4; %Calculate with custom model
        FLAGS_RBF.anova=FLAGS.anova; %Calculate ANOVA based on user preference
        FLAGS_RBF.test_FLAG=0; %Calculate VIF
        calc_channel=ones(1,dimFlag); %Calculate stats for every channel
        nterms_RBF=nterms+max(final_RBFs_added)*dimFlag; %New number of terms to solve for
        
        %Trim comIN
        comIN0_RBF(:,nterms_RBF+1:nterms+u*dimFlag)=[];
        %Calculate Algebraic and RBF coefficients with calc_xcalib function
        [xcalib_RBF, ANOVA_GRBF] = calc_xcalib(comIN0_RBF,targetMatrix0,series0,...
            nterms_RBF,nseries0,dimFlag,FLAGS_RBF,customMatrix_RBF,anova_pct,loadlist,'Direct w RBF',calc_channel);
        
        %Extract RBF coefficients
        coeff_algRBFmodel=xcalib_RBF(1:nterms_RBF,:); %Algebraic and RBF coefficient matrix
        coeff_algRBFmodel_alg=xcalib_RBF(1:nterms,:); %new algebraic coefficients
        coeff_algRBFmodel_RBF_diag=xcalib_RBF(nterms+1:nterms_RBF,:); %new RBF coefficients, spaced on diagonals
        %Extract only RBF coefficients in compact matrix
        coeff_algRBFmodel_RBF=zeros(max(final_RBFs_added),dimFlag);
        for i=1:max(final_RBFs_added)
            coeff_algRBFmodel_RBF(i,:)=diag(coeff_algRBFmodel_RBF_diag(1+dimFlag*(i-1):dimFlag*i,:));
        end
        
        %Update basis parameters in Hist variables
        cHist_tot{u+1} = coeff_algRBFmodel;
        for s=1:dimFlag
            epsHist(final_RBFs_added(s)+1:end,s) = 0;
            centerIndexHist(final_RBFs_added(s)+1:end,s) = 0;
            center_daHist(final_RBFs_added(s)+1:end,:,s)=0; %Variable stores the voltages of the RBF centers.
            %Dim 1= RBF #
            %Dim 2= Channel for voltage
            %Dim 3= Dimension center is placed in ( what load channel it is helping approximate)
        end
        %Trim Variables
        epsHist(max(final_RBFs_added)+1:end,:) = [];
        centerIndexHist(max(final_RBFs_added)+1:end,:) = [];
        center_daHist(max(final_RBFs_added)+1:end,:,:)=[];
        
        %update the approximation
        aprxIN2=comIN0_RBF*xcalib_RBF;
        aprxIN2_Hist{u+1} = aprxIN2;
        
        taresGRBF = -xcalib_RBF(nterms_RBF+1:end,:);
        taretalRBF=taresGRBF(series0,:);
        aprxINminGZ2=aprxIN2+taretalRBF; %Approximation that does not include intercept terms
        tareGRBFHist{u+1} = taresGRBF;
        
        %    QUESTION: JRP; IS THIS NECESSARY/USEFUL?
        [~,taresGRBF_STDDEV_all] = meantare(series0,aprxINminGZ2-targetMatrix0);
        taresGRBFSTDEV = taresGRBF_STDDEV_all(s_1st0,:);
        
        %Calculate tare corrected load approximation
        aprxINminTARE2=aprxINminGZ2-taretalRBF;
        
        %Calculate and store residuals
        targetRes2 = targetMatrix0-aprxINminTARE2;
        newRes2 = targetRes2'*targetRes2;
        resSquare2 = diag(newRes2);
        resSquareHist(u+1,:) = resSquare2;
        resStdHist(u+1,:)=std(targetRes2);
    end
    
    
    %OUTPUT FUNCTION
    %Function creates all outputs for calibration, GRBF section
    section={'Calibration GRBF'};
    newStruct=struct('aprxINminTARE2',aprxINminTARE2,...
        'epsHist',epsHist,...
        'coeff_algRBFmodel_RBF',coeff_algRBFmodel_RBF,...
        'centerIndexHist',centerIndexHist,...
        'center_daHist',center_daHist,...
        'ANOVA',ANOVA,...
        'ANOVA_GRBF', ANOVA_GRBF,...
        'coeff_algRBFmodel_alg',coeff_algRBFmodel_alg,...
        'h_GRBF',h_GRBF,...
        'numBasis',max(final_RBFs_added),...
        'nterms',nterms+max(final_RBFs_added)*dimFlag,...
        'coeff_algRBFmodel',coeff_algRBFmodel,...
        'coeff',coeff);
    uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
        [fieldnames(uniqueOut); fieldnames(newStruct)],1);
    output(section,FLAGS,targetRes2,loadCapacities,fileName,numpts0,nseries0,...
        taresGRBF,taresGRBFSTDEV,loadlist,series0,excessVec0,dimFlag,voltagelist,...
        reslist,numBasis,pointID0,series20,file_output_location,REPORT_NO,algebraic_model,uniqueOut)
    %END CALIBRATION GRBF SECTION
    
    %%
    if FLAGS.balVal == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                    RBF SECTION FOR VALIDATION                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %goal to use centers, width and coefficients to validate parameters against
        %independent data
        
        fprintf('\n ********** Starting Validation GRBF Calculations **********\n')
        %Initialize structure for unique outputs for section
        uniqueOut=struct();
        
        comINvalid_RBF=create_comIN_RBF(dainputsvalid,epsHist,center_daHist,h_GRBF); %Generate comIN for RBFs
        comINvalid_algRBF=[comINvalid, comINvalid_RBF]; %Combine comIN from algebraic terms and RBF terms to multiply by coefficients
        
        aprxINminGZ2valid=comINvalid_algRBF*coeff_algRBFmodel; %find approximation with alg and RBF Coefficients
        
        % SOLVE FOR TARES BY TAKING THE MEAN
        [~,s_1st,~] = unique(seriesvalid);
        [taresAllPointsvalid2,taretalstdvalid2] = meantare(seriesvalid,aprxINminGZ2valid-targetMatrixvalid);
        taresGRBFvalid = taresAllPointsvalid2(s_1st,:);
        taresGRBFSTDEVvalid = taretalstdvalid2(s_1st,:);
        
        %Calculate tare corrected load approximation
        aprxINminTARE2valid=aprxINminGZ2valid-taresAllPointsvalid2;
        
        %Residuals
        targetRes2valid = targetMatrixvalid-aprxINminTARE2valid;      %0=b-Ax
        newRes2valid = targetRes2valid'*targetRes2valid;
        resSquare2valid = diag(newRes2valid);
        
        %CALCULATE PREDICTION INTERVAL FOR POINTS
        if FLAGS.loadPI==1
            
            [loadPI_valid_GRBF]=calc_PI(ANOVA_GRBF,anova_pct,comINvalid_algRBF,aprxINminTARE2valid); %Calculate prediction interval for loads
            
            newStruct=struct('loadPI_valid_GRBF',loadPI_valid_GRBF);
            uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
                [fieldnames(uniqueOut); fieldnames(newStruct)],1);
        end
        
        %OUTPUT FUNCTION
        %Function creates all outputs for validation, GRBF section
        section={'Validation GRBF'};
        newStruct=struct('aprxINminTARE2valid',aprxINminTARE2valid);
        uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
            [fieldnames(uniqueOut); fieldnames(newStruct)],1);
        output(section,FLAGS,targetRes2valid,loadCapacitiesvalid,fileNamevalid,numptsvalid,nseriesvalid,...
            taresGRBFvalid,taresGRBFSTDEVvalid,loadlist,seriesvalid,excessVecvalid,dimFlagvalid,voltagelist,...
            reslist,numBasis,pointIDvalid,series2valid,file_output_location,REPORT_NO,algebraic_model,uniqueOut)
    end
    %END GRBF SECTION FOR VALIDATION
end
%END GRBF SECTION

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
        GRBF.coeff_algRBFmodel=coeff_algRBFmodel;
        GRBF.center_daHist=center_daHist;
        GRBF.h_GRBF=h_GRBF;
        GRBF.ANOVA=ANOVA_GRBF;
    else
        GRBF='GRBFS NOT PLACED';
    end
    
    %Function that performs all ANOVA calculations and outputs
    [aprxINminGZapprox,loadPI_approx]=AOX_approx_funct(coeff,natzerosapprox,excessVecapprox,FLAGS,seriesapprox,...
        series2approx,pointIDapprox,loadlist,file_output_location,GRBF,ANOVA,anova_pct);
    
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
fprintf('%s',strcat('Check '," ",file_output_location,' for output files.'))
fprintf('\n \n');
runTime=toc;
if isdeployed % Optional, use if you want the non-deployed version to exit immediately
    input('Press enter to finish and close');
end