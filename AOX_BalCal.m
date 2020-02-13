% Main Driving AOX_BalCal program
% Copyright 2019 Andrew Meade, Ali Arya Mokhtarzadeh, Javier Villarreal and John Potthoff.  All Rights Reserved.
%
% Required files to run:
%   anova.m
%   AOX_approx_funct.m
%   AOX_GUI.m
%   balCal_algEqns.m
%   balCal_meritFunction2.m
%   calc_PI.m
%   calc_xcalib.m
%   correlationPlot.m
%   create_comIN_RBF.m
%   customMatrix_builder.m
%   customMatrix_labels.m
%   ID_outliers.m
%   load_and_PI_file_output.m
%   meantare.m
%   output.m
%   plotResPages.m
%   print_approxcsv.m
%   print_dlmwrite.m
%   termSelect_GUI.m
%   AOX_GUI.fig
%   termSelect_GUI.fig
%   vif_dl.m
%   nasa.png
%   rice.png

%initialize the workspace
clc;
clearvars;
close all;

fprintf('Copyright 2019 Andrew Meade, Ali Arya Mokhtarzadeh, Javier Villarreal, and John Potthoff.  All Rights Reserved.\n')
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       USER INPUT SECTION


out = AOX_GUI; %Call GUI
if out.cancel == 1
    return
end
tic;

FLAGS.mode=out.mode; %mode==1 for Balance Calibration, mode==2 for general approximation
%TO SELECT Algebraic Model                                  set FLAGS.balCal = 1;
%TO SELECT Algebraic and GRBF Model                         set FLAGS.balCal = 2;
FLAGS.balCal = out.grbf;
%DEFINE THE NUMBER OF BASIS FUNCTIONS
numBasis = out.basis;
%GRBF EPSILON (WIDTH CONTROL)
min_eps=out.min_eps; %Fasshauer pg 234, large epsilon= 'spiky'
max_eps=out.max_eps;
%SET SELF TERMINATE OPTION FOR RBFS
pos_str={'No Early Termination','Validation Error Termination','Prediction Interval Termination','VIF + Prediction Interval Termination'}; %Possible self-termination options
match=strcmp(pos_str,out.selfTerm_str);
FLAGS.valid_selfTerm=match(2);
FLAGS.PI_selfTerm=match(3);
FLAGS.VIF_selfTerm=match(4);

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
%ANOVA OPTIONS
FLAGS.anova = out.anova;
FLAGS.loadPI = out.anova; %Previous separate option in GUI, now PI is calculated for valid/approx automatically if ANOVA is performed
if FLAGS.mode==1
    FLAGS.BALFIT_Matrix=out.BALFIT_Matrix;
    FLAGS.BALFIT_ANOVA=out.BALFIT_ANOVA;
else
    FLAGS.BALFIT_Matrix=0;
    FLAGS.BALFIT_ANOVA=0;
end
FLAGS.Rec_Model=out.Rec_Model;
anova_pct=out.anova_pct;
FLAGS.approx_and_PI_print=out.approx_and_PI_print;

%ALG Model Refinement Options
zero_threshold=out.zero_threshold; %In Balfit: MATH MODEL SELECTION THRESHOLD IN % OF CAPACITY. This variable is used in performing SVD. Datapoints where the gage output (voltage) is less
%... then the threshold as a percentage of gage capacity are set to zero
%for constructing comIN and performing SVD
FLAGS.high_con=out.high_con; %Flag for enforcing term hierarchy constraint
VIFthresh=out.VIF_thresh; %Threshold for max allowed VIF
FLAGS.search_metric=out.search_metric; %Search metric for recommended math model optimization
sig_pct=out.sig_pct; %Percent confidence for designating terms as significant

FLAGS.svd=0; %Flag for performing SVD for permitted math model
FLAGS.sugEqnLeg=0; %Flag from performing search for legacy suggested equation
FLAGS.sugEqnNew=0; %Flag from performing search for updated suggested equation
FLAGS.back_recEqn=0; %Flag from performing search for recommended equation
FLAGS.forward_recEqn=0; %Flag from performing search for recommended equation
if out.AlgModel_opt>1
    FLAGS.svd=1;
end
if out.AlgModel_opt==3
    FLAGS.sugEqnLeg=1;
elseif out.AlgModel_opt==4
    FLAGS.sugEqnNew=1;
elseif out.AlgModel_opt==5
    FLAGS.forward_recEqn=1;
elseif out.AlgModel_opt==6
    FLAGS.back_recEqn=1;
end

if out.AlgModel_opt<3
    FLAGS.high_con=0; %Not in mode where hierarchy is enforced
end

%Intercept Options
if FLAGS.mode==1 %Intercept options for Balance Calibration Mode
    if out.intercept==1 %Include series intercepts
        FLAGS.glob_intercept=0;
        FLAGS.tare_intercept=1;
    elseif out.intercept==2 %Include global intercept
        FLAGS.glob_intercept=1;
        FLAGS.tare_intercept=0;
    elseif out.intercept==3 %Include no intercepts
        FLAGS.glob_intercept=0;
        FLAGS.tare_intercept=0;
    end
else %If in general approximation mode
    FLAGS.tare_intercept=0;
    if out.intercept==1 %Include global intercept
        FLAGS.glob_intercept=1;
    else %Include no intercepts
        FLAGS.glob_intercept=0;
    end
end

REPORT_NO=out.REPORT_NO;
file_output_location=out.output_location;

%TO SAVE .MAT FILE OF CALIBRATION MODEL
FLAGS.calib_model_save=out.calib_model_save_FLAG;
%TO SAVE INPUT .CAL, .VAL, .APP FILE IN OUTPUT LOCATION
FLAGS.input_save=out.input_save_FLAG;
%                       END USER INPUT SECTION
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       INITIALIZATION SECTION
fprintf('\n Starting Test: '); fprintf(out.REPORT_NO); fprintf('\n');
fprintf('\nWorking ...\n')

% Load data and characterize series
load(out.savePathcal,'-mat');
if FLAGS.mode~=1 %If not in balance calibration mode
    series=ones(size(excessVec0,1),1);
    series2=ones(size(excessVec0,1),1);
end
if exist( 'pointID', 'var')==0
    pointID=cellstr([repmat('P-',size(excessVec0,1),1),num2str((1:size(excessVec0,1))')]);
end

%Check if gage capacities are provided
if exist('gageCapacities','var')==0 || any(gageCapacities==0)
    gageCapacities=max(abs(excessVec0),[],1);
    if FLAGS.mode==1
        warning('Unable to read gage capacities for calibration data. Using maximum absolute value of voltage as gage capacity.');
    end
end

%Check if load capacities are provided
if FLAGS.mode==1 && (exist('loadCapacities','var')==0 || any(loadCapacities==0))
    loadCapacities=max(abs(targetMatrix0),[],1);
    warning('Unable to read load capacities for calibration data. Using maximum absolute value of applied loads as load capacity.');
end

series0 = series;
series20=series2;
pointID0=pointID;
[~,s_1st0,~] = unique(series0);
nseries0 = length(s_1st0);
[numpts0, voltdimFlag] = size(excessVec0); %Size of voltage input (input variables)
loaddimFlag=size(targetMatrix0,2); %Dimension of load input (desired output variable)

if loaddimFlag==voltdimFlag && FLAGS.BALFIT_Matrix==1 %If number of voltage channels is equal to load channels
    FLAGS.calc_balfit=1; %Possible to calculate balfit statistics
else
    FLAGS.calc_balfit=0; %If mismatch, not possible
end

% Loads:
% loadlabes, voltlabels (if they exist)
% loadCapacities, natzeros, targetMatrix0, excessVec0, series0
if FLAGS.model==6 %If user has selected a custom model
    termInclude=out.termInclude;
    %Assemble custom matrix
    customMatrix=customMatrix_builder(voltdimFlag,termInclude,loaddimFlag,FLAGS.glob_intercept);
    
    %Proceed through code with custom equation
    FLAGS.model = 4;
    algebraic_model={'CUSTOM TERM SELECTION'};
elseif FLAGS.model==5
    %Build bustom equation matrix based on the balance type selected
    balanceType=out.balanceEqn;
    %Select the terms to be included
    %Terms are listed in following order:
    % (INTERCEPT), F, |F|, F*F, F*|F|, F*G, |F*G|, F*|G|, |F|*G, F*F*F, |F*F*F|, F*G*G, F*G*H
    termInclude=zeros(12,1); %Tracker for terms to be included, not including intercept
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
        algebraic_model={'BALANCE TYPE 2-D'};
    elseif balanceType==9
        termInclude([1,2,4,5])=1;
        algebraic_model={'BALANCE TYPE 2-E'};
    elseif balanceType==10
        termInclude([1,2,5])=1;
        algebraic_model={'BALANCE TYPE 2-F'};
    end
    %Assemble custom matrix
    customMatrix=customMatrix_builder(voltdimFlag,termInclude,loaddimFlag,FLAGS.glob_intercept);
    %Proceed through code with custom equation
    FLAGS.model = 4;
elseif FLAGS.model == 4
    % Load the custom equation matrix if using a custom algebraic model
    % SEE: CustomEquationMatrixTemplate.csv
    customMatrix = out.customMatrix;
    algebraic_model={'CUSTOM INPUT FILE'};
else
    %Standard Full, truncated, linear model, or no algebraic model
    %Select the terms to be included
    %Terms are listed in following order:
    % (INTERCEPT), F, |F|, F*F, F*|F|, F*G, |F*G|, F*|G|, |F|*G, F*F*F, |F*F*F|, F*G*G, F*G*H
    termInclude=zeros(12,1);
    if FLAGS.model==3 %Linear eqn
        termInclude(1)=1; %Include only linear terms
        algebraic_model={'LINEAR'};
    elseif FLAGS.model==2 %Truncated Eqn type
        termInclude([1,3,5])=1;
        algebraic_model={'TRUNCATED (BALANCE TYPE 1-A)'};
    elseif FLAGS.model==1 %Full Eqn type
        termInclude(1:12)=1;
        algebraic_model={'FULL'};
    elseif FLAGS.model==0 %No Algebraic Model
        algebraic_model={'NO ALGEBRAIC MODEL'};
    end
    %Assemble custom matrix
    customMatrix=customMatrix_builder(voltdimFlag,termInclude,loaddimFlag,FLAGS.glob_intercept);
    %Proceed through code with custom equation
    FLAGS.model = 4;
end

if FLAGS.tare_intercept==1
    customMatrix = [customMatrix; ones(nseries0,loaddimFlag)];
end

% Load (output) data labels if present, otherwise use default values.
if exist('loadlabels','var')==0 || isempty(loadlabels)==1
    if FLAGS.mode==1 %If in balance calibration mode
        loadlist = {'NF','BM','S1','S2','RM','AF','PLM', 'PCM', 'MLM', 'MCM'};
        voltagelist = {'rNF','rBM','rS1','rS2','rRM','rAF','rPLM','rPCM','rMLM','rMCM'};
    else %If in general approximation mode
        loadlist=cell(1,loaddimFlag);
        voltagelist=cell(1,voltdimFlag);
        for i=1:loaddimFlag
            loadlist{i} = strcat('OUT',num2str(i));
        end
        for i=1:voltdimFlag
            voltagelist{i} = strcat('INPUT',num2str(i));
        end
    end
else
    loadlist = loadlabels;
    voltagelist = voltlabels;
end
reslist = strcat('res',loadlist);
% Voltage data labels if present, otherwise use default values.
if exist('loadunits','var')==0 
    if FLAGS.mode==1 && loaddimFlag<=10 && voltdimFlag<=10 %If in balance calibration mode
        loadunits = {'lbs','in-lbs','lbs','lbs','in-lbs','lbs','in-lbs', 'in-lbs', 'in-lbs', 'in-lbs'};
        voltunits = {'microV/V','microV/V','microV/V','microV/V','microV/V','microV/V','microV/V','microV/V','microV/V','microV/V'};
    else %If in general approaximation mode
        loadunits = repmat({'-'},1,loaddimFlag);
        voltunits = repmat({'-'},1,voltdimFlag);
    end
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
if FLAGS.mode==1 %If in balance calibration mode
    globalZeros = mean(natzeros,1);
    % Subtracts global zeros from signal.
    dainputs0 = excessVec0 - globalZeros;
else %If in general function approximation mode
    dainputs0 = excessVec0; %No 'natural zeros'
end

% The Custom model calculates all terms, and then excludes them in
% the calibration process as determined by the customMatrix.
nterms = 2*voltdimFlag*(voltdimFlag+2)+factorial(voltdimFlag)/factorial(voltdimFlag-2)+factorial(voltdimFlag)/(factorial(3)*factorial(voltdimFlag-3))+1;

% Creates the algebraic combination terms of the inputs.
% Also creates intercept terms; a different intercept for each series.
[comIN0,high,high_CELL] = balCal_algEqns(FLAGS.model,dainputs0,series0,FLAGS.tare_intercept,voltagelist);
if FLAGS.tare_intercept==1
    dainputs_zero=zeros(1,voltdimFlag); %Artificial datapoint at all zero voltage
    comIN_zero = balCal_algEqns(FLAGS.model,dainputs_zero,1,0,voltagelist); %Algebraic combination terms of zero point
end

%Define 'required' custom Matrix: minimum terms that must be included
customMatrix_req=zeros(size(comIN0,2),loaddimFlag);
if FLAGS.mode==1 %If in balance calibration mode, model must include linear voltage from corresponding channel
    lin_ind=sub2ind(size(customMatrix_req),1+(1:loaddimFlag),1:loaddimFlag); %Indices for linear voltage terms in each channel
    customMatrix_req(lin_ind)=customMatrix(lin_ind); %Must include linear voltage from channel if included in provided terms
end
if FLAGS.glob_intercept==1 %If global intercept term is used
    %Custom Matrix must include linear voltage from each respective
    %channel, and global intercept term
    customMatrix_req(1,:)=1; %Must include global intercept term
elseif FLAGS.tare_intercept==1
    %Custom Matrix must include linear voltage from each respective
    %channel, and all series intercepts (for tares)
    customMatrix_req(nterms+1:end,:)=1; %Must include series intercepts
end

%%% Balfit Stats and Regression Coeff Matrix
if FLAGS.calc_balfit==1
    balfitdainputs0 = targetMatrix0;
    balfittargetMatrix0 = balCal_algEqns(3,dainputs0,series0,0);
    balfittargetMatrix0=balfittargetMatrix0(:,2:end);
    balfitcomIN0 = balCal_algEqns(FLAGS.model,balfitdainputs0,series0,FLAGS.tare_intercept);
    %%% Balfit Stats and Regression Coeff Matrix
end

%Creates vectors that will not have outliers removed for balfit
series = series0;
targetMatrix = targetMatrix0;
comIN = comIN0;

%% Use SVD to test for permitted math model
if FLAGS.svd==1
    [customMatrix_permitted, FLAGS]=SVD_permittedEqn(customMatrix, customMatrix_req, voltdimFlag, loaddimFlag, dainputs0, FLAGS, targetMatrix0, series0, voltagelist, zero_threshold, gageCapacities); %Call function to determine permitted eqn
    customMatrix_orig=customMatrix; %Store original customMatrix
    customMatrix=customMatrix_permitted; %Proceed with permitted custom eqn
end

%Check for if permitted math model is hierarchically supported.  If not,
%unable to enforce hierarchy constraint in further model optimization
if FLAGS.high_con>0 %If enforcing hierarchy constraint
    for i=1:loaddimFlag
        incTerms=customMatrix(1:nterms,i); %Terms included in permitted model
        supTerms=any(high(boolean(customMatrix(1:nterms,i)),:),1); %Terms needed for variable support
        
        if any(~incTerms(boolean(supTerms))) %If any terms needed for support are not included in the permitted model
            FLAGS.high_con=0; %Cannot enforce hierarchy constraint
            warning('Permitted Math Model is not hierarchically supported. Unable to enforce hierarchy constraint in ALG Model Refinement.');
            break
        end
    end
    
end
%% Find suggested Eqn using Reference Balfit B29 method
if FLAGS.sugEqnLeg==1
    [customMatrix_sug, FLAGS]=modelOpt_suggested(VIFthresh, customMatrix, customMatrix_req, loaddimFlag, nterms, comIN0, sig_pct, targetMatrix0, high, FLAGS);
    customMatrix=customMatrix_sug;
end

%% Find suggested Eqn using new method
if FLAGS.sugEqnNew==1
    [customMatrix_sug, FLAGS]=modelOpt_suggestedNew(VIFthresh, customMatrix, customMatrix_req, loaddimFlag, nterms, comIN0, sig_pct, targetMatrix0, high, FLAGS);
    customMatrix=customMatrix_sug;
end
%% Find recommended Eqn using 'backward elimination' method
if FLAGS.back_recEqn==1
    [customMatrix_rec,FLAGS]=modelOpt_backward(VIFthresh, customMatrix, customMatrix_req, loaddimFlag, nterms, comIN0, sig_pct, targetMatrix0, high, FLAGS);
    customMatrix=customMatrix_rec;
end

%% Find recommended Eqn using 'forward selection' method
%User preferences
FLAGS.VIF_stop=1; %Terminate search once VIF threshold is exceeded
if FLAGS.forward_recEqn==1
    [customMatrix_rec,FLAGS]=modelOpt_forward(VIFthresh, customMatrix, customMatrix_req, loaddimFlag, nterms, comIN0, sig_pct, targetMatrix0, high, FLAGS);
    customMatrix=customMatrix_rec;
end

%% Resume calibration
%Calculate xcalib (coefficients)
[xcalib, ANOVA] = calc_xcalib(comIN       ,targetMatrix       ,series,...
    nterms,nseries0,loaddimFlag,FLAGS,customMatrix,anova_pct,loadlist,'Direct');

if FLAGS.calc_balfit==1
    %Calculate balfit coefficeints
    [balfitxcalib, balfitANOVA] = calc_xcalib(balfitcomIN0,balfittargetMatrix0,series,...
        nterms,nseries0,voltdimFlag,FLAGS,customMatrix,anova_pct,voltagelist,'BALFIT');
end

% APPROXIMATION
% define the approximation for inputs minus global zeros (includes
% intercept terms)
aprxIN = comIN0*xcalib;

% RESIDUAL
targetRes = targetMatrix0-aprxIN;

% Identify Outliers After Filtering
% (Threshold approach)
if FLAGS.balOut == 1
    fprintf('\n Identifying Outliers....')
    
    %Identify outliers based on residuals
    [OUTLIER_ROWS,num_outliers,prcnt_outliers,rowOut,colOut] = ID_outliers(targetRes,numpts0,numSTD,FLAGS);
    
    %Store outlier specific variables for output
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
            nterms,nseries0,loaddimFlag,FLAGS,customMatrix,anova_pct,loadlist,'Direct');
        
        if FLAGS.calc_balfit==1
            %%% Balfit Stats and Regression Coeff Matrix
            [balfitxcalib, balfitANOVA] = calc_xcalib(balfitcomIN0,balfittargetMatrix0,series,...
                nterms,nseries0,voltdimFlag,FLAGS,customMatrix,anova_pct,voltagelist,'BALFIT');
            %%% Balfit Stats and Matrix
        end
        
        % APPROXIMATION
        % define the approximation for inputs minus global zeros (includes
        % intercept terms)
        aprxIN = comIN0*xcalib;
        
        % RESIDUAL
        targetRes = targetMatrix0-aprxIN;
    end
end

% Splits xcalib into Coefficients and Intercepts (which are negative Tares)
coeff = xcalib(1:nterms,:);
if FLAGS.tare_intercept==1 %If tare loads were included in regression
    tares = -xcalib(nterms+1:end,:);
else
    tares=zeros(nseries0,loaddimFlag); %Else set to zero (no series intercepts)
end
intercepts=-tares;
taretal=tares(series0,:);
aprxINminGZ=aprxIN+taretal; %Approximation that does not include intercept terms

%    QUESTION: JRP; IS THIS NECESSARY/USEFUL?
if FLAGS.tare_intercept==1
    [~,tares_STDDEV_all] = meantare(series0,aprxINminGZ-targetMatrix0);
else
    tares_STDDEV_all=zeros(size(targetMatrix0));
end
tares_STDDEV = tares_STDDEV_all(s_1st0,:);

if FLAGS.calc_balfit==1
    %%% Balfit Stats and Regression Coeff Matrix
    balfit_C1INV = xcalib((1:voltdimFlag), :);
    balfit_D1 = zeros(voltdimFlag,voltdimFlag);
    balfit_INTERCEPT = globalZeros;
    balfit_C1INVC2 = balfitxcalib((voltdimFlag+1:nterms), :)*balfit_C1INV;
    balfit_regress_matrix = [globalZeros ; balfit_INTERCEPT ; balfit_C1INV ; balfit_D1 ; balfit_C1INVC2 ];
    %%% Balfit Stats and Matrix
end

if out.model~=0 %If any algebraic terms included
    %OUTPUT FUNCTION
    %Function creates all outputs for calibration, algebraic section
    section={'Calibration Algebraic'};
        newStruct=struct('aprxIN',aprxIN,...
            'coeff',coeff,...
            'nterms',nterms,...
            'ANOVA',ANOVA,...
            'targetMatrix0',targetMatrix0,...
            'loadunits',{loadunits(:)},...
            'voltunits',{voltunits(:)},...
            'description',description,...
            'gageCapacities',gageCapacities);
        uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
            [fieldnames(uniqueOut); fieldnames(newStruct)],1);
        
        if FLAGS.mode==1 %Outputs for balance calibration mode
            newStruct=struct('loadCapacities',loadCapacities,...
                'tares',tares, 'balance_type',balance_type,...
                'tares_STDDEV',tares_STDDEV);
            uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
                [fieldnames(uniqueOut); fieldnames(newStruct)],1);
        end
        
        if FLAGS.calc_balfit==1
            newStruct=struct('balfitcomIN',balfitcomIN0,...
                'balfitxcalib',balfitxcalib,...
                'balfittargetMatrix',balfittargetMatrix0,...
                'balfitANOVA',balfitANOVA,...
                'balfit_regress_matrix',balfit_regress_matrix);
            uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
                [fieldnames(uniqueOut); fieldnames(newStruct)],1);
        end
        
        %Output results from calibration algebraic section
        output(section,FLAGS,targetRes,fileName,numpts0,nseries0,...
            loadlist,series0,excessVec0,voltdimFlag,loaddimFlag,voltagelist,...
            reslist,numBasis,pointID0,series20,file_output_location,REPORT_NO,algebraic_model,uniqueOut)
else
    fprintf('   NO ALGEBRAIC MODEL INCLUDED \n');
end
%END CALIBRATION ALGEBRAIC SECTION

%%
if FLAGS.balVal == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       VALIDATION - ALGEBRAIC SECTION                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n ********** Starting Validation Algebraic Calculations **********\n')
    
    %Initialize structure for unique outputs for section
    uniqueOut=struct();
    
    load(out.savePathval,'-mat'); %Load validation data
    if FLAGS.mode~=1 %If not in balance calibration mode
        %All data is from same series
        seriesvalid=ones(size(excessVecvalid,1),1); 
        series2valid=ones(size(excessVecvalid,1),1);
    end
    [validSeries,s_1stV,~] = unique(seriesvalid); %Define series for validation data
    
    % Dimensions of data
    [numptsvalid,voltdimFlagvalid] = size(excessVecvalid); %Number of datapoints and voltage channels
    loaddimFlagvalid=size(targetMatrixvalid,2); %Dimension of load input (desired output variable)
    
    if exist( 'pointIDvalid', 'var')==0 %Create standard point ID labels if they don't exist
        pointIDvalid=cellstr([repmat('P-',numptsvalid,1),num2str((1:numptsvalid)')]);
    end
    
    if voltdimFlag~=voltdimFlagvalid || loaddimFlag~= loaddimFlagvalid %Check if mismatch between calibration and validation data.  If so, exit program
        fprintf('\n  ');
        fprintf('\n MISMATCH IN CALIBRATION/VALIDATION DATA DIMENSIONS.  UNABLE TO PROCEED.\n');
        fprintf('\n');
        if isdeployed % Optional, use if you want the non-deployed version to not exit immediately
            input('Press enter to finish and close');
        end
        return; %Quit run
    end
    
    if FLAGS.mode==1
        %find the average natural zeros (also called global zeros)
        globalZerosvalid = mean(natzerosvalid,1);
        % Subtract the Global Zeros from the Inputs and Local Zeros
        dainputsvalid = excessVecvalid-globalZerosvalid;
    
        %load capacities
        loadCapacitiesvalid(loadCapacitiesvalid == 0) = realmin;
    else
        dainputsvalid = excessVecvalid;
    end
    
    %find number of series0; this will tell us the number of tares
    nseriesvalid = max(seriesvalid);
          
    % Call the Algebraic Subroutine
    comINvalid = balCal_algEqns(FLAGS.model,dainputsvalid,seriesvalid,0); %Generate term combinations
    
    %VALIDATION APPROXIMATION
    %define the approximation for inputs minus global zeros
    aprxINvalid = comINvalid*coeff;        %to find approximation, JUST USE COEFF FOR VALIDATION (NO ITERCEPTS)
    
    %%%%% 3/23/17 Zap intercepts %%%
    aprxINminGZvalid = aprxINvalid;
    checkitvalid = aprxINminGZvalid-targetMatrixvalid;
    
    if FLAGS.tare_intercept==1 %If including Tare loads
        % SOLVE FOR TARES BY TAKING THE MEAN
        [taresAllPointsvalid,taretalstdvalid] = meantare(seriesvalid,checkitvalid);
    else
        taresAllPointsvalid=zeros(size(checkitvalid));
        taretalstdvalid=zeros(size(checkitvalid));
    end
    
    taresvalid     = taresAllPointsvalid(s_1stV,:);
    tares_STDEV_valid = taretalstdvalid(s_1stV,:);
    
    %Tare corrected approximation
    aprxINminTAREvalid=aprxINminGZvalid-taresAllPointsvalid;
    
    %RESIDUAL
    targetResvalid = targetMatrixvalid-aprxINminTAREvalid;
    std_targetResvalid=std(targetResvalid);
    
    if out.model~=0 %If any algebraic terms included
        %CALCULATE PREDICTION INTERVAL FOR POINTS
        if FLAGS.loadPI==1
            [loadPI_valid]=calc_PI(ANOVA,anova_pct,comINvalid,aprxINvalid); %Calculate prediction interval for loads
            
            %Save variables for output
            newStruct=struct('loadPI_valid',loadPI_valid);
            uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
                [fieldnames(uniqueOut); fieldnames(newStruct)],1);
        end
        
        %OUTPUT FUNCTION
        %Function creates all outputs for validation, algebraic section
        newStruct=struct('aprxINminTAREvalid',aprxINminTAREvalid);
        uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
            [fieldnames(uniqueOut); fieldnames(newStruct)],1);
        
        if FLAGS.mode==1
            newStruct=struct('loadCapacitiesvalid',loadCapacitiesvalid,...
                'tares',taresvalid,'tares_STDDEV',tares_STDEV_valid);
            uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
                [fieldnames(uniqueOut); fieldnames(newStruct)],1);
        end
        
        section={'Validation Algebraic'};
        output(section,FLAGS,targetResvalid,fileNamevalid,numptsvalid,nseriesvalid,...
            loadlist, seriesvalid ,excessVecvalid,voltdimFlag,loaddimFlag,voltagelist,...
            reslist,numBasis,pointIDvalid,series2valid,file_output_location,REPORT_NO,algebraic_model,uniqueOut)
    else
        fprintf('   NO ALGEBRAIC MODEL INCLUDED \n');
    end
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
    
    %Initialize Variables
    targetRes2=targetRes;
    aprxINminGZ2 = aprxINminGZ;
    aprxIN2_Hist = cell(numBasis,1);
    tareGRBFHist = cell(numBasis,1);
    cHist_tot=cell(numBasis,1);
    centerIndexLoop=zeros(1,loaddimFlag);
    eta=zeros(length(excessVec0(:,1)),loaddimFlag);
    eps=zeros(1,loaddimFlag);
    rbfINminGZ=zeros(length(excessVec0(:,1)),numBasis,loaddimFlag);
    rbfc_INminGZ=zeros(length(excessVec0(:,1)),loaddimFlag);
    epsHist=zeros(numBasis,loaddimFlag);
    centerIndexHist=zeros(numBasis,loaddimFlag);
    center_daHist=zeros(numBasis,voltdimFlag,loaddimFlag);
    resSquareHist=zeros(numBasis,loaddimFlag);
    resStdHist=zeros(numBasis,loaddimFlag);
    dist=zeros(size(dainputs0,1),size(dainputs0,1),size(dainputs0,2));
    
    
    for i=1:size(dainputs0,2)
        dist(:,:,i)=dainputs0(:,i)'-dainputs0(:,i); %solve distance between each datapoint in each dimension, Eqn 16 from Javier's notes
    end
    R_square=sum(dist.^2,3); %Eqn 17 from Javier's notes: squared distance between each point
    R_square_find=R_square; %Save copy of distance matrix
    R_square_find(R_square_find==0)=NaN; %Eliminate zero values (on diagonal)
    min_R_square=min(R_square_find); %Find distance to closest point
    %Set limits on width (shape factor)
    h_GRBF=sqrt(max(min(R_square_find))); %Point spacing parameter
    %     eps_min=0.1; %Fasshauer pg 234, large epsilon= 'spiky'
    %     eps_max=1.0;
    %     min_eps=0.07; %Fasshauer pg 234, large epsilon= 'spiky'
    %     max_eps=.1;
    
    max_mult=5; %CHANGE
    maxPer=ceil(max_mult*numBasis/size(dainputs0,1)); %Max number of RBFs that can be placed at any 1 location: max_mult* each point's true 'share' or RBFs
    %     maxPer=ceil(0.05*numBasis); %Max number of RBFs that can be placed at any 1 location
    maxPer=1;
    
    %Initialize self terminate variables:
    self_Terminate=false(1,loaddimFlag); %Logical vector that stores if RBF addition has been terminated in each channel
    RBFs_added=zeros(1,loaddimFlag); %Count of how many RBFs have been added in each channel
    if FLAGS.valid_selfTerm==1 %If RBF addition will be terminated based on validation data error
        resStdHistvalid=zeros(numBasis,loaddimFlag); %History of valid data residual standard deviation vs RBF number
        period_change=zeros(numBasis,loaddimFlag); %Storge for Change in validation standard deviation over last 'n' additions
        period_length=max([10,0.1*numBasis]); %CHANGE?: Length of period for self termination
        comINvalid_RBF=zeros(size(dainputsvalid,1),numBasis*loaddimFlag);
    end
    if (FLAGS.PI_selfTerm==1 || FLAGS.VIF_selfTerm==1) %If RBF addition will be terminated based on VIF or Prediction interval
        calib_PI_rms_Hist=zeros(numBasis,loaddimFlag); %History of prediction interval RMS vs RBF number
        period_change=zeros(numBasis,loaddimFlag);  %Storge for Change in PI over last 'n' additions
        period_length=max([10,0.1*numBasis]); %CHANGE?: Length of period for self termination
        if out.model~=0
            [loadPI_ALG]=calc_PI(ANOVA,anova_pct,comIN0(:,1:nterms),aprxIN); %Calculate prediction interval for loads with algebraic model
        else
            loadPI_ALG=Inf(size(targetMatrix0));
        end
        calib_ALG_PI_rms=sqrt(sum((loadPI_ALG).^2,1)/numpts0); %RMS for calibration PI
    end
    
    if FLAGS.VIF_selfTerm==1 %Initialize variables for self terminating based on VIF
        max_VIF_hist=zeros(numBasis,loaddimFlag); %History variable of VIF as RBFs are added
        comIN0_RBF=comIN0; %Initialize 'X' matrix for RBF predictor variables
        VIF_lim=9.95; %Limit for acceptable VIF
        
        for i=1:loaddimFlag %Check Algebraic model max VIF in each channel
            if out.model~=0
                max_VIF_alg=max(ANOVA(i).VIF);
            else
                max_VIF_alg=1;
            end
            if max_VIF_alg>VIF_lim %If Algebraic model already exceeds max VIF
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
            customMatrix_RBF=[ones(nterms,loaddimFlag);ones(nseries0,loaddimFlag)];
        end
    end
    %END SELF-TERMINATION INITIALIZATION
    
    count=zeros(size(dainputs0)); %Initialize matrix to count how many RBFs have been placed at each location
    for u=1:numBasis
        RBFs_added(not(self_Terminate))=u; %Counter for how many RBFs have been placed in each channel
        if FLAGS.VIF_selfTerm==1 %If self terminating based on VIF
            comIN0_RBF_VIFtest=[comIN0_RBF,zeros(numpts0,1)]; %Initialize
        end
        for s=1:loaddimFlag %Loop places center and determines width for RBF in each channel
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
                    eps(s) = fminbnd(@(eps) balCal_meritFunction2(eps,targetRes2(:,s),eta(:,s),h_GRBF,voltdimFlag),min_eps,max_eps );

                    %DEFINE RBF W/O COEFFFICIENT FOR MATRIX ('X') OF PREDICTOR VARIABLES
                    rbfINminGZ_temp=exp(-((eps(s)^2)*(eta(:,s)))/h_GRBF^2); %From 'Iterated Approximate Moving Least Squares Approximation', Fasshauer and Zhang, Equation 22
%                     rbfINminGZ_temp=((eps(s)^voltdimFlag)/(sqrt(pi^voltdimFlag)))*exp(-((eps(s)^2)*(eta(:,s)))/h_GRBF^2); %From 'Iterated Approximate Moving Least Squares Approximation', Fasshauer and Zhang, Equation 22
                    %                     rbfINminGZ_temp=rbfINminGZ_temp-mean(rbfINminGZ_temp); %Bias is mean of RBF
                    
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
                                    calib_PI_rms_Hist(u:end,:)=[]; %Trim history variable
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
        RBF_custom=repmat(eye(loaddimFlag,loaddimFlag),u,1);
        for i=1:loaddimFlag
            RBF_custom(loaddimFlag*RBFs_added(i)+1:end,i)=0;
        end
        if FLAGS.model==4
            customMatrix_RBF=[customMatrix(1:nterms,:);RBF_custom;customMatrix(nterms+1:end,:)];
        else
            customMatrix_RBF=[ones(nterms,loaddimFlag);RBF_custom;ones(nseries0,loaddimFlag)];
        end
        
        %Add RBFs to comIN0 variable to solve with alg coefficients
        comIN0_RBF=[comIN0(:,1:nterms),zeros(size(comIN0,1),u*loaddimFlag),comIN0(:,nterms+1:end)];
        for i=1:u
            comIN0_RBF(:,nterms+1+loaddimFlag*(i-1):nterms+loaddimFlag*(i))=rbfINminGZ(:,i,:);
        end
        
        %New flag structure for calc_xcalib
        FLAGS_RBF=FLAGS; %Initialize as global flag structure
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
                calc_channel=ones(1,loaddimFlag); %Calculate every channel
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
        nterms_RBF=nterms+u*loaddimFlag; %New number of terms to solve for
        
        %Calculate Algebraic and RBF coefficients with calc_xcalib function
        [xcalib_RBF, ANOVA_GRBF] = calc_xcalib(comIN0_RBF,targetMatrix0,series0,...
            nterms_RBF,nseries0,loaddimFlag,FLAGS_RBF,customMatrix_RBF,anova_pct,loadlist,'Direct w RBF',calc_channel);
        
        if u>1 && any(self_Terminate) %Coefficients for self terminated channels are retained from previous run for channels that have self terminated
            xcalib_RBF(1:size(xcalib_RBF_last,1)-nseries0,self_Terminate)=xcalib_RBF_last(1:end-nseries0,self_Terminate);
            xcalib_RBF(end-nseries0:end,self_Terminate)=xcalib_RBF_last(end-nseries0:end,self_Terminate);
        end
        xcalib_RBF_last=xcalib_RBF;
        
        %Store basis parameters in Hist variables
        epsHist(u,not(self_Terminate)) = eps(not(self_Terminate));
%         cHist_tot{u} = coeff_algRBFmodel;
        centerIndexHist(u,not(self_Terminate)) = centerIndexLoop(not(self_Terminate));
        for s=1:loaddimFlag
            if self_Terminate(s)==0
                center_daHist(u,:,s)=dainputs0(centerIndexLoop(s),:); %Variable stores the voltages of the RBF centers.
                %Dim 1= RBF #
                %Dim 2= Channel for voltage
                %Dim 3= Dimension center is placed in ( what load channel it is helping approximate)
            end
        end
                
        %Find and Store tares
        if FLAGS.tare_intercept==1 %If tare loads were included in regression
            [xcalib_RBF,taresGRBF]=RBF_tareCalc(xcalib_RBF,nterms_RBF,dainputs_zero,comIN_zero,epsHist(1:u,:),center_daHist(1:u,:,:),h_GRBF); %Calculate tares from series specific intercepts           
        else
            taresGRBF=zeros(nseries0,loaddimFlag); %Else set to zero (no series intercepts)
        end
        taretalRBF=taresGRBF(series0,:);
        tareGRBFHist{u} = taresGRBF;
        
        %update the approximation
        aprxIN2=comIN0_RBF*xcalib_RBF; %Approximation including series intercepts
        aprxIN2_Hist{u} = aprxIN2;
        aprxINminGZ2=aprxIN2+taretalRBF; %Approximation that does not include series intercept terms
        %Calculate tare corrected load approximation
        aprxINminTARE2=aprxINminGZ2-taretalRBF;
        
        %    QUESTION: JRP; IS THIS NECESSARY/USEFUL?
        % Find Standard Deviation of mean tares
        if FLAGS.tare_intercept==1 %If tare loads were included in regression
            [~,taresGRBF_STDDEV_all] = meantare(series0,aprxINminGZ2-targetMatrix0);
        else
            taresGRBF_STDDEV_all=zeros(size(targetMatrix0));
        end
        taresGRBFSTDEV = taresGRBF_STDDEV_all(s_1st0,:);
        
        %Extract RBF coefficients
        coeff_algRBFmodel=xcalib_RBF(1:nterms_RBF,:); %Algebraic and RBF coefficient matrix
        coeff_algRBFmodel_alg=xcalib_RBF(1:nterms,:); %new algebraic coefficients
        coeff_algRBFmodel_RBF_diag=xcalib_RBF(nterms+1:nterms_RBF,:); %new RBF coefficients, spaced on diagonals
        %Extract only RBF coefficients in compact matrix
        coeff_algRBFmodel_RBF=zeros(u,loaddimFlag);
        for i=1:u
            coeff_algRBFmodel_RBF(i,:)=diag(coeff_algRBFmodel_RBF_diag(1+loaddimFlag*(i-1):loaddimFlag*i,:));
        end        
        
        %Store basis parameters in Hist variables
        cHist_tot{u} = coeff_algRBFmodel;
        
        %Calculate and store residuals
        targetRes2 = targetMatrix0-aprxINminTARE2;
        newRes2 = targetRes2'*targetRes2;
        resSquare2 = diag(newRes2);
        resSquareHist(u,:) = resSquare2;
        resStdHist(u,:)=std(targetRes2);
        
        %Validation Error Self-Termination Check: use new ALG+RBF model to
        %determine standard deviation of residuals for validation data
        if FLAGS.valid_selfTerm==1
            %Test on validation data
            comINvalid_RBF(:,(u-1)*loaddimFlag+1:u*loaddimFlag)=create_comIN_RBF(dainputsvalid,epsHist(u,:),center_daHist(u,:,:),h_GRBF); %Generate comIN for RBFs
            comINvalid_algRBF=[comINvalid, comINvalid_RBF(:,1:u*loaddimFlag)]; %Combine comIN from algebraic terms and RBF terms to multiply by coefficients
            
            aprxINminGZ2valid=comINvalid_algRBF*coeff_algRBFmodel; %find approximation with alg and RBF Coefficients
            
            % SOLVE FOR TARES BY TAKING THE MEAN
            [~,s_1st,~] = unique(seriesvalid);
            if FLAGS.tare_intercept==1 %If tare loads were included in regression
                [taresAllPointsvalid2,taretalstdvalid2] = meantare(seriesvalid,aprxINminGZ2valid-targetMatrixvalid);
            else
                taresAllPointsvalid2=zeros(size(targetMatrixvalid));
                taretalstdvalid2 = zeros(size(targetMatrixvalid));
            end
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
            for i=1:loaddimFlag
                if period_change(u,i)>0 && self_Terminate(i)==0
                    fprintf(strcat('\n Channel'," ", string(i), ' Reached validation period change termination criteria, # RBF=',string(u)));
                    self_Terminate(i)=1;
                end
            end
        end
        
        %Prediction Error Self-Termination Check: caculate PI for
        %calibration data and store RMS of all calibration points
        if (FLAGS.PI_selfTerm==1 || FLAGS.VIF_selfTerm==1) && any(~self_Terminate) %If terminating based on PI or VIF and not all the channels have been self terminated
            [loadPI_GRBF_iter]=calc_PI(ANOVA_GRBF,anova_pct,comIN0_RBF(:,1:nterms_RBF),aprxIN2,calc_channel); %Calculate prediction interval for loads
            for i=1:loaddimFlag
                if self_Terminate(i)==1 %If channel is self terminated, use PI RMS from previous iteration
                    calib_PI_rms_Hist(u,i)=calib_PI_rms_Hist(u-1,i);
                else
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
            for i=1:loaddimFlag
                if period_change(u,i)>0 && self_Terminate(i)==0
                    fprintf(strcat('\n Channel'," ", string(i), ' Reached Prediction Interval period change termination criteria, # RBF=',string(u)));
                    self_Terminate(i)=1;
                    if all(self_Terminate) %If all the channels have now been self-terminated
                        calib_PI_rms_Hist(u+1:end,:)=[];  %Trim storage variable
                    end
                end
            end
        end
        
        if all(self_Terminate) %Check if all channels have self terminated
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
    final_RBFs_added=RBFs_added; %Initialize count of # RBFs/channel for final model
    
    %If validation self-termination selected, recalculate for RBF number of min
    %validation STD
    if FLAGS.valid_selfTerm==1
        fprintf(strcat('\n Trimming RBFs for minimum validation STD'));
        %Find RBF number for lowest Validation STD
        min_validSTD_num=zeros(1,loaddimFlag);
        for i=1:loaddimFlag
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
        min_calibPI_num=zeros(1,loaddimFlag);
        final_calibPI_rms=zeros(1,loaddimFlag);
        for i=1:loaddimFlag
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
        RBF_custom=repmat(eye(loaddimFlag,loaddimFlag),max(final_RBFs_added),1);
        for i=1:loaddimFlag
            RBF_custom(loaddimFlag*final_RBFs_added(i)+1:end,i)=0;
        end
        if FLAGS.model==4
            customMatrix_RBF=[customMatrix(1:nterms,:);RBF_custom;customMatrix(nterms+1:end,:)];
        else
            customMatrix_RBF=[ones(nterms,loaddimFlag);RBF_custom;ones(nseries0,loaddimFlag)];
        end
        
        %Trim comIN
        comIN0_RBF(:,nterms_RBF+1:nterms+u*loaddimFlag)=[];
        %Zero out parameters for trimmed RBFs in each channel
        for s=1:loaddimFlag 
            epsHist(final_RBFs_added(s)+1:end,s) = 0;
            centerIndexHist(final_RBFs_added(s)+1:end,s) = 0;
            center_daHist(final_RBFs_added(s)+1:end,:,s)=0; %Variable stores the voltages of the RBF centers.
            %Dim 1= RBF #
            %Dim 2= Channel for voltage
            %Dim 3= Dimension center is placed in ( what load channel it is helping approximate)
        end
        %Trim RBF property variables
        epsHist(max(final_RBFs_added)+1:end,:) = [];
        centerIndexHist(max(final_RBFs_added)+1:end,:) = [];
        center_daHist(max(final_RBFs_added)+1:end,:,:)=[];
             
        %New flag structure for calc_xcalib
        FLAGS_RBF=FLAGS; %Initialize as global flag structure
        FLAGS_RBF.model=4; %Calculate with custom model
        FLAGS_RBF.anova=FLAGS.anova; %Calculate ANOVA based on user preference
        FLAGS_RBF.test_FLAG=0; %Calculate VIF
        calc_channel=ones(1,loaddimFlag); %Calculate stats for every channel
        nterms_RBF=nterms+max(final_RBFs_added)*loaddimFlag; %New number of terms to solve for
        
        %Calculate Algebraic and RBF coefficients with calc_xcalib function
        [xcalib_RBF, ANOVA_GRBF] = calc_xcalib(comIN0_RBF,targetMatrix0,series0,...
            nterms_RBF,nseries0,loaddimFlag,FLAGS_RBF,customMatrix_RBF,anova_pct,loadlist,'Direct w RBF',calc_channel);
        
        %Find and Store tares
        if FLAGS.tare_intercept==1 %If tare loads were included in regression
            [xcalib_RBF,taresGRBF]=RBF_tareCalc(xcalib_RBF,nterms_RBF,dainputs_zero,comIN_zero,epsHist,center_daHist,h_GRBF); %Calculate tares from series specific intercepts
        else
            taresGRBF=zeros(nseries0,loaddimFlag); %Else set to zero (no series intercepts)
        end
        taretalRBF=taresGRBF(series0,:);
        tareGRBFHist{u+1} = taresGRBF;
        
        %update the approximation
        aprxIN2=comIN0_RBF*xcalib_RBF; %Approximation including series intercepts
        aprxIN2_Hist{u+1} = aprxIN2;
        aprxINminGZ2=aprxIN2+taretalRBF; %Approximation that does not include series intercept terms
        %Calculate tare corrected load approximation
        aprxINminTARE2=aprxINminGZ2-taretalRBF;
        
        %    QUESTION: JRP; IS THIS NECESSARY/USEFUL?
        % Find Standard Deviation of mean tares
        if FLAGS.tare_intercept==1 %If tare loads were included in regression
            [~,taresGRBF_STDDEV_all] = meantare(series0,aprxINminGZ2-targetMatrix0);
        else
            taresGRBF_STDDEV_all=zeros(size(targetMatrix0));
        end
        taresGRBFSTDEV = taresGRBF_STDDEV_all(s_1st0,:);
                
        %Extract RBF coefficients
        coeff_algRBFmodel=xcalib_RBF(1:nterms_RBF,:); %Algebraic and RBF coefficient matrix
        coeff_algRBFmodel_alg=xcalib_RBF(1:nterms,:); %new algebraic coefficients
        coeff_algRBFmodel_RBF_diag=xcalib_RBF(nterms+1:nterms_RBF,:); %new RBF coefficients, spaced on diagonals
        %Extract only RBF coefficients in compact matrix
        coeff_algRBFmodel_RBF=zeros(max(final_RBFs_added),loaddimFlag);
        for i=1:max(final_RBFs_added)
            coeff_algRBFmodel_RBF(i,:)=diag(coeff_algRBFmodel_RBF_diag(1+loaddimFlag*(i-1):loaddimFlag*i,:));
        end
        
        %Update basis parameters in Hist variables
        cHist_tot{u+1} = coeff_algRBFmodel;

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
            'nterms',nterms+max(final_RBFs_added)*loaddimFlag,...
            'coeff_algRBFmodel',coeff_algRBFmodel,...
            'coeff',coeff);
        uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
            [fieldnames(uniqueOut); fieldnames(newStruct)],1);
        
        if FLAGS.mode==1
        newStruct=struct('loadCapacities',loadCapacities,'tares',taresGRBF,'tares_STDDEV',taresGRBFSTDEV);
                uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
            [fieldnames(uniqueOut); fieldnames(newStruct)],1);
        end
        
        output(section,FLAGS,targetRes2,fileName,numpts0,nseries0,...
            loadlist,series0,excessVec0,voltdimFlag,loaddimFlag,voltagelist,...
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
        if FLAGS.tare_intercept==1 %If tare loads were included in regression
            [taresAllPointsvalid2,taretalstdvalid2] = meantare(seriesvalid,aprxINminGZ2valid-targetMatrixvalid);
        else
            taresAllPointsvalid2=zeros(size(targetMatrixvalid));
            taretalstdvalid2=zeros(size(targetMatrixvalid));
        end
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
            
            %Store Variables for output
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
        
        if FLAGS.mode==1
            newStruct=struct('loadCapacitiesvalid',loadCapacitiesvalid,'tares',taresGRBFvalid,'tares_STDDEV',taresGRBFSTDEVvalid);
            uniqueOut = cell2struct([struct2cell(uniqueOut); struct2cell(newStruct)],...
                [fieldnames(uniqueOut); fieldnames(newStruct)],1);
        end
            
            output(section,FLAGS,targetRes2valid,fileNamevalid,numptsvalid,nseriesvalid,...
                loadlist,seriesvalid,excessVecvalid,voltdimFlagvalid,loaddimFlagvalid,voltagelist,...
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
    if FLAGS.mode~=1
        seriesapprox=ones(size(excessVecapprox,1),1);
        series2approx=ones(size(excessVecapprox,1),1);
        natzerosapprox=0;
    end
    
    if voltdimFlag~=size(excessVecapprox,2) %Check if mismatch between calibration and approximation data.  If so, exit program
        fprintf('\n  ');
        fprintf('\n MISMATCH IN CALIBRATION/APPROXIMATION DATA DIMENSIONS.  UNABLE TO PROCEED.\n');
        fprintf('\n');
        if isdeployed % Optional, use if you want the non-deployed version to not exit immediately
            input('Press enter to finish and close');
        end
        return; %Quit run
    end
    
    if exist( 'pointIDvalid', 'var')==0
        pointIDapprox=cellstr([repmat('P-',size(excessVecapprox,1),1),num2str((1:size(excessVecapprox,1))')]);
    end
    
    if FLAGS.balCal == 2 %If RBFs were placed, put parameters in structure
        GRBF.epsHist=epsHist;
        GRBF.coeff_algRBFmodel=coeff_algRBFmodel;
        GRBF.center_daHist=center_daHist;
        GRBF.h_GRBF=h_GRBF;
        GRBF.ANOVA=ANOVA_GRBF;
    else
        GRBF='GRBFS NOT PLACED';
    end
    
    %Function that performs all Approximation calculations and outputs
    [aprxINminGZapprox,loadPI_approx]=AOX_approx_funct(coeff,natzerosapprox,excessVecapprox,FLAGS,seriesapprox,...
        series2approx,pointIDapprox,loadlist,file_output_location,GRBF,ANOVA,anova_pct);
    
end
%END APPROXIMATION SECTION

%File Cleanup
if FLAGS.input_save==0  %If user did not want to save input data files, delete
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
if isdeployed % Optional, use if you want the non-deployed version to not exit immediately
    input('Press enter to finish and close');
end





function [xcalib_RBF,taresGRBF]=RBF_tareCalc(xcalib_RBF,nterms_RBF,dainputs_zero,comIN_zero,epsHist,center_daHist,h_GRBF)
%Function calculates tare loads for load model including RBFs. This is
%accomplished by calculating all coefficients simultaneously including
%series specific intercepts.  The tares are then extracted from the series
%specific intercepts.

%Each series intercept includes 2 components: 1 portion shifts
%the surface to the 'reality' of 0 load=0 voltage.  The second
%portion provides a shift for the tare loads. The 'reality
%shift' is a global intercept that must be applied to each
%series.  The second portion of the shift is different in each
%series based on the tare load applied. Therefore, by finding
%the shift required for our RBF surface alone to match the
%condition of 0 load at 0 voltage, we can split the shift into
%its 2 portions.

%INPUTS:
%  xcalib_RBF = Coefficient matrix including series specific intercepts
%  nterms_RBF  =  Number of predictor variable terms including RBFs
%  dainputs_zero  =  Vector of zero voltages
%  comIN_zero  =  Matrix of algebraic predictor variables at zero voltage
%  epsHist  =  Epsilon (width control) values for RBFs
%  center_daHist = Center locations for RBFs
%  h_GRBF = h values for RBFs (controls width)

%OUTPUTS:
%  y = Merit Value for success in RBF fitting residuals.  Object of optimization is to minimize y

seriesShift = xcalib_RBF(nterms_RBF+1:end,:);
comIN_zero_RBF=create_comIN_RBF(dainputs_zero,epsHist,center_daHist,h_GRBF); %Generate comIN for RBFs at zero voltage
comIN_zero_algRBF=[comIN_zero, comIN_zero_RBF]; %Combine comIN from algebraic terms and RBF terms to multiply by coefficients
currentZero=comIN_zero_algRBF*xcalib_RBF(1:nterms_RBF,:); %Current load predicted for zero voltage without series shift
realityShift=-currentZero; %Extract portion of each series shift that shifts to the reality of 0 voltage=0 voltage
tareShift=seriesShift-realityShift; %Remainder of series shift accounts for tares
taresGRBF=-tareShift; %Negative of tare shifts are tare loads

%Update xcalib by splitting intercepts into global reality
%shift and series specific tare shift
xcalib_RBF(1,:)=realityShift; %Include global intercept for reality shift
xcalib_RBF(nterms_RBF+1:end,:)=tareShift; %Series specific intercepts are for tares
end



