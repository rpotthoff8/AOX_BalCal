% Copyright 2018 Andrew Meade, Ali Arya Mokhtarzadeh and Javier Villarreal.  All Rights Reserved.
%
%balanceCalibration_with_RBF_8D.m
%requires "balCal_meritFunction.m" to run
%input file: "BuffetBalance-CalDataOfOct2015-3F3M.csv"
%output file: "balCal_output_ALGB.xls"
%output file: "balCal_output_GRBF.xls"
%
%%
%initialize the workspace
clc;
clearvars;
close all;
workspace;
disp('Copyright 2018 Andrew Meade, Ali Arya Mokhtarzadeh and Javier Villarreal.  All Rights Reserved.')
% The mean of the approximation residual (testmatrix minus local approximation) for each section is taken as the tare for that channel and section. The tare is subtracted from the global values to make the local loads. The accuracy of the validation, when compared to the known loads, is very sensitive to the value of the tares (which is unknown) and NOT the order of the calibration equations.
% Because of measurement noise in the voltage the APPROXIMATION tare is computed by post-processing. The average and stddev is taken of all channels per section. If the stddev is less than 0.25% of the capacity for any station the tare is equal to the average for that channel. If the stddev is greater than 0.25% then the approximation at the local zero is taken as the tare for that channel and section. The global approximation is left alone but the tare is subtracted from the values to make the local loads. Line 3133.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       USER INPUT SECTION
%
out = AOX_GUI;
if out.cancel == 1
    return
end
%TO SELECT Algebraic Model                         set balCal_FLAG = 1;
%TO SELECT Algebraic and GRBF Model                set balCal_FLAG = 2;
balCal_FLAG = out.grbf;
%
%TO SELECT INDIRECT APPROACH         set approach_FLAG = 1;
approach_FLAG = 0;%out.approach;
%
%SELECT ALGEBRAIC MODEL     set model_FLAG = 1 (full), 2 (trun), or 3 (linear)
model_FLAG = out.model;
%
%TO OUTPUT LOAD PERFORMANCE PARAMETERS TO SCREEN        set print_FLAG = 1;
%TO OUTPUT INTERMEDIATE PERFORMANCE PARAMETERS TO SCREEN set print_FLAG = 2;
print_FLAG = out.tables;
%
%TO OUTPUT DATA TO EXCEL (CSV)                          set excel_FLAG = 1;
excel_FLAG = out.excel;
%
%TO OUTPUT CORRELATION PLOTS                       set corr_FLAG = 1;
corr_FLAG = out.corr;
%
%TO OUTPUT CORRELATION PLOTS                       set rescorr_FLAG = 1;
rescorr_FLAG = out.rescorr;
%
%TO OUTPUT RESIDUALS                               set rest_FLAG = 1;
res_FLAG = out.res;
%
%TO OUTPUT HISTOGRAMS                              set hist_FLAG = 1;
hist_FLAG = out.hist;
%
%TO SELECT Validation of the Model                 set balVal_FLAG = 1;
balVal_FLAG = out.valid;
%
%TO SELECT Approximation from Cal Data             set balApprox_FLAG = 1;
balApprox_FLAG = out.approx;
%
%DEFINE THE NUMBER OF BASIS FUNCTIONS
numBasis = out.basis;
%
%TO IDENTIFY OUTLIERS                              set balOut_FLAG = 1;
balOut_FLAG = out.outlier;
numSTD = out.numSTD;          %Number of standard deviations for outlier threshold.
%
%TO USE THE INPUT-OUTPUT MATRICES WITHOUT OUTLIERS   set zeroed_FLAG = 1;
zeroed_FLAG = out.zeroed;
%
%TO USE LATIN HYPERCUBE SAMPLING
LHS_Flag = out.lhs;
numLHS = out.numLHS; %Number of times to iterate.
LHSp = out.LHSp; %Percent of data used to create sample.
% A mean of the coefficients is calculated afterwards.
%
%                       END USER INPUT SECTION
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(out.savePathcal,'-mat');
series0 = series;
% Loads:
% loadlabes, voltlabels (if they exist)
% loadCapacities, natzeros, targetMatrix0, excessVec0, series0

disp('Starting Calculations')
% if approach_FLAG == 1
%     disp(' ');
%     disp('Using the Indirect Approach for Calibration');
% else
%     disp(' ');
%     disp('Using the Direct Approach for Calibration');
% end

if exist('loadlabels','var')
    loadlist = loadlabels;
    voltagelist = voltlabels;
else
    loadlist = {'NF','BM','S1','S2','RM','AF','PLM', 'PCM', 'MLM', 'MCM'};
    voltagelist = {'rNF','rBM','rS1','rS2','rRM','rAF','rPLM','rPCM','rMLM','rMCM'};
end

if corr_FLAG == 1
    figure('Name','Correlation plot','NumberTitle','off');
    correlationPlot(targetMatrix0, excessVec0, loadlist, voltagelist);
end

%find the average natural zeros (also called global zeros)
globalZeros = mean(natzeros);

[~,s_1st0,s_id0] = unique(series0);
nseries0 = length(s_1st0);
[numpts0, dimFlag] = size(excessVec0);

[localZeros0,localZerosAllPoints0] = localzeros(series0,excessVec0);
globalZerosAllPoints0 = ones(numpts0,1)*globalZeros;

dainputs0 = excessVec0 - globalZerosAllPoints0;

switch model_FLAG
    case {1,4}
        % Full Algebraic Model
        nterms = 2*dimFlag*(dimFlag+2);
    case 2
        % Truncated Algebraic Model
        nterms = dimFlag*(dimFlag+3)/2;
    case 3
        % Linear Algebraic Model
        nterms = dimFlag;
end

%% Temporary; will be user input in the future
if 1
    model_FLAG = 4;
    nterms = 2*dimFlag*(dimFlag+2);
    customMatrix = randi(2,nterms,dimFlag)-1;
end
%%
if model_FLAG == 4
    customMatrix = [customMatrix; ones(nseries0,dimFlag)];
end

comIN0 = balCal_algEqns(model_FLAG,dainputs0,series0);

if LHS_Flag == 0
    numLHS = 1;
end
lhs_check = zeros(length(excessVec0),1);
ind = 1:length(excessVec0(:,1));

if LHS_Flag == 1
    disp('  ')
    disp('Number of LHS Iterations Selected =')
    disp(numLHS);
end

for lhs = 1:numLHS
    
    if LHS_Flag == 1
        sample = AOX_LHS(series0,excessVec0,LHSp);
        lhs_check(sample) = 1;
        ind(find(lhs_check-1)); % This line outputs which data points haven't been sampled yet
        pct_sampled = sum(lhs_check)/length(lhs_check); % This line outputs what percentage of points have been sampled
    else
        sample = [1:length(series0)]';
    end
    
    series = series0(sample);
    targetMatrix = targetMatrix0(sample,:);
    comIN = comIN0(sample,:);
    
    [~,s_1st,~] = unique(series);
    nseries = length(s_1st);
    
    disp('  ')
    disp('Working ...')
    % Indirect approach uses the modeled voltage
%     if approach_FLAG == 1
%         if balCal_FLAG == 2
%             excessVec0 = qtaprxINminGZ2 + globalZerosAllPoints;
%         else
%             excessVec0 = qtaprxINminGZ + globalZerosAllPoints;
%         end
%     end  

    xcalib = zeros(nterms+nseries0,dimFlag);
    for k = 1:dimFlag
        comIN_k = comIN;
        
        if model_FLAG == 4
            comIN_k(:,customMatrix(:,k)==0) = [];
        end
        
        %SOLUTION
        xcalib_k = pinv(comIN_k)*targetMatrix(:,k);
        
        if model_FLAG == 4
            xcalib(customMatrix(:,k)==1,k) = xcalib_k;
        else
            xcalib(:,k) = xcalib_k;
        end
        
    end
    
    if LHS_Flag == 1
        x_all(:,:,lhs) = xcalib;
    end
end

if LHS_Flag == 1
    xcalib = mean(x_all,3);
    xcalib_std = std(x_all,[],3);
end
coeff = xcalib(1:nterms,:);
tares = -xcalib(nterms+1:end,:);
taresAllPoints = tares(s_id0,:);

%  Creates Matrix for the volts to loads
APPROX_AOX_COEFF_MATRIX = coeff;
xapprox = coeff;
if excel_FLAG == 1
    filename = 'APPROX_AOX_COEFF_MATRIX.csv';
    csvwrite(filename,xapprox)
end

%APPROXIMATION
%define the approximation for inputs minus global zeros
aprxIN = comIN0*xcalib;

%RESIDUAL
targetRes = targetMatrix0-aprxIN;      %0=b-Ax

reslist = {'resNF','resBM','resS1','resS2','resRM','resAF','resPLM',...
    'resPCM', 'resMLM', 'resMCM'};
if rescorr_FLAG == 1
    figure('Name','Residual correlation plot','NumberTitle','off');
    correlationPlot(excessVec0, targetRes, voltagelist, reslist);
end

%find the sum of squares of the residual using the dot product
resSquare = dot(targetRes,targetRes)';
%AAM note to self - in matlab, diag(A'*A) is the same as dot(A,A)'

% Identify Outliers After Filtering
% (Threshold approach) ajm 8/2/17
if balOut_FLAG == 1
    detect_targetRes = targetRes;
    
    % Use the modeled input for the rest of the calculations
    for n = 1:dimFlag
        normtargetRes(:,n) = detect_targetRes(:,n)/loadCapacities(n);
    end
    out_meanValue = mean(normtargetRes);
    
    % Identify outliers. They are considered outliers if the residual
    % is more than 3 standard deviations as % of capacity from the mean.
    out_standardDev = std(normtargetRes);
    numSTD = 3.0; % Whatever you want.
    thresholdValue = numSTD * (out_standardDev) - out_meanValue;
    
    for n = 1:dimFlag
        if thresholdValue(1,n) <= 0.0025
            thresholdValue(1,n) = 0.0025;
        end
    end
    
    outlierIndices = abs(normtargetRes) > thresholdValue;
    
    % ID outlier rows :
    zero_counter = 1;
    for k1 = 1:numpts
        for k4 = 1:dimFlag
            if outlierIndices(k1,k4) == 1
                outlier_values(zero_counter,1) = k1;
                zero_counter = zero_counter + 1;
            end
        end
    end
    
    OUTLIER_ROWS = unique(outlier_values,'rows');
    
    num_outliers = length(OUTLIER_ROWS);
    prcnt_outliers = 100.0*num_outliers/numpts;
    
    % Use the reduced input and target files
    if zeroed_FLAG == 1
        disp(' ************************************************************************ ');
        disp('Find the reduced data in zeroed_targetMatrix and zeroed_excessVec');
        disp(' ************************************************************************ ');
        
        % Mark the outliers values, store and use for recalculation:
        zeroed_targetMatrix = targetMatrix0;
        zeroed_excessVec = excessVec0;
        zeroed_series = series0;
        zeroed_numpts = numpts;
        
        zeroed_targetMatrix(OUTLIER_ROWS,:) = [];
        zeroed_excessVec(OUTLIER_ROWS,:) = [];
        
        targetMatrix0 = unique(zeroed_targetMatrix,'rows');
        excessVec0 = unique(zeroed_excessVec,'rows');
        
        zeroed_series(OUTLIER_ROWS) = [];
        
        numpts =  zeroed_numpts - num_outliers;
        
        targetMatrix0 = zeroed_targetMatrix;
        excessVec0 = zeroed_excessVec;
        series0 = zeroed_series;
        
        dainputs = zeros(numpts,dimFlag);
        dalz = zeros(numpts,dimFlag);
        localZerosAllPoints = zeros(numpts,dimFlag);
        aprxINminGZ = zeros(numpts,dimFlag);
        aprxLZminGZ = zeros(numpts,dimFlag);
        zeroed_checkit = zeros(numpts,dimFlag);
        taresAllPoints = zeros(numpts,dimFlag);
        globalZerosAllPoints = zeros(numpts,dimFlag);
        eta = zeros(numpts,dimFlag);
        zeroed_zoop = zeros(numpts,dimFlag);
        
        rbfc_INminLZ = zeros(numpts,dimFlag);
        rbfc_INminGZ = zeros(numpts,dimFlag);
        rbfc_LZminGZ = zeros(numpts,dimFlag);
        rbfINminLZ = zeros(numpts,dimFlag);
        rbfINminGZ = zeros(numpts,dimFlag);
        rbfLZminGZ = zeros(numpts,dimFlag);
        
        [localZeros,localZerosAllPoints] = localzeros(series0,excessVec0);
        globalZerosAllPoints = ones(numpts,1)*globalZeros;
        
        % Subtract the Global Zeros from the Inputs and Local Zeros
        for i=1:dimFlag
            dainputs(:,i) = excessVec0(:,i) - globalZerosAllPoints(i);
            dalz(:,i) = localZerosAllPoints(:,i) - globalZerosAllPoints(i);
        end
        
        for i=1:dimFlag
            biggee(:,i) = 0;
        end
        [comIN,comLZ]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numpts,series0,nseries0, dainputs, dalz, biggee);
        
        comINminLZ = comIN-comLZ;
        
        %SOLUTION
        xcalib = pinv(comINminLZ')*targetMatrix0;
        
        %APPROXIMATION
        %define the approximation for inputs minus local zeros
        aprxINminLZ = comINminLZ'*xcalib;
        
        %A DIFFERENT APPROXIMATION
        %define the approximation for inputs minus global zeros
        intercepts = -(comGZ'*xcalib);
        aprxIN = (xcalib'*comIN)';
        aprxLZ = (xcalib'*comLZ)';       %to find tares AAM042016
        for m=1:length(aprxIN)
            aprxINminGZ(m,:) = aprxIN(m,:);
            aprxLZminGZ(m,:) = aprxLZ(m,:);%to find tares
            
            %subtracts the targets from the global approx
            %This will be averaged over the individual series0 to give the tares
            checkit(m,:) = aprxINminGZ(m,:)-targetMatrix0(m,:);
        end
        
        taresAllPoints = meantare(series0,checkit);
        %RESIDUAL
        targetRes = targetMatrix0+taresAllPoints-aprxINminGZ;      %0=b-Ax
    end
end

%find the sum of squares of the residual using the dot product
resSquare = dot(targetRes,targetRes)';
%AAM note to self - in matlab, diag(A'*A) is the same as dot(A,A)'

%OUTPUTS FOR ALGEBRAIC SECTION
for k=1:length(targetRes(1,:))
    [goop(k),kstar(k)] = max(abs(targetRes(:,k)));
    goopVal(k) = abs(targetRes(kstar(k),k));
    xCent(k) = excessVec0(kstar(k),k);
    maxTargets(k) = max(targetRes(:,k));
    minTargets(k) = min(targetRes(:,k));
end
perGoop = 100*(goop./loadCapacities);
davariance = var(targetRes);
gee = mean(targetRes);
standardDev10 = std(targetRes);
standardDev = standardDev10';
stdDevPercentCapacity = 100*(standardDev'./loadCapacities);
ratioGoop = goop./standardDev';
ratioGoop(isnan(ratioGoop)) = realmin;

%    theminmaxband = abs(maxTargets + minTargets);
theminmaxband = 100*(abs(maxTargets + minTargets)./loadCapacities);

[~,s_1st,~] = unique(series0);
taresALGB = taresAllPoints(s_1st,:);

%OUTPUT HISTOGRAM PLOTS
if hist_FLAG == 1
    figure('Name','Calibration - ALGB','NumberTitle','off')
    for k0=1:length(targetRes(1,:))
        subplot(2,3,k0)
        binWidth = 0.25;
        edges = [-4.125:binWidth:4.125];
        h = histogram(targetRes(:,k0)/standardDev(k0,:),edges,'Normalization','probability');
        centers = edges(1:end-1)+.125;
        values = h.Values*100;
        bar(centers,values,'barwidth',1)
        ylabel('% Data Pts');
        xlim([-4 4]);
        ylim([0 50]);
        hold on
        plot(linspace(-4,4,100),binWidth*100*normpdf(linspace(-4,4,100),0,1),'r')
        hold off
        xlabel(['\Delta',loadlist{k0},'/\sigma']);
    end
end

%START PRINT OUT PERFORMANCE INFORMATION TO THE SCREEN
if print_FLAG == 1
    %
    %% END Direct or Indirect Approach to Calibration
    
    %% Identify the Possible Outliers
    if balOut_FLAG == 1
        disp(' ***** ');
        disp(' ');
        disp('Number of Outliers =');
        disp(num_outliers);
        disp('Outliers % of Data =');
        disp(prcnt_outliers);
    end
    
    % Recalculated Calibration with Reduced Matrices
    if zeroed_FLAG == 1
        disp(' ************************************************************************ ');
        disp('Find the reduced data in zeroed_targetMatrix and zeroed_excessVec');
        disp(' ************************************************************************ ');
    end
    
    
    %%%%%%% 6_14_18 ajm
    calib_twoSigmaALGB = standardDev'.*2;
    calib_algebraic_2Sigma = array2table(calib_twoSigmaALGB,'VariableNames',loadlist(1:dimFlag))
    
    %Should I use strtrim()  ? -AAM 042116
    calib_algebraic_Tares = array2table(taresALGB,'VariableNames',loadlist(1:dimFlag))
    
    coefficientsALGB = [xcalib;intercepts];
    for coeffCount=1:size(coefficientsALGB(:,1),1)
        numCombin(coeffCount) = coeffCount;
    end
    numCombinName = cellstr(num2str(numCombin'));
    numCombinName(nterms+nseries0+1) = cellstr('Intercept');
    
    calib_mean_algebraic_Resids_sqrd = array2table(resSquare'./numpts,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Pcnt_Capacity_Max_Mag_Load_Resids = array2table(perGoop,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Std_Dev_pcnt = array2table(stdDevPercentCapacity,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Max_Load_Resids = array2table(maxTargets,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Min_Load_Resids = array2table(minTargets,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Ratio_Max_Mag_Load_Resid_and_Std_Dev = array2table(ratioGoop,'VariableNames',loadlist(1:dimFlag))
    
    % Prints the minmaxband
    calib_alg_per_minmaxband = array2table(theminmaxband,'VariableNames',loadlist(1:dimFlag))
    %%%%%%%%%
    
end

if excel_FLAG == 1
    % Output results to an excel file
    disp('  ');
    disp('ALG CALIBRATION MODEL GLOBAL LOAD APPROXIMATION FILE: CALIB_AOX_GLOBAL_ALG_RESULT.csv');
    % CALIB_AOX_GLOBAL_ALG_RESULT = aprxINminGZ;
    disp(' ');
    
    filename = 'CALIB_AOX_GLOBAL_ALG_RESULT.csv';
    csvwrite(filename,aprxINminGZ)
end

if res_FLAG == 1
    figure('Name','Algebraic Model Calibration; Residuals of Load Versus Data Point Index','NumberTitle','off')
    plotResPages(series0, targetRes, loadCapacities, stdDevPercentCapacity, loadlist)
    %    hold off
end

if balCal_FLAG == 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             VOLTAGE TO LOAD (DIRECT) - RBF SECTION                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %goal to minimize: minimize the sum of the squares (dot product) of each of the 8
    %residual vectors 'targetRes' 'target1' ... 'target8'
    %dt1 = dot(target1,target1);
    %find centers by finding the index of max residual, using that index to
    %subtract excess(counter)-excess(indexMaxResid) and then taking the dot
    %product of the resulting column vector
    for i=1:dimFlag
        for s=1:length(series0)
            targetRes2(s,i) = targetRes(s,i);
        end
    end
    aprxINminGZ2 = aprxINminGZ;
    
    etaHist = cell(numBasis,1);
    aprxINminGZ_Hist = cell(numBasis,1);
    tareHist = cell(numBasis,1);
    
    for i=1:dimFlag
        dainputscalib(:,i) = excessVec0(:,i)-globalZeros(i);
        dalzcalib(:,i) = localZerosAllPoints(:,i)-globalZeros(i);
    end
    
    %    localZeroMatrix = localZerosAllPoints;
    globalZerosAllPoints = zeros(length(excessVec0(:,1)),dimFlag); % ajm 6_2_18
    
    etaLZ = dot(dalzcalib-dainputscalib,dalzcalib-dainputscalib);
    etaGZ = dot(globalZerosAllPoints-dainputscalib,globalZerosAllPoints-dainputscalib);
    
    for u=1:numBasis
        for s=1:dimFlag
            [goopLoop(s),centerIndexLoop(s)] = max(abs(targetRes2(:,s)));
            
            for r=1:length(excessVec0(:,1))
                eta(r,s) = dot(dainputscalib(r,:)-dainputscalib(centerIndexLoop(s),:),dainputscalib(r,:)-dainputscalib(centerIndexLoop(s),:));
            end
            
            %find widths 'w' by optimization routine
            w(s) = fminbnd(@(w) balCal_meritFunction2(w,targetRes2(:,s),eta(:,s)),0,1 );
            
            rbfINminLZ(:,s)=exp(eta(:,s)*log(abs(w(s)))) - exp(etaLZ(:,s)*log(abs(w(s))));
            rbfINminGZ(:,s)=exp(eta(:,s)*log(abs(w(s))));
            rbfLZminGZ(:,s)=exp(etaLZ(:,s)*log(abs(w(s))));%to find tares AAM042016
            
            coeff(s) = dot(rbfINminGZ(:,s),targetRes2(:,s)) / dot(rbfINminGZ(:,s),rbfINminGZ(:,s));
            
            rbfc_INminLZ(:,s) = coeff(s)*rbfINminLZ(:,s);
            rbfc_INminGZ(:,s) = coeff(s)*rbfINminGZ(:,s);
            rbfc_LZminGZ(:,s) = coeff(s)*rbfLZminGZ(:,s); %to find tares AAM042016
        end
        
        wHist(u,:) = w;
        cHist(u,:) = coeff;
        centerIndexHist(u,:) = centerIndexLoop;
        etaHist{u} = eta;
        
        %update the approximation
        aprxINminGZ2 = aprxINminGZ2+rbfc_INminGZ;
        aprxINminGZ_Hist{u} = aprxINminGZ2;
        
        % SOLVE FOR TARES BY TAKING THE MEAN
        taresAllPointsGRBF = meantare(series0,aprxINminGZ2-targetMatrix0);
        
        taresGRBF = taresAllPointsGRBF(s_1st,:);
        
        tareGRBFHist{u} = taresGRBF;
        
        targetRes2 = targetMatrix0-aprxINminGZ2+taresAllPointsGRBF;      %0=b-Ax
        newRes2 = targetRes2'*targetRes2;
        resSquare2 = diag(newRes2);
        resSquareHist(u,:) = resSquare2;
    end
    
    for k2=1:length(targetRes(1,:))
        [goop2(k2),kstar2(k2)] = max(abs(targetRes2(:,k2)));
        goopVal2(k2) = abs(targetRes2(kstar2(k2),k2));
        xCent2(k2) = dainputscalib(kstar2(k2),k2);
        maxTargets2(k2) = max(targetRes2(:,k2));
        minTargets2(k2) = min(targetRes2(:,k2));
    end
    perGoop2 = 100*(goop2./loadCapacities);
    standardDev20 = std(targetRes2);
    standardDev2 = standardDev20';
    stdDevPercentCapacity2 = 100*(standardDev2'./loadCapacities);
    ratioGoop2 = goop2./standardDev2';
    ratioGoop2(isnan(ratioGoop2)) = realmin;
    
    theminmaxband2 = 100*(abs(maxTargets2 + minTargets2)./loadCapacities);
    
    %***************** ajm 5/12/18
    if res_FLAG == 1
        figure('Name','GRBF + Algebraic Model Calibration; Residuals of Load Versus Data Point Index','NumberTitle','off')
        plotResPages(series0, targetRes2, loadCapacities, stdDevPercentCapacity2, loadlist)
        %    hold off
    end
    
    %%% Diagnostics %%%%% 5/16/18
    %
    %justtherbfs = aprxINminGZ2-aprxINminGZ;
    %
    %if res_FLAG == 1
    %    figure('Name','Looking at GRBF Distribution in Calibration','NumberTitle','off')
    %    plotResPages(series0, justtherbfs, loadCapacities, stdDevPercentCapacity2)
    %    hold off
    %end
    
    %OUTPUT HISTOGRAM PLOTS
    if hist_FLAG == 1
        figure('Name','Calibration - GRBF','NumberTitle','off')
        for k0=1:length(targetRes2(1,:))
            subplot(2,3,k0)
            binWidth = 0.25;
            edges = [-4.125:binWidth:4.125];
            h = histogram(targetRes2(:,k0)/standardDev2(k0,:),edges,'Normalization','probability');
            centers = edges(1:end-1)+.125;
            values = h.Values*100;
            bar(centers,values,'barwidth',1)
            ylabel('% Data Pts');
            xlim([-4 4]);
            ylim([0 50]);
            hold on
            plot(linspace(-4,4,100),binWidth*100*normpdf(linspace(-4,4,100),0,1),'r')
            hold off
            xlabel(['\Delta',loadlist{k0},'/\sigma']);
        end
    end
    
    if excel_FLAG == 1
        disp(' ***** ');
        disp('  ');
        disp('ALG+GRBF CALIBRATION MODEL GLOBAL LOAD APPROXIMATION: Check CALIB_AOX_GLOBAL_GRBF_RESULT.csv file');
        
        filename = 'CALIB_AOX_GLOBAL_GRBF_RESULT.csv';
        csvwrite(filename,aprxINminGZ2)
        
        filename = 'APPROX_AOX_GRBF_ws.csv';
        csvwrite(filename,wHist)
        
        filename = 'APPROX_AOX_GRBF_coeffs.csv';
        csvwrite(filename,cHist)
        
        filename = 'APPROX_AOX_GRBF_Centers.csv';
        csvwrite(filename,centerIndexHist)
    end
    
    
    if print_FLAG == 1
        
        disp(' ');
        disp('Number of GRBFs =');
        disp(numBasis);
        disp(' ');
        
        twoSigmaGRBF = standardDev'.*2;
        calib_GRBF_2Sigma = array2table(twoSigmaGRBF,'VariableNames',loadlist(1:dimFlag))
        
        %Should I use strtrim()  ? -AAM 042116
        calib_GRBF_Tares = array2table(taresGRBF,'VariableNames',loadlist(1:dimFlag))
        
        calib_mean_GRBF_Resids_sqrd = array2table(resSquare2'./numpts,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Pcnt_Capacity_Max_Mag_Load_Resids = array2table(perGoop2,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Std_Dev_pcnt = array2table(stdDevPercentCapacity2,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Max_Load_Resids = array2table(maxTargets2,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Min_Load_Resids = array2table(minTargets2,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Ratio_Max_Mag_Load_Resid_and_Std_Dev = array2table(ratioGoop2,'VariableNames',loadlist(1:dimFlag))
        
        % Prints the GRBF minmax
        calib_GRBF_minmaxband_per_capacity = array2table(theminmaxband2,'VariableNames',loadlist(1:dimFlag))
    end
    
end

if balVal_FLAG == 1
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       VALIDATION SECTION      AJM 7/1/17                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %DEFINE THE VALIDATION CSV INPUT FILE AND SELECT THE RANGE OF DATA VALUES TO READ
    %     inputFile_balCal = 'MK14C-ChkLds-Ames2011-Meade-8D.csv';
    %     loadCapacitiesvalid =    csvread(inputFile_balCal,9,4,'E10..L10');
    %     natzerosvalid =          csvread(inputFile_balCal,12,12,'M13..T16');
    %     seriesvalid =            csvread(inputFile_balCal,19,2,'C20..C139');
    %     targetMatrixvalid =      csvread(inputFile_balCal,19,4,'E20..L139');
    %     excessVecvalid =         csvread(inputFile_balCal,19,12,'M20..T139');
    load(out.savePathval,'-mat');
    [~,s_1st,s_id] = unique(seriesvalid);
    xvalid = xcalib;
    
    excessVecvalid0 = excessVecvalid;
    % num of data points
    numptsvalid = length(seriesvalid);
    dimFlagvalid = length(excessVecvalid(1,:));
    
    %find the average natural zeros (also called global zeros)
    globalZerosvalid = mean(natzerosvalid);
    
    %load capacities
    loadCapacitiesvalid(loadCapacitiesvalid == 0) = realmin;
    
    %find number of series0; this will tell us the number of tares
    nseriesvalid = max(seriesvalid);
    
    %find zero points of each series0 and number of points in a series0
    %localZerosAllPoints is the same as localZeroMatrix defined in the RBF
    %section
    [localZerosvalid,localZerosAllPointsvalid] = localzeros(seriesvalid,excessVecvalid0);
    globalZerosAllPointsvalid = ones(numptsvalid,1)*globalZerosvalid;
    
    % Subtract the Global Zeros from the Inputs and Local Zeros
    
    for k=1:dimFlagvalid
        dainputsvalid(:,k) = excessVecvalid0(:,k)-globalZerosvalid(k);
        dalzvalid(:,k) = localZerosAllPointsvalid(:,k)-globalZerosvalid(k);
        dagzvalid(:,k) = 0;
    end
    
    %%% 5/16/18
    %Remember that  excessVec0 = excessVec0_complete - globalZerosAllPoints;
    excessVecvalidkeep = excessVecvalid0  - globalZerosAllPointsvalid;
    %%%
    
    % Build the Algebraic Model
    lasttarevalid = seriesvalid(numptsvalid);
    % Full Algebraic Model
    if model_FLAG == 1
        nterms = 2*dimFlag*(dimFlag+2);
    end
    % Truncated Algebraic Model
    if model_FLAG == 2
        nterms = dimFlag*(dimFlag+3)/2;
    end
    % Linear Algebraic Model
    if model_FLAG == 3
        nterms = dimFlag;
    end
    
    % Call the Algebraic Subroutine
    comGZvalid = zeros(nterms+1,1);
    [comINvalid,comLZvalid,comGZvalid]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numptsvalid,seriesvalid,nseriesvalid,dainputsvalid,dalzvalid,dagzvalid);
    
    comINminLZvalid = comINvalid-comLZvalid;
    
    %VALIDATION APPROXIMATION
    %define the approximation for inputs minus global zeros
    %    interceptsvalid = -(comGZvalid'*xvalid);  % ajm 5/17/18
    aprxINvalid = (xvalid'*comINvalid)';        %to find approximation AJM111516
    aprxLZvalid = (xvalid'*comLZvalid)';       %to find tares AAM042016
    
    aprxINminLZvalid = comINminLZvalid'*xvalid;
    
    for m=1:length(aprxINvalid(:,1))
        %%%%% 3/23/17 Zap intercepts %%%
        %        aprxINminGZvalid(m,:) = aprxINvalid(m,:)+interceptsvalid;
        aprxINminGZvalid(m,:) = aprxINvalid(m,:);
        %%%%%%
        
        checkitvalid(m,:) = aprxINminGZvalid(m,:)-targetMatrixvalid(m,:);
    end
    
    % SOLVE FOR TARES BY TAKING THE MEAN
    taresAllPointsvalid = meantare(seriesvalid,checkitvalid);
    zapvalid     = taresAllPointsvalid(s_1st,:);
    %RESIDUAL
    targetResvalid = targetMatrixvalid-aprxINminGZvalid+taresAllPointsvalid;
    
    %targetMatrixGlobalvalid = targetMatrixvalid + taresAllPointsvalid; % just for testing ajm 5_7_18
    
    resSquarevalid = dot(targetResvalid,targetResvalid)';
    
    aprxINminGZvalidprime = targetMatrixvalid+taresAllPointsvalid;
    
    %OUTPUTS FOR VALIDATION ALGEBRAIC SECTION
    for k=1:length(targetResvalid(1,:))
        [goopvalid(k),kstarvalid(k)] = max(abs(targetResvalid(:,k)));
        goopValvalid(k) = abs(targetResvalid(kstarvalid(k),k));
        xCentvalid(k) = excessVecvalidkeep(kstarvalid(k),k);
        maxTargetsvalid(k) = max(targetResvalid(:,k));
        minTargetsvalid(k) = min(targetResvalid(:,k));
    end
    perGoopvalid = 100*(goopvalid./loadCapacitiesvalid);
    davariancevalid = var(targetResvalid);
    geevalid = mean(targetResvalid);
    standardDevZ = std(targetResvalid);
    standardDevvalid = standardDevZ';
    stdDevPercentCapacityvalid = 100*(standardDevvalid'./loadCapacitiesvalid);
    ratioGoopvalid = goopvalid./standardDevvalid';
    ratioGoopvalid(isnan(ratioGoopvalid)) = realmin;
    
    theminmaxbandvalid = 100*(abs(maxTargetsvalid + minTargetsvalid)./loadCapacitiesvalid);
    
    %OUTPUT HISTOGRAM PLOTS
    if hist_FLAG == 1
        figure('Name','Validation - ALGB','NumberTitle','off')
        for k0=1:length(targetResvalid(1,:))
            subplot(2,3,k0)
            binWidth = 0.25;
            edges = [-4.125:binWidth:4.125];
            h = histogram(targetResvalid(:,k0)/standardDevvalid(k0,:),edges,'Normalization','probability');
            centers = edges(1:end-1)+.125;
            values = h.Values*100;
            bar(centers,values,'barwidth',1)
            ylabel('% Data Pts');
            xlim([-4 4]);
            ylim([0 50]);
            hold on
            plot(linspace(-4,4,100),binWidth*100*normpdf(linspace(-4,4,100),0,1),'r')
            hold off
            xlabel(['\Delta',loadlist{k0},'/\sigma']);
        end
    end
    
    if print_FLAG == 1
        %
        % Full Algebraic Model
        if model_FLAG == 1
            disp(' ');
            disp('%%%%%%%%%%%%%%%%%');
            disp(' ');
            disp('VALIDATION RESULTS: Full Algebraic Model');
        end
        % Truncated Algebraic Model
        if model_FLAG == 2
            disp(' ');
            disp('%%%%%%%%%%%%%%%%%');
            disp(' ');
            disp('VALIDATION RESULTS: Truncated Algebraic Model');
        end
        % Linear Algebraic Model
        if model_FLAG == 3
            disp(' ');
            disp('%%%%%%%%%%%%%%%%%');
            disp(' ');
            disp('VALIDATION RESULTS: Linear Algebraic Model');
        end
        
        disp('  ');
        disp('Validation data file read =');
        disp(out.savePathval);
        disp('  ');
        disp('Number of validation data points =');
        disp(numptsvalid);
        disp('  ');
        
        alg_Tares_valid = array2table(zapvalid,'VariableNames',loadlist(1:dimFlag))
        
        mean_alg_Resids_sqrd_valid = array2table(resSquarevalid'./numptsvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Pcnt_Capacity_Max_Mag_Load_Resids_valid = array2table(perGoopvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Std_Dev_pcnt_valid = array2table(stdDevPercentCapacityvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Max_Load_Resids_valid = array2table(maxTargetsvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Min_Load_Resids_valid = array2table(minTargetsvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Ratio_Max_Mag_Load_Resid_and_Std_Dev_valid = array2table(ratioGoopvalid,'VariableNames',loadlist(1:dimFlag))
        
        % Prints the minmaxband
        alg_per_minmaxband_valid = array2table(theminmaxbandvalid,'VariableNames',loadlist(1:dimFlag))
    end
    
    if excel_FLAG == 1
        %%%%
        disp('ALG VALIDATION MODEL GLOBAL LOAD APPROXIMATION: VALID_AOX_GLOBAL_ALG_RESULT in Workspace');
        disp(' ');
        
        filename = 'VALID_AOX_GLOBAL_ALG_RESULT.csv';
        csvwrite(filename,aprxINminGZvalid)
    end
    
    if res_FLAG == 1
        figure('Name','Algebraic Model Validation; Residuals of Load Versus Data Point Index','NumberTitle','off')
        plotResPages(seriesvalid, targetResvalid, loadCapacities, stdDevPercentCapacityvalid, loadlist)
        %    hold off
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    RBF SECTION FOR VALIDATION     AJM 12/10/16                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %goal to use centers, width and coefficients to validate parameters against
    %independent data
    
    for k=1:dimFlagvalid
        for s=1:length(seriesvalid)
            targetResvalidX(s,k) = targetResvalid(s,k);
        end
    end
    targetRes2valid = targetResvalidX;
    aprxINminGZ2valid = aprxINminGZvalid;
    
    % Subtract the Global Zeros from the Inputs and Local Zeros
    for k=1:dimFlagvalid
        dainputsvalid(:,k) = excessVecvalid0(:,k)-globalZerosvalid(k);
        dalzvalid(:,k) = localZerosAllPointsvalid(:,k)-globalZerosvalid(k);
        dagzvalid(:,k) = 0;
    end
    
    for k=1:dimFlag % ajm 6_8_18
        dainputscalib(:,k) = excessVec0(:,k) - globalZeros(k);
    end
    
    if balCal_FLAG == 2
        etaHistvalid = cell(numBasis,1);
        aprxINminGZ_Histvalid = cell(numBasis,1);
        tareHistvalid = cell(numBasis,1);
        
        [~,~,s_id] = unique(seriesvalid);
        localZeroMatrixvalid = localZerosvalid(s_id,:);
        
        etaLZvalid = dot(localZeroMatrixvalid-dainputsvalid,localZeroMatrixvalid-dainputsvalid);
        etaGZvalid = dot(globalZerosAllPointsvalid-dainputsvalid,globalZerosAllPointsvalid-dainputsvalid);
        
        for u=1:numBasis
            for s=1:length(excessVec0(1,:)) % loops through the components
                
                centerIndexLoop(s) = centerIndexHist(u,s); %Have to use the history or it gets overwritten
                
                for r=1:length(excessVecvalid(:,1))
                    adiffervalid(r,:) = dainputscalib(centerIndexHist(u,s),:)-dainputsvalid(r,:);
                    etavalid(r,s) = dot(adiffervalid(r,:),adiffervalid(r,:));
                end
                
                w(s) = wHist(u,s); % Have to use the history or it gets overwritten
                
                rbfINminLZvalid(:,s)=exp(etavalid(:,s)*log(abs(wHist(u,s)))) - exp(etaLZvalid(:,s)*log(abs(wHist(u,s))));
                rbfINminGZvalid(:,s)=exp(etavalid(:,s)*log(abs(wHist(u,s))));
                rbfLZminGZvalid(:,s)=exp(etaLZvalid(:,s)*log(abs(wHist(u,s))));%to find tares AAM042016
                coeffvalid(s) = cHist(u,s); %Have to use the history or it gets overwritten
                
                rbfc_INminLZvalid(:,s) = cHist(u,s)*rbfINminLZvalid(:,s);
                rbfc_INminGZvalid(:,s) = cHist(u,s)*rbfINminGZvalid(:,s);
                rbfc_LZminGZvalid(:,s) = cHist(u,s)*rbfLZminGZvalid(:,s);%to find tares AAM042016
            end
            
            %update the approximation
            aprxINminGZ2valid = aprxINminGZ2valid+rbfc_INminGZvalid;
            aprxINminGZ_Histvalid{u} = aprxINminGZ2valid;
            
            rbfc_INminGZ_Histvalid{u} = rbfc_INminGZvalid;  % temp ajm 6_7_18
            
            rbf_etavalid_Hist{u} = etavalid;   % temp ajm 6_7_18
            rbf_excessvec_center_Hist{u} = excessVec0(centerIndexHist(u,:),:);  % temp ajm 6_7_18
            rbf_excessvec_valid_Hist{u} = excessVecvalid;  % temp ajm 6_7_18
            
            % SOLVE FOR TARES BY TAKING THE MEAN
            [~,s_1st,~] = unique(seriesvalid);
            taresAllPointsvalid2 = meantare(seriesvalid,aprxINminGZ2valid-targetMatrixvalid)
            
            taresGRBFvalid = taresAllPointsvalid2(s_1st,:);
            tareHistvalid{u} = taresGRBFvalid;
            
            targetRes2valid = targetMatrixvalid+taresAllPointsvalid2-aprxINminGZ2valid;      %0=b-Ax
            targetMatrixGlobalGRBFvalid = targetMatrixvalid+taresAllPointsvalid2;  % temp for ajm 6_7_18
            
            newRes2valid = targetRes2valid'*targetRes2valid;
            resSquare2valid = diag(newRes2valid);
            resSquareHistvalid(u,:) = resSquare2valid;
        end
        
        for b=1:nseriesvalid
            for c=1:length(excessVecvalid(1,:))
                for a=1:numBasis
                    tempvalid(a) = tareHistvalid{a}(b,c);
                end
                tare3Sigvalid(b,c) = 3*std(tempvalid);
            end
        end
        
        total_rbfc_INminGZ_valid = rbfc_INminGZ_Histvalid{1};            % temp ajm 6_7_18
        
        for u=2:numBasis
            total_rbfc_INminGZ_valid = total_rbfc_INminGZ_valid + rbfc_INminGZ_Histvalid{u};            % temp ajm 6_7_18
        end
        
        %OUTPUTS FOR GRBF VALIDATION SECTION
        for k2=1:length(targetResvalid(1,:))
            [goop2valid(k2),kstar2valid(k2)] = max(abs(targetRes2valid(:,k2)));
            goop2valid(k2) = abs(targetRes2valid(kstar2valid(k2),k2));
            %        xCent2valid(k2) = excessVec0(centerIndexLoop(s)(k2),k2);%%??
            maxTargets2valid(k2) = max(targetRes2valid(:,k2));
            minTargets2valid(k2) = min(targetRes2valid(:,k2));
        end
        perGoop2valid = 100*(goop2valid./loadCapacitiesvalid);
        
        standardDev20valid = std(targetRes2valid);
        standardDev2valid = standardDev20valid';
        stdDevPercentCapacity2valid = 100*(standardDev2valid'./loadCapacitiesvalid);
        ratioGoop2valid = goop2valid./standardDev2valid';
        ratioGoop2valid(isnan(ratioGoop2valid)) = realmin;
        
        theminmaxband2valid = 100*(abs(maxTargets2valid + minTargets2valid)./loadCapacitiesvalid);
        
        %OUTPUT HISTOGRAM PLOTS
        if hist_FLAG == 1
            figure('Name','Validation - GRBF','NumberTitle','off')
            for k0=1:length(targetRes2valid(1,:))
                subplot(2,3,k0)
                binWidth = 0.25;
                edges = [-4.125:binWidth:4.125];
                h = histogram(targetRes2valid(:,k0)/standardDev2valid(k0,:),edges,'Normalization','probability');
                centers = edges(1:end-1)+.125;
                values = h.Values*100;
                bar(centers,values,'barwidth',1)
                ylabel('% Data Pts');
                xlim([-4 4]);
                ylim([0 50]);
                hold on
                plot(linspace(-4,4,100),binWidth*100*normpdf(linspace(-4,4,100),0,1),'r')
                hold off
                xlabel(['\Delta',loadlist{k0},'/\sigma']);
            end
        end
        if print_FLAG == 1
            %
            disp(' ***** ');
            disp(' ');
            disp('Number of GRBFs =');
            disp(numBasis);
            
            
            twoSigmaGRBFvalid = standardDevvalid'.*2;
            GRBF_2Sigmavalid = array2table(twoSigmaGRBFvalid,'VariableNames',loadlist(1:dimFlag))
            
            %Should I use strtrim()  ? -AAM 042116
            GRBF_Taresvalid = array2table(taresGRBFvalid,'VariableNames',loadlist(1:dimFlag))
            
            mean_GRBF_Resids_sqrdvalid = array2table(resSquare2valid'./numptsvalid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Pcnt_Capacity_Max_Mag_Load_Resid_valid = array2table(perGoop2valid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Std_Dev_pcnt_valid = array2table(stdDevPercentCapacity2valid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Max_Load_Resids_valid = array2table(maxTargets2valid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Min_Load_Resids_valid = array2table(minTargets2valid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Ratio_Max_Mag_Load_Resid_and_Std_Dev_valid = array2table(ratioGoop2valid,'VariableNames',loadlist(1:dimFlag))
            
            % Prints the GRBF minmax
            GRBF_minmaxband_per_capacity_valid = array2table(theminmaxband2valid,'VariableNames',loadlist(1:dimFlag))
        end
        
        if res_FLAG == 1
            figure('Name','GRBF + Algebraic Model Validation; Residuals of Load Versus Data Point Index','NumberTitle','off')
            plotResPages(seriesvalid, targetRes2valid, loadCapacities, stdDevPercentCapacity2valid, loadlist)
            %    hold off
        end
        
        % Diagnostics %%%%% 5/16/18
        %justtherbfsvalid = aprxINminGZ2valid-aprxINminGZvalid;
        %
        %%deltatherbfsvalid = justtherbfs - justtherbfsvalid;
        %
        %if res_FLAG == 1
        %    figure('Name','Looking at GRBF Distribution in Validation','NumberTitle','off')
        %    plotResPages(seriesvalid, justtherbfsvalid, loadCapacities, stdDevPercentCapacity2valid)
        %    hold off
        %end
        
        if excel_FLAG == 1
            disp(' ');
            disp('ALG+GRBF VALIDATION MODEL GLOBAL LOAD APPROXIMATION: Check VALID_AOX_GLOBAL_GRBF_RESULT.csv file');
            disp(' ');
            
            filename = 'VALID_AOX_GLOBAL_GRBF_RESULT.csv';
            csvwrite(filename,aprxINminGZ2valid)
        end
    end
end

if balApprox_FLAG == 1
    %
    %
    % Copyright ©2016 Andrew Meade and Ali Arya Mokhtarzadeh.  All Rights Reserved.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                        APPROXIMATION SECTION      AJM 6/29/17           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %DEFINE THE PRODUCTION CSV INPUT FILE AND SELECT THE RANGE OF DATA VALUES TO READ
    %
    
    load(out.savePathapp,'-mat');
    %
    
    % num of data points
    numptsapprox = length(excessVecapprox);
    %
    
    
    %natural zeros (also called global zeros)
    globalZerosapprox = mean(natzerosapprox);
    
    %%% make an array out of the globalZerosapprox vector
    for i=1:numptsapprox
        globalZerosAllPointsapprox(i,:) = globalZerosapprox;
    end
    %%%
    
    
    
    %% Subtract the Global Zeros from the Inputs %%%%%%%%%%
    
    for k=1:dimFlag
        
        dainputsapprox(:,k) = excessVecapprox(:,k)-globalZerosAllPointsapprox(:,k);
        
        dalzapprox(:,k) = globalZerosAllPointsapprox(:,k)-globalZerosAllPointsapprox(:,k);
        
    end
    %%%%%%%%%%%%
    
    
    %%
    %% Build the Algebraic Model
    %%
    
    %     n(1) = 2*dimFlag*(dimFlag+2);
    %     n(2) = dimFlag*(dimFlag+3)/2;
    %     n(3) = dimFlag;
    %     model_FLAG = find(n==size( xapproxer,1)-1);
    
    
    %% Full Algebraic Model
    if model_FLAG == 1
        nterms = 2*dimFlag*(dimFlag+2);
    end
    
    %% Truncated Algebraic Model
    if model_FLAG == 2
        nterms = dimFlag*(dimFlag+3)/2;
    end
    
    %% Linear Algebraic Model
    if model_FLAG == 3
        nterms = dimFlag;
    end
    
    
    
    % Call the Algebraic Subroutine
    %
    
    comGZapprox= zeros(nterms,1);
    
    
    for i=1:dimFlag
        biggee(:,i) = 0;
    end
    
    [comINapprox,comLZapprox,comGZapprox]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numptsapprox,0,0,dainputsapprox,dalzapprox,biggee);
    
    %model_FLAG
    %nterms
    %dimFlag
    %numptsapprox
    
    
    %%
    %%
    
    
    %LOAD APPROXIMATION
    %define the approximation for inputs minus global zeros
    interceptsapprox = -(comGZapprox'*xapprox);
    aprxINapprox = ( xapprox'*comINapprox)';        %to find ?? AJM111516
    %%
    %%
    for m=1:length(aprxINapprox)
        aprxINminGZapprox(m,:) = aprxINapprox(m,:);
    end
    %%
    %%
    
    %%%%%%
    if excel_FLAG == 1
        disp(' ');
        disp('%%%%%%%%%%%%%%%%%');
        disp(' ');
        disp('ALG MODEL APPROXIMATION RESULTS: Check Global_ALG_Approx.csv in file');
        
        filename = 'Global_ALG_Approx.csv';
        csvwrite(filename,aprxINminGZapprox)
        
    else
        
        disp(' ');
        disp('%%%%%%%%%%%%%%%%%');
        disp(' ');
        disp('ALG MODEL APPROXIMATION RESULTS: Check aprxINminGZapprox in Workspace');
    end
    %%%%%%
    
    
    
    
    %
    
    
    %
    %
    %
    % Copyright ©2016 Andrew Meade and Ali Arya Mokhtarzadeh.  All Rights Reserved.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    RBF SECTION FOR APPROXIMATION     AJM 6/29/17                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %goal to use centers, width and coefficients to approxate parameters against
    %independent data
    
    %%
    
    aprxINminGZ2approx = aprxINminGZapprox;
    
    
    %%
    
    if balCal_FLAG == 2
        
        etaHistapprox = cell(numBasis,1);
        aprxINminGZ_Histapprox = cell(numBasis,1);
        
        
        etaGZapprox = dot(globalZerosAllPointsapprox-excessVecapprox,globalZerosAllPointsapprox-excessVecapprox);
        
        for u=1:numBasis
            for s=1:length(excessVec0(1,:)) % loops through the 8 components
                
                centerIndexLoop(s) = centerIndexHist(u,s); %Have to use the history or it gets overwritten
                
                for r=1:length(excessVecapprox(:,1))
                    etaapprox(r,s) = dot(excessVecapprox(r,:)-excessVec0(centerIndexLoop(s),:),excessVecapprox(r,:)-excessVec0(centerIndexLoop(s),:));
                end
                
                w(s) = wHist(u,s); % Have to use the history or it gets overwritten
                
                
                rbfINminGZapprox(:,s)=exp(etaapprox(:,s)*log(abs(w(s))));
                coeffapprox(s) = cHist(u,s); %Have to use the history or it gets overwritten
                
                rbfc_INminGZapprox(:,s) = coeffapprox(s)*rbfINminGZapprox(:,s);
                
                
            end
            
            
            
            wHistapprox(u,:) = w;
            cHistapprox(u,:) = coeffapprox;
            centerIndexHist(u,:) = centerIndexLoop;
            etaHistapprox{u} = etaapprox;
            
            %UPDATE THE RESIDUAL
            
            %update the approximation
            
            aprxINminGZ2approx = aprxINminGZ2approx+rbfc_INminGZapprox;
            aprxINminGZ_Histapprox{u} = aprxINminGZ2approx;
            
        end
        
        
        
        if excel_FLAG == 1
            disp(' ');
            disp('%%%%%%%%%%%%%%%%%');
            disp(' ');
            disp('ALG + GRBF MODEL APPROXIMATION RESULTS: Check Global_ALG+GRBF_Approx.csv in file');
            
            filename = 'Global_ALG+GRBF_Approx.csv';
            csvwrite(filename,aprxINminGZ2approx)
        else
            disp(' ');
            disp('GRBF MODEL APPROXIMATION RESULTS: Check aprxINminGZ2approx in Workspace');
            disp(' ');
        end
        
    end
    
    %End Approximation Option
end

disp('  ')
disp('Calculations Complete.')

%
% Copyright ©2016 Andrew Meade and Ali Arya Mokhtarzadeh.  All Rights Reserved.
%

%toc
