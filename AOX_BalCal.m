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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
load(out.savePathcal,'-mat');
series0 = series;
%
%                       END USER INPUT SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6_14_18 ajm
warning('off','all');
%

disp('Starting Calculations')
if approach_FLAG == 1
    disp(' ');
    disp('Using the Indirect Approach for Calibration');
else
    disp(' ');
    disp('Using the Direct Approach for Calibration');
end

% 6_14_18 ajm
matrixcolumnlabels = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};

loadlist = {'NF','BM','S1','S2','RM','AF','PLM', 'PCM', 'MLM', 'MCM'};
voltagelist = {'rNF','rBM','rS1','rS2','rRM','rAF','rPLM','rPCM','rMLM','rMCM'};

if corr_FLAG == 1
    figure('Name','Correlation plot','NumberTitle','off');
    correlationPlot(targetMatrix0, excessVec0, loadlist, voltagelist);
end

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
        [series,targetMatrix,excessVec,sample] = AOX_LHS(series0,targetMatrix0,excessVec0,LHSp);
        lhs_check(sample) = 1;
        ind(find(lhs_check-1)); % This line outputs which data points haven't been sampled yet
        pct_sampled = sum(lhs_check)/length(lhs_check); % This line outputs what percentage of points have been sampled
    else
        excessVec = excessVec0;
        targetMatrix = targetMatrix0;
        series = series0;
    end

    [numpts, dimFlag] = size(excessVec);

    %find the average natural zeros (also called global zeros)
    globalZeros = mean(natzeros);

    [~,s_1st,s_id] = unique(series);
    %find number of series; this will tell us the number of tares
    nseries = length(s_1st);

    %find zero points of each series and number of points in a series
    %localZerosAllPoints is the same as localZeroMatrix defined in the RBF
    %section
    [localZeros,localZerosAllPoints] = localzeros(series,excessVec);
    globalZerosAllPoints = ones(numpts,1)*globalZeros;

    disp('  ')
    disp('Working ...')
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 VOLTAGE TO LOAD (DIRECT) SECTION      AJM 8/3/17           %
    %Use the measured to best voltages mapping to find the calibration coefficients
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %

    % Indirect approach uses the modeled voltage
    if approach_FLAG == 1
        if balCal_FLAG == 2
            excessVec = qtaprxINminGZ2 + globalZerosAllPoints;
        else
            excessVec = qtaprxINminGZ + globalZerosAllPoints;
        end
    end

    %%% Subtract the Global Zeros from the Inputs and Local Zeros %%%%%%%%%%
    dainputs = excessVec - globalZerosAllPoints;
    dalz = localZerosAllPoints - globalZerosAllPoints;
    for i=1:dimFlag
        dainputs(:,i) = excessVec(:,i)-globalZeros(i);
        %
        dalz(:,i) = localZerosAllPoints(:,i)-globalZeros(i);
        %
        biggee(:,i) = 0;
        %
    end

    %%% Build the Algebraic Model

    lasttare = series(numpts);

    % Full Algebraic Model
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
    comIN = balCal_algEqns(model_FLAG,dainputs);
    [comIN,comLZ,comGZ,uncert_comIN]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numpts,series,lasttare,dainputs,dalz,biggee);


    % Effectively removes the original tare values so we can calculate the averages
    for i=1:lasttare
        comIN(nterms+i,:) = 0;
        comLZ(nterms+i,:) = 0;
        comGZ(nterms+i,:) = 0;
    end

    %%
    %       %CHANGED JRP
    %     for loopk=1:numpts
    %         comLZ(nterms+series(loopk),loopk) = 1.0;
    %     end
    %
    %
    %%


    comINminLZ = comIN-comLZ;

    %START: Changes made for bootstrap CI for coefficients, tare from mean
    %(zap)
    %START issues with resampling
    %Moved outside of loop
    prevVal = 1;
    for i=1:length(series)
        if series(i)==prevVal
            localZeros(prevVal,:) = excessVec0(i,:);
            localZerosAllPoints(i,:) = localZeros(prevVal,:);
            prevVal = prevVal+1;
        else
            counter(prevVal) = i;
            localZerosAllPoints(i,:) = localZeros(prevVal-1,:);
        end
        %        globalZerosAllPoints(i,:) = globalZeros;
    end

    %add '1' to every value of counter to find the indices of the local zeros
    indexLocalZero = counter(:)+1;

    %%start function
    bootalpha=.05;
    nbootstrap=200;
    f=@zapFinder;
    [fout]=f(comINminLZ',targetMatrix,series,excessVec0,targetMatrix0,globalZerosAllPoints,localZerosAllPoints,dimFlag,model_FLAG,nterms,numpts,lasttare,nseries);
    fzap=fout(1:nseries,:);
    faprxLZminGZ_series=fout(size(fzap,1)+1:size(fzap,1)+nseries,:);
    fxcalib=fout(size(fzap,1)+size(faprxLZminGZ_series)+1:size(fout,1),:);

    f_ci=bootci(nbootstrap,{f,comINminLZ',targetMatrix,series,excessVec0,targetMatrix0,globalZerosAllPoints,localZerosAllPoints,dimFlag,model_FLAG,nterms,numpts,lasttare,nseries}, 'type', 'cper','alpha',bootalpha);
    fzap_ci=f_ci(:,1:nseries,:);
    faprxLZminGZ_ci=f_ci(:,size(fzap,1)+1:size(fzap,1)+nseries,:);
    fxcalib_ci=f_ci(:,size(fzap,1)+size(faprxLZminGZ_series)+1:size(fout,1),:);

    for loopk=1:numpts
        comLZ(nterms+series(loopk),loopk) = 1.0;
    end

    comINminLZ = comIN-comLZ;


    % END: bootstrap section

    %%
    %SOLUTION
    xcalib = comINminLZ'\targetMatrix;                    % '\' solution of Ax=b


    %    xcalib = lsqminnorm(comINminLZ',targetMatrix);             % alternate solution method
    %    xcalib = pinv(comINminLZ')*targetMatrix;             % alternate solution method
    %    xcalib = pinv(comINminLZ',1e-2)*targetMatrix;             % alternate solution method
    % 5/17/18
    for i=1:nterms+1
        xvalid(i,:) = xcalib(i,:);
    end

    %  Creates Matrix for the volts to loads
    %  ajm for the users to view 5/11/18
    APPROX_AOX_COEFF_MATRIX = xcalib;

    for i=1:nterms
        xapproxer(i,:) = xcalib(i,:);
    end

    xapproxer(nterms+1,:) = 0.0;

    if excel_FLAG == 1
        filename = 'APPROX_AOX_COEFF_MATRIX.csv';
        Z = xapproxer;
        %     xlRange = 'matrixcolumnlabels(1)1:matrixcolumnlabels(dimFlag)nterms';
        xlRange=char(strcat(matrixcolumnlabels(1),'1:',matrixcolumnlabels(dimFlag),num2str(nterms)));
        xlswrite(filename,Z,xlRange)
        %   'made it here'
        %%
    end

    if LHS_Flag == 1
        x_all(:,:,lhs) = xcalib;
    end
end

if LHS_Flag == 1
    xcalib = mean(x_all,3);
    xcalib_std = std(x_all,[],3);
end

excessVec = excessVec;
targetMatrix = targetMatrix0;
series = series0;

[localZeros,localZerosAllPoints] = localzeros(series,excessVec);

for i=1:dimFlag
    dainputs2(:,i) = excessVec(:,i)-globalZeros(i);
    dalz2(:,i) = localZerosAllPoints(:,i)-globalZeros(i);
    biggee(:,i) = 0;
end

% Call the Algebraic Subroutine
[comIN,comLZ,comGZ,uncert_comIN]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numpts,series,lasttare,dainputs2,dalz2,biggee);


% Effectively removes the original tare values so we can calculate the averages
for i=1:lasttare
    comIN(nterms+i,:) = 0;
    comLZ(nterms+i,:) = 0;
    comGZ(nterms+i,:) = 0;
end


%CHANGED 9 JAN 19 JRP
% for loopk=1:numpts
%     comLZ(nterms+series(loopk),loopk) = 1.0;
% end


%%
%%
%ADDED 9 Jan 19 JRP: Monte carlo
carlo=0;
if carlo==1
    nCarlo=10000;
    for i=1:nseries
        localZeros_da(i,:)=localZeros(i,:)-globalZeros(:)';
    end

    % %Define magnitude of noise:
    % percentVoltage=0.01; %percent of max recorded voltage that is noise
    % Maxnoise=percentVoltage*max(abs(dainputs2));
    Maxnoise=1;

    %equally distributed noise
    % sigma=Maxnoise./2; % 95% of the noise points will be within the defined Maxnoise

    sigma=1/2; %Trust all channels down to 1 microvolt, sigma is trust/2

    %Monte carlo for uncert in tares due to uncert in read voltages
    for i=1:nseries
        % Normally distributed noise
        noise=sigma.*randn(nCarlo,size(localZeros_da,2));

        % % Equally distributed Noise
        % noise=-Maxnoise+(2*Maxnoise).*rand(nCarlo,size(localZeros_da,2));

        series_tare=zeros(nCarlo,1);
        series_tare(:,1)=i;
        dalz2_tare=zeros(nCarlo,size(localZeros_da,2));

        for j=1:nCarlo
            dalz2_tare(j,:)=localZeros_da(i,:);
        end
        dainputs2_tare=dalz2_tare+noise;

        [comIN_tare,comLZ_tare,comGZ_tare,uncert_comIN_tare]=balCal_algEquations3(model_FLAG,nterms,dimFlag,nCarlo,series_tare,i,dainputs2_tare,dalz2_tare,biggee);
        for j=1:lasttare
            comIN_tare(nterms+j,:) = 0;
            comLZ_tare(nterms+j,:) = 0;
            comGZ_tare(nterms+j,:) = 0;
        end
        aprxIN_tare = (xcalib'*comIN_tare)';
        aprxLZ_tare = (xcalib'*comLZ_tare)';       %to find tares AAM042016
        tareStd(i,:)=std(aprxIN_tare);
        tare(i,:)=mean(aprxLZ_tare);
    end

    %Monte carlo for uncert in all datapoints due to uncert in read voltages
    for i=1:size(dainputs2,1) %Change 20 to # datapoints size(dainputs2,1)

        % Normally distributed noise
        noise=sigma.*randn(nCarlo,size(localZeros_da,2));

        % Equally distributed Noise
        % noise=-Maxnoise+(2*Maxnoise).*rand(nCarlo,size(localZeros_da,2));

        for j=1:nCarlo
            dalz2_carlo(j,:)=dalz2(i,:);
            dainputs2_carlo(j,:)=dainputs2(i,:);
            series_carlo(j)=series(i);
        end
        dainputs2_carlo=dainputs2_carlo+noise;
        [comIN_carlo,comLZ_carlo,comGZ_carlo,uncert_comIN_carlo]=balCal_algEquations3(model_FLAG,nterms,dimFlag,nCarlo,series_carlo,max(series_carlo),dainputs2_carlo,dalz2_carlo,biggee);
        for j=1:lasttare
            comIN_carlo(nterms+j,:) = 0;
            comLZ_carlo(nterms+j,:) = 0;
            comGZ_carlo(nterms+j,:) = 0;
        end
        aprxIN_carlo = (xcalib'*comIN_carlo)';
        aprxLZ_carlo = (xcalib'*comLZ_carlo)';       %to find tares AAM042016
        load_std_carlo(i,:)=std(aprxIN_carlo);
    end
end
%END added for monte-carlo


comINminLZ = comIN-comLZ;

%APPROXIMATION
%define the approximation for inputs minus local zeros
aprxINminLZ = comINminLZ'*xcalib;

%DIFFERENT APPROXIMATION
%define the approximation for inputs minus global zeros
intercepts = -(comGZ'*xcalib);
aprxIN = (xcalib'*comIN)';
aprxLZ = (xcalib'*comLZ)';       %to find tares AAM042016
for m=1:length(aprxIN)
    %    aprxINminGZ(m,:) = aprxIN(m,:)+intercepts;  %ajm zap 3/23/18
    %    aprxLZminGZ(m,:) = aprxLZ(m,:)+intercepts;  %ajm zap 3/23/18
    aprxINminGZ(m,:) = aprxIN(m,:);
    aprxLZminGZ(m,:) = aprxLZ(m,:); %to find tares AAM042016

    %subtracts the targets from the global approx
    %This will be averaged over the individual series to give the tares
    checkit(m,:) = aprxINminGZ(m,:)-targetMatrix(m,:);
end

taretal = meantare(series,checkit);
%RESIDUAL
targetRes = targetMatrix+taretal-aprxINminGZ;      %0=b-Ax

reslist = {'resNF','resBM','resS1','resS2','resRM','resAF','resPLM',...
    'resPCM', 'resMLM', 'resMCM'};
if rescorr_FLAG == 1
    figure('Name','Residual correlation plot','NumberTitle','off');
    correlationPlot(excessVec, targetRes, voltagelist, reslist);
end

%find the sum of squares of the residual using the dot product
resSquare = dot(targetRes,targetRes)';
%AAM note to self - in matlab, diag(A'*A) is the same as dot(A,A)'

%Run function to calculate uncertainty on loads from calibration
[combined_uncert,tare_uncert, FL_uncert]=uncert_prop(xcalib,fxcalib_ci,comIN,dimFlag,uncert_comIN,indexLocalZero,lasttare,nterms,aprxIN,series);

%---------------------------------------------------------------
%---------------------------------------------------------------
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
        zeroed_targetMatrix = targetMatrix;
        zeroed_excessVec = excessVec;
        zeroed_series = series;
        zeroed_numpts = numpts;

        zeroed_targetMatrix(OUTLIER_ROWS,:) = [];
        zeroed_excessVec(OUTLIER_ROWS,:) = [];

        targetMatrix = unique(zeroed_targetMatrix,'rows');
        excessVec = unique(zeroed_excessVec,'rows');

        zeroed_series(OUTLIER_ROWS) = [];

        numpts =  zeroed_numpts - num_outliers;

        targetMatrix = zeroed_targetMatrix;
        excessVec = zeroed_excessVec;
        series = zeroed_series;

        dainputs = zeros(numpts,dimFlag);
        dalz = zeros(numpts,dimFlag);
        localZerosAllPoints = zeros(numpts,dimFlag);
        aprxINminGZ = zeros(numpts,dimFlag);
        aprxLZminGZ = zeros(numpts,dimFlag);
        zeroed_checkit = zeros(numpts,dimFlag);
        taretal = zeros(numpts,dimFlag);
        globalZerosAllPoints = zeros(numpts,dimFlag);
        eta = zeros(numpts,dimFlag);
        zeroed_zoop = zeros(numpts,dimFlag);

        rbfc_INminLZ = zeros(numpts,dimFlag);
        rbfc_INminGZ = zeros(numpts,dimFlag);
        rbfc_LZminGZ = zeros(numpts,dimFlag);
        rbfINminLZ = zeros(numpts,dimFlag);
        rbfINminGZ = zeros(numpts,dimFlag);
        rbfLZminGZ = zeros(numpts,dimFlag);

        [localZeros,localZerosAllPoints] = localzeros(series,excessVec);
        globalZerosAllPoints = ones(numpts,1)*globalZeros;

        % Subtract the Global Zeros from the Inputs and Local Zeros
        for i=1:dimFlag
            dainputs(:,i) = excessVec(:,i) - globalZerosAllPoints(i);
            dalz(:,i) = localZerosAllPoints(:,i) - globalZerosAllPoints(i);
        end

        for i=1:dimFlag
            biggee(:,i) = 0;
        end
        [comIN,comLZ]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numpts,series,lasttare, dainputs, dalz, biggee);

        comINminLZ = comIN-comLZ;

        %SOLUTION
        xcalib = pinv(comINminLZ')*targetMatrix;

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
            %This will be averaged over the individual series to give the tares
            checkit(m,:) = aprxINminGZ(m,:)-targetMatrix(m,:);
        end

        taretal = meantare(series,checkit);
        %RESIDUAL
        targetRes = targetMatrix+taretal-aprxINminGZ;      %0=b-Ax
    end
end

%find the sum of squares of the residual using the dot product
resSquare = dot(targetRes,targetRes)';
%AAM note to self - in matlab, diag(A'*A) is the same as dot(A,A)'

%OUTPUTS FOR ALGEBRAIC SECTION
for k=1:length(targetRes(1,:))
    [goop(k),kstar(k)] = max(abs(targetRes(:,k)));
    goopVal(k) = abs(targetRes(kstar(k),k));
    xCent(k) = excessVec(kstar(k),k);
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

[~,s_1st,~] = unique(series);
taresALGB = taretal(s_1st,:);

%OUTPUT HISTOGRAM PLOTS
if hist_FLAG == 1
    for k0=1:length(targetRes(1,:))
        figure;
        [histALGB, binValues] = hist(targetRes(:,k0)/standardDev(k0,:),20);
        normalizedCounts = 100 * histALGB / sum(histALGB);
        bar(binValues, normalizedCounts, 'barwidth', 1);
        xlabel('Ratio Between Load Residual and Standard Deviation');
        ylabel('Number of Readings in % of Number of Data Points');
        xlim([-4 4]);
        ylim([0 50]);
        title(strrep(['Histogram of ALG Calibration Model for %s',7],'%s',loadlist(k0)));
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
    numCombinName(nterms+lasttare+1) = cellstr('Intercept');

    calib_mean_algebraic_Resids_sqrd = array2table(resSquare'./numpts,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Pcnt_Capacity_Max_Mag_Load_Resids = array2table(perGoop,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Std_Dev_pcnt = array2table(stdDevPercentCapacity,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Max_Load_Resids = array2table(maxTargets,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Min_Load_Resids = array2table(minTargets,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Ratio_Max_Mag_Load_Resid_and_Std_Dev = array2table(ratioGoop,'VariableNames',loadlist(1:dimFlag))

    % Prints the minmaxband
    calib_alg_per_minmaxband = array2table(theminmaxband,'VariableNames',loadlist(1:dimFlag))
    %%%%%%%%%

    if excel_FLAG == 1
        % Output results to an excel file
        disp('  ');
        disp('ALG CALIBRATION MODEL GLOBAL LOAD APPROXIMATION FILE: CALIB_AOX_GLOBAL_ALG_RESULT.csv');
        % CALIB_AOX_GLOBAL_ALG_RESULT = aprxINminGZ;
        disp(' ');

        filename = 'CALIB_AOX_GLOBAL_ALG_RESULT.csv';
        Z = aprxINminGZ;
        %        xlRange = 'A1:Jnumpts';
        xlRange=char(strcat('A1:J',num2str(numpts)));
        xlswrite(filename,Z,xlRange)
    end

end

if res_FLAG == 1
    figure('Name','Algebraic Model Calibration; Residuals of Load Versus Data Point Index','NumberTitle','off')
    plotResPages(series, targetRes, loadCapacities, stdDevPercentCapacity )
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
        for s=1:length(series)
            targetRes2(s,i) = targetRes(s,i);
        end
    end
    aprxINminGZ2 = aprxINminGZ;

    etaHist = cell(numBasis,1);
    aprxINminGZ_Hist = cell(numBasis,1);
    tareHist = cell(numBasis,1);

    for i=1:dimFlag
        dainputscalib(:,i) = excessVec(:,i)-globalZeros(i);
        dalzcalib(:,i) = localZerosAllPoints(:,i)-globalZeros(i);
    end

    %    localZeroMatrix = localZerosAllPoints;
    globalZerosAllPoints = zeros(length(excessVec(:,1)),dimFlag); % ajm 6_2_18

    etaLZ = dot(dalzcalib-dainputscalib,dalzcalib-dainputscalib);
    etaGZ = dot(globalZerosAllPoints-dainputscalib,globalZerosAllPoints-dainputscalib);

    for u=1:numBasis
        for s=1:dimFlag
            [goopLoop(s),centerIndexLoop(s)] = max(abs(targetRes2(:,s)));

            for r=1:length(excessVec(:,1))
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
        taretalGRBF = meantare(series,aprxINminGZ2-targetMatrix);

        taresGRBF = taretalGRBF(s_1st,:);

        tareGRBFHist{u} = taresGRBF;

        targetRes2 = targetMatrix-aprxINminGZ2+taretalGRBF;      %0=b-Ax
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
        plotResPages(series, targetRes2, loadCapacities, stdDevPercentCapacity2)
        %    hold off
    end

    %%% Diagnostics %%%%% 5/16/18
    %
    %justtherbfs = aprxINminGZ2-aprxINminGZ;
    %
    %if res_FLAG == 1
    %    figure('Name','Looking at GRBF Distribution in Calibration','NumberTitle','off')
    %    plotResPages(series, justtherbfs, loadCapacities, stdDevPercentCapacity2)
    %    hold off
    %end

    %OUTPUT HISTOGRAM PLOTS
    if hist_FLAG == 1 && balCal_FLAG == 2
        for k3=1:length(targetRes2(1,:))
            figure;
            [histGRBF2, binValues2] = hist(targetRes2(:,k3)/standardDev2(k3,:),20);
            normalizedCounts2 = 100 * histGRBF2 / sum(histGRBF2);
            bar(binValues2, normalizedCounts2, 'barwidth', 1);
            xlabel('Ratio Between Load Residual and Standard Deviation')
            ylabel('Number of Readings in % of Number of Data Points')
            xlim([-4 4]);
            ylim([0 50]);
            title(strrep(['Histogram of ALG+GRBF Calibration Model for %s',7],'%s',loadlist(k3)));
        end
    end

    if excel_FLAG == 1 && balCal_FLAG == 2
        disp(' ***** ');
        disp('  ');
        disp('ALG+GRBF CALIBRATION MODEL GLOBAL LOAD APPROXIMATION: Check CALIB_AOX_GLOBAL_GRBF_RESULT.csv file');

        filename = 'CALIB_AOX_GLOBAL_GRBF_RESULT.csv';
        Z = aprxINminGZ2;
        %     xlRange = 'A1:JnumBasis';
        xlRange=char(strcat('A1:J',num2str(numBasis)));
        xlswrite(filename,Z,xlRange)
        %%%%%%
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

    if excel_FLAG == 1 && balCal_FLAG == 2

        filename = 'APPROX_AOX_GRBF_ws.csv';
        Z = wHist;
        %     xlRange = 'A1:JnumBasis';
        xlRange=char(strcat('A1:J',num2str(numBasis)));
        xlswrite(filename,Z,xlRange)

        filename = 'APPROX_AOX_GRBF_coeffs.csv';
        Z = cHist;
        %     xlRange = 'A1:JnumBasis';
        xlRange=char(strcat('A1:J',num2str(numBasis)));
        xlswrite(filename,Z,xlRange)

        filename = 'APPROX_AOX_GRBF_Centers.csv';
        Z = centerIndexHist;
        %     xlRange = 'A1:JnumBasis';
        xlRange=char(strcat('A1:J',num2str(numBasis)));
        xlswrite(filename,Z,xlRange)
        %%%%%%%
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

    excessVecvalid0 = excessVecvalid;
    % num of data points
    numptsvalid = length(seriesvalid);
    dimFlagvalid = length(excessVecvalid(1,:));

    %find the average natural zeros (also called global zeros)
    globalZerosvalid = mean(natzerosvalid);

    %load capacities
    loadCapacitiesvalid(loadCapacitiesvalid == 0) = realmin;

    %find number of series; this will tell us the number of tares
    nseriesvalid = max(seriesvalid);

    %find zero points of each series and number of points in a series
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
    %Remember that  excessVec = excessVec0_complete - globalZerosAllPoints;
    excessVecvalidkeep = excessVecvalid0  - globalZerosAllPointsvalid;
    %%%

    % Build the Algebraic Model
    lasttarevalid = seriesvalid(numptsvalid);
    % Full Algebraic Model
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

    [sharedvals,indexLocalZerovalid]=intersect(seriesvalid, [1:max(seriesvalid)],'stable'); %Create index
    % Call the Algebraic Subroutine
    comGZvalid = zeros(nterms+1,1);

    [comINvalid,comLZvalid,comGZvalid,uncert_comINvalid]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numptsvalid,seriesvalid,1, dainputsvalid,dalzvalid,dagzvalid);
    %%


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
    taretalvalid = meantare(seriesvalid,checkitvalid);
    %RESIDUAL
    targetResvalid = targetMatrixvalid-aprxINminGZvalid+taretalvalid;

    %targetMatrixGlobalvalid = targetMatrixvalid + taretalvalid; % just for testing ajm 5_7_18

    resSquarevalid = dot(targetResvalid,targetResvalid)';

    aprxINminGZvalidprime = targetMatrixvalid+taretalvalid;

    %OUTPUTS FOR VALIDATION ALGEBRAIC SECTION

        %Run function to calculate uncertainty on loads output in approximation
[combined_uncertvalid,tare_uncertvalid, FL_uncertvalid]=uncert_prop(xvalid,fxcalib_ci,comINvalid,dimFlag,uncert_comINvalid,indexLocalZerovalid,lasttarevalid,nterms,aprxINvalid,seriesvalid);


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
    if hist_FLAG == 1 && balCal_FLAG == 2
        for k0=1:length(targetResvalid(1,:))
            figure;
            [histGRBFvalid, binValuesvalid] = hist(targetResvalid(:,k0)/standardDevvalid(k0,:),20);
            normalizedCountsvalid = 100 * histGRBFvalid / sum(histGRBFvalid);
            bar(binValuesvalid, normalizedCountsvalid, 'barwidth', 1);
            xlabel('Ratio Between Load Residual and Standard Deviation')
            ylabel('Number of Readings in % of Number of Data Points')
            xlim([-4 4]);
            ylim([0 50]);
            title(strrep(['Histogram of ALG Validation Model for %s',7],'%s',loadlist(k0)));
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

        if excel_FLAG == 1
            %%%%
            disp('ALG VALIDATION MODEL GLOBAL LOAD APPROXIMATION: VALID_AOX_GLOBAL_ALG_RESULT in Workspace');
            disp(' ');

            filename = 'VALID_AOX_GLOBAL_ALG_RESULT.csv';
            Z = aprxINminGZvalid;
            xlRange = 'matrixcolumnlabels(1)1:matrixcolumnlabels(dimFlag)numpts';
            xlswrite(filename,Z,xlRange)
            %%%%
        end
        zapvalid=taretalvalid(indexLocalZerovalid(1:(numel(indexLocalZerovalid)-1)),:); %tare loads
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

    if res_FLAG == 1
        figure('Name','Algebraic Model Validation; Residuals of Load Versus Data Point Index','NumberTitle','off')
        plotResPages(seriesvalid, targetResvalid, loadCapacities, stdDevPercentCapacityvalid )
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
        dainputscalib(:,k) = excessVec(:,k) - globalZeros(k);
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
            for s=1:length(excessVec(1,:)) % loops through the components

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
            rbf_excessvec_center_Hist{u} = excessVec(centerIndexHist(u,:),:);  % temp ajm 6_7_18
            rbf_excessvec_valid_Hist{u} = excessVecvalid;  % temp ajm 6_7_18

            % SOLVE FOR TARES BY TAKING THE MEAN
            [~,s_1st,~] = unique(seriesvalid);
            taretalvalid2 = meantare(seriesvalid,aprxINminGZ2valid-targetMatrixvalid)

            taresGRBFvalid = taretalvalid2(s_1st,:);
            tareHistvalid{u} = taresGRBFvalid;

            targetRes2valid = targetMatrixvalid+taretalvalid2-aprxINminGZ2valid;      %0=b-Ax
            targetMatrixGlobalGRBFvalid = targetMatrixvalid+taretalvalid2;  % temp for ajm 6_7_18

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
            %        xCent2valid(k2) = excessVec(centerIndexLoop(s)(k2),k2);%%??
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
        if hist_FLAG == 1 && balCal_FLAG == 2
            for k3=1:length(targetRes2valid(1,:))
                figure;
                [histGRBF2valid, binValues2valid] = hist(targetRes2valid(:,k3)/standardDev2valid(k3,:),20);
                normalizedCounts2valid = 100 * histGRBF2valid / sum(histGRBF2valid);
                bar(binValues2valid, normalizedCounts2valid, 'barwidth', 1);
                xlabel('Ratio Between Load Residual and Standard Deviation')
                ylabel('Number of Readings in % of Number of Data Points')
                xlim([-4 4]);
                ylim([0 50]);
                title(strrep(['Histogram of ALG+GRBF Validation Model for %s',7],'%s',loadlist(k3)));
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
            plotResPages(seriesvalid, targetRes2valid, loadCapacities, stdDevPercentCapacity2valid )
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

        if excel_FLAG == 1 && balCal_FLAG == 2
            disp(' ');
            disp('ALG+GRBF VALIDATION MODEL GLOBAL LOAD APPROXIMATION: Check VALID_AOX_GLOBAL_GRBF_RESULT.csv file');
            disp(' ');

            filename = 'VALID_AOX_GLOBAL_GRBF_RESULT.csv';
            Z = aprxINminGZ2valid;
            xlRange = 'matrixcolumnlabels(1)1:matrixcolumnlabels(dimFlag)numpts';
            xlswrite(filename,Z,xlRange)
        end
    end
end

if balApprox_FLAG == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                        APPROXIMATION SECTION      AJM 6/29/17           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %DEFINE THE PRODUCTION CSV INPUT FILE AND SELECT THE RANGE OF DATA VALUES TO READ
    %
    %     inputApprox_balCal = 'MK14C-ChkLds-Ames2011-Meade-8D_voltage.csv';
    %     excessVecapprox =         csvread(inputApprox_balCal,19,12,'M20..T143');
    load(out.savePathapp,'-mat');

    nseriesapprox=max(seriesapprox);
    lasttareapprox=nseriesapprox;
    [sharedvals,indexLocalZeroapprox]=intersect(seriesapprox, [1:max(seriesapprox)],'stable'); %Create index
    
%     % testing
%     % nseriesapprox = nseriesvalid;
%     nseriesapprox = 11;
%     indexLocalZeroapprox(1) = 1;
%     indexLocalZeroapprox(2) = 10;
%     indexLocalZeroapprox(3) = 19;
%     indexLocalZeroapprox(4) = 28;
%     indexLocalZeroapprox(5) = 37;
%     indexLocalZeroapprox(6) = 46;
%     indexLocalZeroapprox(7) = 55;
%     indexLocalZeroapprox(8) = 66;
%     indexLocalZeroapprox(9) = 77 ;
%     indexLocalZeroapprox(10) = 88;
%     indexLocalZeroapprox(11) = 99;
%     indexLocalZeroapprox(12) = 110;
%     loadCapacitiesapprox(1) = 2500;
%     loadCapacitiesapprox(2) = 2500;
%     loadCapacitiesapprox(3) = 1250;
%     loadCapacitiesapprox(4) = 1250;
%     loadCapacitiesapprox(5) = 5000;
%     loadCapacitiesapprox(6) = 700;

    % testing

    % num of data points
    numptsapprox = length(excessVecapprox);

    dimFlagapprox = length(excessVecapprox(1,:));

    %find the average natural zeros (also called global zeros)
    globalZerosapprox = mean(natzerosapprox);

    %%% make an array out of the globalZerosapprox vector
    for i=1:numptsapprox
        globalZerosAllPointsapprox(i,:) = globalZerosapprox;
    end

    % Subtract the Global Zeros from the Inputs %%%%%%%%%%
    for k=1:dimFlagapprox
        dainputsapprox(:,k) = excessVecapprox(:,k)-globalZerosAllPointsapprox(:,k);
        dalzapprox(:,k) = globalZerosAllPointsapprox(:,k)-globalZerosAllPointsapprox(:,k);
    end

    % Build the Algebraic Model
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
    comGZapprox= zeros(nterms+1,1);
    for i=1:dimFlag
        biggeeapprox(:,i) = 0;
    end


    [comINapprox,comLZapprox,comGZapprox,uncert_comINapprox]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numptsapprox,0,1,dainputsapprox,dalzapprox,biggeeapprox);

    for i=1:nterms+1
        xapprox(i,:) = xcalib(i,:);
    end

    %LOAD APPROXIMATION
    %define the approximation for inputs minus global zeros
    interceptsapprox = -(comGZapprox'*xapprox);
    aprxINapprox = (xapprox'*comINapprox)';        %to find ?? AJM111516

    for m=1:length(aprxINapprox)
        %%%% 3/23/18 Remove Intercepts %%%%%
        %        aprxINminGZapprox(m,:) = aprxINapprox(m,:)+interceptsapprox;
        aprxINminGZapprox(m,:) = aprxINapprox(m,:);

        stdevchecktestapprox(m,:) = aprxINminGZapprox(m,:);
    end

    %Run function to calculate uncertainty on loads from approx
[combined_uncertapprox,tare_uncertapprox, FL_uncertapprox]=uncert_prop(xapprox,fxcalib_ci,comINapprox,dimFlag,uncert_comINapprox,indexLocalZeroapprox,lasttareapprox,nterms,aprxINapprox,seriesapprox);





    
%     % SOLVE FOR TARES
%     for i=1:nseriesapprox
%         zoopapprox = zeros(length(excessVecapprox(:,1)),dimFlag);
%         
%         kx=indexLocalZeroapprox(i)-1;
%         
%         daemmlengthapprox = indexLocalZeroapprox(i+1) - indexLocalZeroapprox(i);
%         
%         stdevchecktestapprox = zeros(daemmlengthapprox,dimFlag);   %% ajm 7_18_18
%         
%         for m= indexLocalZeroapprox(i): indexLocalZeroapprox(i+1)-1
%             zoopapprox(m-kx,:) = aprxINminGZapprox(m,:);
%             stdevchecktestapprox(m-kx,:) = aprxINminGZapprox(m,:);
%         end
%         
%         zapapprox(i,:) = mean(zoopapprox)*numptsapprox/(indexLocalZeroapprox(i+1)-indexLocalZeroapprox(i));
%         
%         zapstdevapprox(i,:) =  std(stdevchecktestapprox);  %%% ajm 7_17_18
%     end
%     
%     for i=1:nseriesapprox
%         for j= 1: dimFlag
%             stdevfilterapprox(i,j) = 100.0*zapstdevapprox(i,j)/loadCapacitiesapprox(1,j); %% ajm 7_17_18
%             
%             if stdevfilterapprox(i,j) > 0.25
%                 zapapprox(i,j) = aprxINminGZapprox(indexLocalZeroapprox(i),j);
%             end
%         end
%         for m=indexLocalZeroapprox(i):indexLocalZeroapprox(i+1)-1
%             taretalapprox(m,:) = zapapprox(i,:);
%         end
%     end
    

    disp(' ');
    disp('%%%%%%%%%%%%%%%%%');
    disp('  ');
    disp('Approximation data file read =');
    disp(out.savePathapp);
    disp('  ');

    if excel_FLAG == 1
        disp('  ');
        disp('ALG MODEL GLOBAL LOAD APPROXIMATION: Check APPROX_AOX_GLOBAL_ALG_RESULT.csv file');
        disp(' ');

        filename = 'APPROX_AOX_GLOBAL_ALG_RESULT.csv';
        Z = aprxINminGZapprox;
        xlRange = 'A1:JnumBasis';
        xlswrite(filename,Z,xlRange)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    RBF SECTION FOR APPROXIMATION     AJM 6/29/17                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %goal to use centers, width and coefficients to approxate parameters against
    %independent data

    aprxINminGZ2approx = aprxINminGZapprox;
    if balCal_FLAG == 2

        etaHistapprox = cell(numBasis,1);
        aprxINminGZ_Histapprox = cell(numBasis,1);

        etaGZapprox = dot(dalzapprox-dainputsapprox,dalzapprox-dainputsapprox);

        for u=1:numBasis
            for s=1:length(excessVec(1,:)) % loops through the 8 components

                centerIndexLoop(s) = centerIndexHist(u,s); %Have to use the history or it gets overwritten

                for r=1:length(excessVecapprox(:,1))
                    etaapprox(r,s) = dot(dainputsapprox(r,:)-dainputscalib(centerIndexLoop(s),:),dainputsapprox(r,:)-dainputscalib(centerIndexLoop(s),:));
                end

                rbfINminGZapprox(:,s)=exp(etaapprox(:,s)*log(abs(wHist(u,s))));

                rbfc_INminGZapprox(:,s) = cHist(u,s)*rbfINminGZapprox(:,s);
            end

            %update the approximation
            aprxINminGZ2approx = aprxINminGZ2approx+rbfc_INminGZapprox;
            aprxINminGZ_Histapprox{u} = aprxINminGZ2approx;

            % SOLVE FOR TARES BY TAKING THE MEAN
            for i=1:nseriesapprox
                zoop2approx = zeros(length(excessVecapprox(:,1)),dimFlag);

                kx=indexLocalZeroapprox(i)-1;

                daemmlengthapprox = indexLocalZeroapprox(i+1) - indexLocalZeroapprox(i);

                stdevchecktest2approx = zeros(daemmlengthapprox,dimFlag);   %% ajm 7_18_18

                for m= indexLocalZeroapprox(i): indexLocalZeroapprox(i+1)-1
                    zoop2approx(m-kx,:) = aprxINminGZ2approx(m,:);
                    stdevchecktest2approx(m-kx,:) = aprxINminGZ2approx(m,:);
                end

                zap2approx(i,:) = mean(zoop2approx)*numptsapprox/(indexLocalZeroapprox(i+1)-indexLocalZeroapprox(i));

                zapstdev2approx(i,:) =  std(stdevchecktest2approx);  %%% ajm 7_17_18
            end

            for i=1:nseriesapprox
                for j= 1: dimFlag
                    stdevfilter2approx(i,j) = 100.0*zapstdev2approx(i,j)/loadCapacitiesapprox(1,j); %% ajm 7_17_18

                    if stdevfilter2approx(i,j) > 0.25
                        zap2approx(i,j) = (aprxINminGZ2approx(indexLocalZeroapprox(i),j)+ aprxINminGZ2approx(indexLocalZeroapprox(i+1)-1,j))/2.0;
                    end
                end

                for m=indexLocalZeroapprox(i):indexLocalZeroapprox(i+1)-1
                    taretalGRBFapprox(m,:) = zap2approx(i,:);
                end
            end
        end

        if excel_FLAG == 1 && balCal_FLAG == 2
            disp(' ');
            disp('ALG+GRBF MODEL GLOBAL LOAD APPROXIMATION: Check APPROX_AOX_GLOBAL_GRBF_RESULT.csv file');
            disp(' ');

            filename = 'APPROX_AOX_GLOBAL_GRBF_RESULT.csv';
            Z = aprxINminGZ2approx;
            xlRange = 'A1:JnumBasis';
            xlswrite(filename,Z,xlRange)
        end
    end
end

disp('  ')
disp('Calculations Complete.')

% Tidy up the Workspace
%clearvars -except CALIB_AOX_GLOBAL_ALG_RESULT CALIB_AOX_GLOBAL_GRBF_RESULT VALID_AOX_GLOBAL_ALG_RESULT VALID_AOX_GLOBAL_GRBF_RESULT APPROX_AOX_GLOBAL_ALG_RESULT APPROX_AOX_GLOBAL_GRBF_RESULT APPROX_AOX_COEFF_MATRIX  APPROX_AOX_GRBF_Centers APPROX_AOX_GRBF_coeffs APPROX_AOX_GRBF_ws BALFIT_DATA_REDUCTION_MATRIX_IN_AMES_FORMAT  OUTLIER_ROWS OUTLIER_ROWS2 zeroed_series zeroed_targetMatrix zeroed_excessVec outlierIndices   % 5_11_18 ajm

%
% Copyright ©2016 Andrew Meade and Ali Arya Mokhtarzadeh.  All Rights Reserved.
%

%toc
function [out]=zapFinder(comINminLZ,target,series_complete,excessVec0_complete,targetMatrix0_complete,globalZerosAllPoints0,localZerosAllPoints,dimFlag,model_FLAG,nterms,numpts,lasttare,nseries)
xcalib = comINminLZ\target;                    % '\' solution of Ax=b

%    xcalib = lsqminnorm(comINminLZ',targetMatrix);             % alternate solution method
%    xcalib = pinv(comINminLZ')*targetMatrix;             % alternate solution method
%    xcalib = pinv(comINminLZ',1e-2)*targetMatrix;             % alternate solution method

%%  Creates Matrix for the volts to loads
%  ajm for the users to view 5/11/18
%      APPROX_AOX_COEFF_MATRIX = xcalib;
%
%     if LHS_Flag == 1
%         x_all(:,:,lhs) = xcalib;
%     end


% if LHS_Flag == 1
%     xcalib = mean(x_all,3);
%     xcalib_std = std(x_all,[],3);
% end

excessVec0 = excessVec0_complete;
excessVec = excessVec0;
targetMatrix = targetMatrix0_complete;
series = series_complete;

% %START issues with resampling
% %Moved outside of loop
% prevVal = 1;
% for i=1:length(series)
%     if series(i)==prevVal
%         localZeros(prevVal,:) = excessVec0(i,:);
%         localZerosAllPoints(i,:) = localZeros(prevVal,:);
%         prevVal = prevVal+1;
%     else
%         counter(prevVal) = i;
%         localZerosAllPoints(i,:) = localZeros(prevVal-1,:);
%     end
%     %        globalZerosAllPoints(i,:) = globalZeros;
% end
%
% %add '1' to every value of counter to find the indices of the local zeros
% indexLocalZero = counter(:)+1;

for i=1:dimFlag
    %     dainputs2(:,i) = excessVec0(:,i)-globalZeros(i);
    %     %
    %     dalz2(:,i) = localZerosAllPoints(:,i)-globalZeros(i);
    %
    biggee(:,i) = 0;
    %
end
dainputs2=excessVec0-globalZerosAllPoints0;
dalz2=localZerosAllPoints-globalZerosAllPoints0;
% Call the Algebraic Subroutine
%

%
[comIN,comLZ,comGZ]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numpts,series,lasttare,dainputs2,dalz2,biggee);
%%
%%
%%


%% Effectively removes the original tare values so we can calculate the averages
for i=1:lasttare
    comIN(nterms+i,:) = 0;
    comLZ(nterms+i,:) = 0;
    comGZ(nterms+i,:) = 0;
end


%%
%QUESTION
% for loopk=1:numpts
%     comLZ(nterms+series(loopk),loopk) = 1.0;
% end

%%
%%


comINminLZ = comIN-comLZ;

%APPROXIMATION
%define the approximation for inputs minus local zeros
aprxINminLZ = comINminLZ'*xcalib;


%DIFFERENT APPROXIMATION
%define the approximation for inputs minus global zeros
intercepts = -(comGZ'*xcalib);
aprxIN = (xcalib'*comIN)';
aprxLZ = (xcalib'*comLZ)';       %to find tares AAM042016
for m=1:length(aprxIN)
    %    aprxINminGZ(m,:) = aprxIN(m,:)+intercepts;  %ajm zap 3/23/18
    %    aprxLZminGZ(m,:) = aprxLZ(m,:)+intercepts;  %ajm zap 3/23/18

    aprxINminGZ(m,:) = aprxIN(m,:);
    aprxLZminGZ(m,:) = aprxLZ(m,:); %to find tares AAM042016

    %subtracts the targets from the global approx
    %This will be averaged over the individual series to give the tares

    checkit(m,:) = aprxINminGZ(m,:)-targetMatrix(m,:);
end

%create new vectors of series and checkit sorted by series in assending
%order
[series_sort,seriesI]=sort(series);
checkit_sort=checkit(seriesI,:);
localZerosAllPoints_sort=localZerosAllPoints(seriesI,:);
aprxLZminGZ_sort=aprxLZminGZ(seriesI,:);


%count number of entries for each series
% SOLVE FOR TARES BY TAKING THE MEAN
for i=1:nseries
    index_count(i)=sum(series==i); %Find number of datapoints in each series
    fzap(i,:)=mean(checkit_sort(sum(index_count(1:i-1))+1:sum(index_count),:)); %Take mean of checkit for each series
    aprxLZminGZ_series(i,:)=aprxLZminGZ_sort(sum(index_count(1:i-1))+1,:); %Find matrix of Local zero approximations, 1 line for each series
end
out=[fzap;aprxLZminGZ_series;xcalib];
end

function [combined_uncert,tare_uncert, FL_uncert]=uncert_prop(xcalib,fxcalib_ci,comIN,dimFlag,uncert_comIN,indexLocalZero,lasttare,nterms,aprxIN,series)
%Function calculates uncertainty in load from uncertainty in coefficients and input voltages

%Inputs:
%xcalib: Coefficients from calibration
%fxcalib_ci: Confidence interval for coefficients
%comIN: Produced by balCal_algEquations3: Voltages provided put into matrix
%to multiply by coefficients for load output
%dimFlag: Number of dimensions (channels)
%uncert_comIN: Produced by balCal_algEquations3: Voltages provided put into
%matrix to multiply by coefficients for partial derivatives
%indexLocalZero: Index of where which datapoints are tare loads
%lasttare: Number of the final series
%nterms: number of coefficients used for model
%aprxIN: approximation matrix of loads

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
for i=1:size(xcalib,1)
    for j=1:size(xcalib,2)
        for k=1:2
            error(k,i,j)=abs(fxcalib_ci(k,i,j)-xcalib(i,j));
        end
        xcalib_error(i,j)=max(error(:,i,j));
    end
end
comIN_square=comIN.^2; %Square input matrix
xcalib_error_square=xcalib_error.^2; %square error matrix
coeff_uncert_square=(xcalib_error_square'*comIN_square)'; %Matrix for error in each calibration point from uncertainty in coefficients(+/-)
%END:  coeff uncertainty propagation
%ADDED 23 Jan 19 JRP: Analytical calc of uncert in loads due to uncertainty in read
%voltages
uncert=1; %Trust all channels down to 1 microvolt
for i=1:lasttare
    for j=1:dimFlag
        uncert_comIN(nterms+i,:,j) = 0;
    end
end
uncert_comIN_use=uncert_comIN(1:size(xcalib,1),:,:);
volt_uncert_square=zeros(size(aprxIN));
for i=1:dimFlag
    partial(:,:,i)=uncert_comIN_use(:,:,i)'*xcalib;
    uncert_channel_square(:,:,i)=((partial(:,:,i)).^2).*uncert^2;
    volt_uncert_square=volt_uncert_square+(uncert_channel_square(:,:,i));
end

%Combine uncertainty from coeff and input voltages for 1 total uncertainty
%value for every datapoint:
combined_uncert=(coeff_uncert_square+volt_uncert_square).^(.5);
tare_uncert=combined_uncert(indexLocalZero(1:(numel(indexLocalZero)-1)),:); %Uncert for tare loads

%FL Uncert:
for i=1:numel(series)
    if i==indexLocalZero(series(i))
        FL_uncert(i,:)=combined_uncert(i,:);
    else
        FL_uncert(i,:)=(combined_uncert(i,:).^2+tare_uncert(series(i)).^2).^.5;
    end

end

%END added for uncert

end
