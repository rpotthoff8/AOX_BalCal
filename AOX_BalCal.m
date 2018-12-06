% 10/02/18
%
% Copyright 2017 Andrew Meade, Ali Arya Mokhtarzadeh and Javier Villarreal.  All Rights Reserved.
%
%balanceCalibration_with_RBF_8D.m
%requires "balCal_meritFunction.m" to run
%input file: "BuffetBalance-CalDataOfOct2015-3F3M.csv"
%output file: "balCal_output_ALGB.xls"
%output file: "balCal_output_GRBF.xls"
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize the workspace
clc;
clearvars;
close all;
workspace;
disp('Copyright 2017 Andrew Meade, Ali Arya Mokhtarzadeh and Javier Villarreal.  All Rights Reserved.')
% The mean of the approximation residual (testmatrix minus local approximation) for each section is taken as the tare for that channel and section. The tare is subtracted from the global values to make the local loads. The accuracy of the validation, when compared to the known loads, is very sensitive to the value of the tares (which is unknown) and NOT the order of the calibration equations.
% Because of measurement noise in the voltage the APPROXIMATION tare is computed by post-processing. The average and stddev is taken of all channels per section. If the stddev is less than 0.25% of the capacity for any station the tare is equal to the average for that channel. If the stddev is greater than 0.25% then the approximation at the local zero is taken as the tare for that channel and section. The global approximation is left alone but the tare is subtracted from the values to make the local loads. Line 3133.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       USER INPUT SECTION
%
out = AOX_GUI;
%out = AOX_GUIv12;
if out.cancel == 1
    return
end
%TO SELECT Algebraic Model                         set balCal_FLAG = 1;
%TO SELECT Algebraic and GRBF Model                set balCal_FLAG = 2;
balCal_FLAG = out.grbf;
%
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
%TO OUTPUT DATA TO EXCEL                           set excel_FLAG = 1;
excel_FLAG = 1;
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
%DEFINE THE CSV INPUT FILE AND SELECT THE RANGE OF DATA VALUES TO READ
% inputFile_balCal = 'MK14C-CalData-ABC10-Meade-Iter-8D.csv';
% loadCapacities =    csvread(inputFile_balCal,5,4,'E6..L6');
% natzeros =          csvread(inputFile_balCal,8,12,'M9..T12');
% series =            csvread(inputFile_balCal,15,2,'C16..C1804');
% targetMatrix0 =      csvread(inputFile_balCal,15,4,'E16..L1804');
% excessVec0 =         csvread(inputFile_balCal,15,12,'M16..T1804');
load(out.savePathcal,'-mat');
%
%csvread Read a comma separated value file.
%Adapted from published MathWorks documentation by Mokhtarzadeh, 042516
%
%   M = csvread('FILENAME',R,C,RNG) reads a comma separated value formatted
%   file FILENAME, starting at row R and column C where R=0 and C=0
%   specifies the first value in the file (a zero based index), and reads
%   only the range specified by RNG using spreadsheet notation as in
%   RNG = 'A1..B2' where A1 is the upper-left corner of the data to be
%   read and B2 is the lower-right corner.
%   The file can only contain numeric values.
%
%                       END USER INPUT SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

excessVec0_complete = excessVec0;
targetMatrix0_complete = targetMatrix0;
series_complete = series;

iloadCapacities = max(abs(excessVec0_complete));

%%%%%%%% 6_14_18 ajm
warning('off','all');
warning;
%%%%%%%%

disp('  ')
disp('Starting Calculations')

if approach_FLAG == 1
    disp(' ');
    disp('Using the Indirect Approach for Calibration');
    
elseif approach_FLAG == 0
    disp(' ');
    disp('Using the Direct Approach for Calibration');
end

%%%% 6_14_18 ajm
matrixcolumnlabels = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};
matrixcolumnlabels(:) = strrep(matrixcolumnlabels(:),'''',''); %get rid of single quotes

loadlist = {'N1','N2','S1','S2','RM','AF','PLM', 'PCM', 'MLM', 'MCM'};
voltagelist = {'R1','R2','R3','R4','R5','R6','rPLM','rPCM','rMLM','rMCM'};
%%%%


if LHS_Flag == 0
    numLHS = 1;
end
lhs_check = zeros(length(excessVec0_complete),1);
ind = 1:length(excessVec0_complete(:,1));


if LHS_Flag == 1
    disp('  ')
    disp('Number of LHS Iterations Selected =')
    disp(numLHS);
end

for lhs = 1:numLHS
    
    if LHS_Flag == 1
        [series,targetMatrix0,excessVec0,sample] = AOX_LHS(series_complete,targetMatrix0_complete,excessVec0_complete,LHSp);
        lhs_check(sample) = 1;
        ind(find(lhs_check-1)); % This line outputs which data points haven't been sampled yet
        pct_sampled = sum(lhs_check)/length(lhs_check); % This line outputs what percentage of points have been sampled
    end
    
    %% Keep the excessVec0 in reserve for the other subroutines
    % excessVec = excessVec0;
    %%
    
    dimFlag = length(excessVec0(1,:));
    
    % num of data points
    numpts = length(excessVec0(:,1));
    
    
    %find the average natural zeros (also called global zeros)
    globalZeros = mean(natzeros);
    %globalZeros = natzeros;
    
    
    %load capacities
    loadCapacities(loadCapacities == 0) = realmin;
    
    %find number of series; this will tell us the number of tares
    nseries = max(series);
    
    
    
    %find zero points of each series and number of points in a series
    %localZerosAllPoints is the same as localZeroMatrix defined in the RBF
    %section - see if this can be cleaned up! is if/else necessary?  AAM042016
    prevVal = 1;
    for i=1:length(series)
        if series(i)==prevVal
            localZeros0(prevVal,:) = excessVec0(i,:);
            localZerosAllPoints0(i,:) = localZeros0(prevVal,:);
            prevVal = prevVal+1;
        else
            counter(prevVal) = i;
            localZerosAllPoints0(i,:) = localZeros0(prevVal-1,:);
        end
        globalZerosAllPoints0(i,:) = globalZeros;
    end
    
    %add '1' to every value of counter to find the indices of the local zeros
    indexLocalZero = counter(:)+1;
    
    
    %%% 5/16/18
    excessVec = excessVec0 - globalZerosAllPoints0;
    localZerosAllPoints =  localZerosAllPoints0 - globalZerosAllPoints0;
    %%%
    
    
    disp('  ')
    disp('Working ...')
    
    
    
    % @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    %%
    %
    %
    % Copyright ©2016 Andrew Meade and Ali Arya Mokhtarzadeh.  All Rights Reserved.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 VOLTAGE TO LOAD (DIRECT) SECTION      AJM 8/3/17           %
    %Use the measured to best voltages mapping to find the calibration coefficients
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    targetMatrix = targetMatrix0;
    excessVec = excessVec0; % Direct approach uses the measured voltage
    
    % Indirect approach uses the modeled voltage
    if approach_FLAG == 1
        
        if balCal_FLAG == 2
            excessVec0 = qtaprxINminGZ2 + globalZerosAllPoints0;
        else
            excessVec0 = qtaprxINminGZ + globalZerosAllPoints0;
        end
        
    end
    
    %
    %
    
    %find the average natural zeros (also called global zeros)
    globalZeros = mean(natzeros);
    %globalZeros = zeros(length(excessVec(1,:)),1);
    
    
    %load capacities
    loadCapacities(loadCapacities == 0) = realmin;
    
    %find number of series; this will tell us the number of tares
    nseries = max(series);
    
    %find zero points of each series and number of points in a series
    %localZerosAllPoints is the same as localZeroMatrix defined in the RBF
    %section - see if this can be cleaned up! is if/else necessary?  AAM042016
    %
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
    
    
    %% Subtract the Global Zeros from the Inputs and Local Zeros %%%%%%%%%%
    
    for i=1:dimFlag
        dainputs(:,i) = excessVec0(:,i)-globalZeros(i);
        %
        dalz(:,i) = localZerosAllPoints(:,i)-globalZeros(i);
        %
        biggee(:,i) = 0;
        %
    end
    %%%%%%%%%%%%
    
    
    %find the number of points in each series and the number of series
    for j = 1:length(counter)-1
        numPoints(j) = counter(j+1)-counter(j);
        numSeries(j,:) = j;
    end
    
    
    %%
    %% Build the Algebraic Model
    %%
    
    lasttare = series(numpts);
    
    %% Full Algebraic Model
    if model_FLAG == 1
        nterms = 2*dimFlag*(dimFlag+2);
    end
    
    %% Truncated Algebraic Model
    if model_FLAG == 2;
        nterms = dimFlag*(dimFlag+3)/2;
    end
    
    %% Linear Algebraic Model
    if model_FLAG == 3;
        nterms = dimFlag;
    end
    
    
    
    % Call the Algebraic Subroutine
    %
    
    %
    [comIN,comLZ,comGZ]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numpts,series,lasttare,dainputs,dalz,biggee);
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
    
    for loopk=1:numpts
        comLZ(nterms+series(loopk),loopk) = 1.0;
    end
    
    %%
    %%
    
    
    comINminLZ = comIN-comLZ;
    
    %SOLUTION
    xcalib = comINminLZ'\targetMatrix;                    % '\' solution of Ax=b
    
    %BOOTSTRAP COEFFICIENTS
    coeff_solve=@(comIN,target)comIN\target;
    xcalib_ci=bootci(2000,coeff_solve,comINminLZ',targetMatrix);
    % END Bootstrap
    
%     jackstat_xcalib=jackknife(coeff_solve,comINminLZ',targetMatrix);
%     xbar=mean(jackstat_xcalib);
% %     V=(1/(size(jackstat_xcalib,2)-1))*(sum((jackstat_xcalib-xbar).^2));
%     sum=0;
%     for i=1:size(jackstat_xcalib,2)
%         sum=sum+(jackstat_xcalib(i,:)-xbar).^2;
%     end
%     V2=(1/(size(jackstat_xcalib,2)-1))*sum;
%     jackstat_lower=xbar-1.960*sqrt((1/size(jackstat_xcalib,2))*V2);
    
    

    
    %    xcalib = lsqminnorm(comINminLZ',targetMatrix);             % alternate solution method
    %    xcalib = pinv(comINminLZ')*targetMatrix;             % alternate solution method
    %    xcalib = pinv(comINminLZ',1e-2)*targetMatrix;             % alternate solution method
    %% 5/17/18
    for i=1:nterms+1
        xvalid(i,:) = xcalib(i,:);
    end
    
    %%  Creates Matrix for the volts to loads
    %  ajm for the users to view 5/11/18
    %      APPROX_AOX_COEFF_MATRIX = xcalib;
    
    for i=1:nterms
        xapproxer(i,:) = xcalib(i,:);
    end
    
    xapproxer(nterms+1,:) = 0.0;
    
    filename = 'APPROX_AOX_COEFF_MATRIX.csv';
    Z = xapproxer;
    %     xlRange = 'matrixcolumnlabels(1)1:matrixcolumnlabels(dimFlag)nterms';
    xlRange=char(strcat(matrixcolumnlabels(1),'1:',matrixcolumnlabels(dimFlag),num2str(nterms)));
    %xlswrite(filename,Z,xlRange)
    %%
    
    
    if LHS_Flag == 1
        x_all(:,:,lhs) = xcalib;
    end
end

if LHS_Flag == 1
    xcalib = mean(x_all,3);
    xcalib_std = std(x_all,[],3);
end

excessVec0 = excessVec0_complete;
excessVec = excessVec0;
targetMatrix = targetMatrix0_complete;
series = series_complete;

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

for i=1:dimFlag
    dainputs2(:,i) = excessVec0(:,i)-globalZeros(i);
    %
    dalz2(:,i) = localZerosAllPoints(:,i)-globalZeros(i);
    %
    biggee(:,i) = 0;
    %
end

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

for loopk=1:numpts
    comLZ(nterms+series(loopk),loopk) = 1.0;
end

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




%%
%% SOLVE FOR TARES BY TAKING THE MEAN
for i=1:nseries
    zoop = zeros(length(excessVec(:,1)),dimFlag);
    kx=indexLocalZero(i)-1;
    
    for m=indexLocalZero(i):indexLocalZero(i+1)-1
        zoop(m-kx,:) = checkit(m,:);
    end
    
    zap(i,:) = mean(zoop)*numpts/(indexLocalZero(i+1)-indexLocalZero(i));
    %Bootstrap residuals
    meanResidual=@(z)mean(z);
    zap_ci(2*i-1:2*i,:)=bootci(10000,{meanResidual,zoop(1:indexLocalZero(i+1)-indexLocalZero(i),:)},'alpha',.05);
    %end bootstrap
    %jackknife residuals
    jackstat_zap=jackknife('mean',zoop(1:indexLocalZero(i+1)-indexLocalZero(i),:));
    jack_zap=mean(jackstat_zap);
    sum=0;
    for j=1:size(jackstat_zap,2)
        sum=sum+(jackstat_zap(j,:)-jack_zap).^2;
    end
    V2=(1/(size(jackstat_zap,2)-1))*sum;
    jackstat_lower(i,:)=jack_zap-1.960*sqrt((1/size(jackstat_zap,2))*V2);
    %end jackknife
end

%%
%%
for i=1:nseries
    for m=indexLocalZero(i):indexLocalZero(i+1)-1
        taretal(m,:) = zap(i,:);
    end
end
%%
%%


%RESIDUAL
targetRes = targetMatrix+taretal-aprxINminGZ;      %0=b-Ax

%find the sum of squares of the residual using the dot product
resSquare = dot(targetRes,targetRes)';
%AAM note to self - in matlab, diag(A'*A) is the same as dot(A,A)'


%---------------------------------------------------------------
%---------------------------------------------------------------
% Identify Outliers After Filtering
% (Threshold approach) ajm 8/2/17
%

if balOut_FLAG == 1 % !!!
    
    
    detect_targetRes = targetRes;
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    
    % Use the modeled input for the rest of the calculations
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    
    for n = 1:dimFlag
        normtargetRes(:,n) = detect_targetRes(:,n)/loadCapacities(n);
    end
    
    
    out_meanValue = mean(normtargetRes);
    
    
    % Identify outliers. They are considered outliers if the residual
    % is more than 3 standard deviations as % of capacity from the mean.
    
    out_standardDev = std(normtargetRes);
    
    numSTD = 3.0; % Whatever you want.
    
    thresholdValue = numSTD * (out_standardDev) - out_meanValue;
    
    %      if approach_FLAG == 0
    
    for n = 1:dimFlag
        
        if thresholdValue(1,n) <= 0.0025
            thresholdValue(1,n) = 0.0025;
        end
        
    end
    
    %      end
    
    
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
    
    
    %
    %% Use the reduced input and target files
    %
    if zeroed_FLAG == 1
        
        
        disp(' ************************************************************************ ');
        disp('Find the reduced data in zeroed_targetMatrix and zeroed_excessVec');
        disp(' ************************************************************************ ');
        
        
        % Mark the outliers values, store and use for recalculation:
        
        
        zeroed_targetMatrix = targetMatrix;
        zeroed_excessVec = excessVec0;
        zeroed_series = series;
        zeroed_numpts = numpts;
        
        
        %        for k1 = 1:num_outliers;
        %                   zeroed_targetMatrix(OUTLIER_ROWS(k1,1),:) = 10000;
        %                   zeroed_excessVec(OUTLIER_ROWS(k1,1),:) = 10000;
        %        end
        
        zeroed_targetMatrix(OUTLIER_ROWS,:) = [];
        zeroed_excessVec(OUTLIER_ROWS,:) = [];
        
        targetMatrix = unique(zeroed_targetMatrix,'rows');
        excessVec0 = unique(zeroed_excessVec,'rows');
        
        %        for k1 = 1:num_outliers;
        %                   zeroed_series(OUTLIER_ROWS(k1,1)) = [];
        %       end
        
        zeroed_series(OUTLIER_ROWS) = [];
        
        
        numpts =  zeroed_numpts - num_outliers;
        
        targetMatrix = zeroed_targetMatrix;
        excessVec0 = zeroed_excessVec;
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
            globalZerosAllPoints(i,:) = globalZeros;
        end
        
        %add '1' to every value of counter to find the indices of the local zeros
        indexLocalZero = counter(:)+1;
        
        
        %% Subtract the Global Zeros from the Inputs and Local Zeros %%%%%%%%%%
        
        for i=1:dimFlag
            dainputs(:,i) = excessVec0(:,i) - globalZerosAllPoints(i);
            %
            dalz(:,i) = localZerosAllPoints(:,i) - globalZerosAllPoints(i);
            %
        end
        %%%%%%%%%%%%
        
        for i=1:dimFlag
            biggee(:,i) = 0;
        end
        [comIN,comLZ]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numpts,series,lasttare, dainputs, dalz, biggee);
        
        comINminLZ = comIN-comLZ;
        
        %SOLUTION
        xcalib = pinv(comINminLZ')*targetMatrix;             % alternate solution method
 
        
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
        
        
        
        %%
        %% SOLVE FOR TARES BY TAKING THE MEAN
        for i=1:nseries
            zoop = zeros(numpts,dimFlag);
            kx=indexLocalZero(i)-1;
            
            for m=indexLocalZero(i):indexLocalZero(i+1)-1
                zoop(m-kx,:) = checkit(m,:);
            end
            
            zap(i,:) = mean(zoop)*numpts/(indexLocalZero(i+1)-indexLocalZero(i));
            
        end
        
        %%
        %%
        for i=1:nseries
            for m=indexLocalZero(i):indexLocalZero(i+1)-1
                taretal(m,:) = zap(i,:);
            end
        end
        
        %%
        %%
        
        
        %RESIDUAL
        targetRes = targetMatrix+taretal-aprxINminGZ;      %0=b-Ax
        
    end
    %
    %% End of zeroed input and target files
    %
    
    
    
end %!!!

%
%---------------------------------------------------------------
%---------------------------------------------------------------


%

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




% *****

for n = 1:length(indexLocalZero)-1
    taresALGB(n,:) = taretal(indexLocalZero(n),:);
end

% *****


%%%%%%%%


%OUTPUT HISTOGRAM PLOTS
if hist_FLAG == 1;
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


%%%%%%%%



%
%START PRINT OUT PERFORMANCE INFORMATION TO THE SCREEN
if print_FLAG == 1;
    %
    
    
    %% END Direct or Indirect Approach to Calibration
    
    
    
    
    %% Identify the Possible Outliers
    if balOut_FLAG == 1;
        disp(' ***** ');
        disp(' ');
        disp('Number of Outliers =');
        disp(num_outliers);
        disp('Outliers % of Data =');
        disp(prcnt_outliers);
    end
    
    
    
    %% Recalculated Calibration with Reduced Matrices
    if zeroed_FLAG == 1
        
        disp(' ************************************************************************ ');
        disp('Find the reduced data in zeroed_targetMatrix and zeroed_excessVec');
        disp(' ************************************************************************ ');
        
    end
    
    %
    disp('  ');
    disp('ALG CALIBRATION MODEL GLOBAL LOAD APPROXIMATION FILE: CALIB_AOX_GLOBAL_ALG_RESULT.csv');
    % CALIB_AOX_GLOBAL_ALG_RESULT = aprxINminGZ;
    disp(' ');
    
    
    %%%%%
    filename = 'CALIB_AOX_GLOBAL_ALG_RESULT.csv';
    Z = aprxINminGZ;
    %        xlRange = 'A1:Jnumpts';
    xlRange=char(strcat('A1:J',num2str(numpts)));
    %xlswrite(filename,Z,xlRange)
    %%%%%
    
    
    
    
    
    %%%%%%% 6_14_18 ajm
    
    calib_twoSigmaALGB = standardDev'.*2;
    calib_algebraic_2Sigma = array2table(calib_twoSigmaALGB,'VariableNames',loadlist(1:dimFlag))
    
    numSeriesNames = num2str(numSeries);
    numSeriesNameCell = cellstr(numSeriesNames);
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
    
    %% Prints the minmaxband
    calib_alg_per_minmaxband = array2table(theminmaxband,'VariableNames',loadlist(1:dimFlag))
    
    %%%%%%%%%
    
    
    
    if excel_FLAG == 1;
        % Output results to an excel file
        
    end
    
end


if res_FLAG == 1
    figure('Name','Algebraic Model Calibration; Residuals of Load Versus Data Point Index','NumberTitle','off')
    plotResPages(series, targetRes, loadCapacities, stdDevPercentCapacity )
    %    hold off
end


%END PRINT OUT PERFORMANCE INFORMATION TO THE SCREEN


% Start GRBF Option
%
if balCal_FLAG == 2
    %
    % Copyright ©2016 Andrew Meade and Ali Arya Mokhtarzadeh.  All Rights Reserved.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             VOLTAGE TO LOAD (DIRECT) - RBF SECTION                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %goal to minimize: minimize the sum of the squares (dot product) of each of the 8
    %residual vectors 'targetRes' 'target1' ... 'target8'
    %dt1 = dot(target1,target1);
    %find centers by finding the index of max residual, using that index to
    %subtract excess(counter)-excess(indexMaxResid) and then taking the dot
    %product of the resulting column vector
    %
    %%
    
    for i=1:dimFlag
        for s=1:length(series)
            targetRes2(s,i) = targetRes(s,i);
        end
    end
    
    aprxINminGZ2 = aprxINminGZ;
    
    %%
    
    
    
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
        
        
        %%%%%%%%
        %%
        %% SOLVE FOR TARES BY TAKING THE MEAN
        for i=1:nseries
            zoop = zeros(length(excessVec0(:,1)),dimFlag);
            kx=indexLocalZero(i)-1;
            
            for m=indexLocalZero(i):indexLocalZero(i+1)-1
                zoop(m-kx,:) = aprxINminGZ2(m,:)-targetMatrix(m,:);
            end
            
            zap(i,:) = mean(zoop)*numpts/(indexLocalZero(i+1)-indexLocalZero(i));
            
        end
        %%
        %%
        for i=1:nseries
            for m=indexLocalZero(i):indexLocalZero(i+1)-1
                taretalGRBF(m,:) = zap(i,:);
            end
        end
        %%
        %%
        taresGRBF = zap;
        
        tareGRBFHist{u} = taresGRBF;
        
        
        targetRes2 = targetMatrix-aprxINminGZ2+taretalGRBF;      %0=b-Ax
        newRes2 = targetRes2'*targetRes2;
        resSquare2 = diag(newRes2);
        resSquareHist(u,:) = resSquare2;
        
    end
    
    
    % *************
    
    
    %OUTPUTS FOR GRBF SECTION
    
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
    %%    hold off
    %end
    
    %*****************
    
    
    %OUTPUT HISTOGRAM PLOTS
    if hist_FLAG == 1 && balCal_FLAG == 2;
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
    
    
    
    
    
    %%%%%%%%
    %% Number of bases
    %numBasis
    %%
    
    if print_FLAG == 1;
        
        
        %
        disp(' ***** ');
        disp('  ');
        disp('ALG+GRBF CALIBRATION MODEL GLOBAL LOAD APPROXIMATION: Check CALIB_AOX_GLOBAL_GRBF_RESULT.csv file');
        
        filename = 'CALIB_AOX_GLOBAL_GRBF_RESULT.csv';
        Z = aprxINminGZ2;
        %     xlRange = 'A1:JnumBasis';
        xlRange=char(strcat('A1:J',num2str(numBasis)));
        %xlswrite(filename,Z,xlRange)
        %%%%%%
        
        
        disp(' ');
        disp('Number of GRBFs =');
        disp(numBasis);
        disp(' ');
        
        %%%%%%
        
        twoSigmaGRBF = standardDev'.*2;
        calib_GRBF_2Sigma = array2table(twoSigmaGRBF,'VariableNames',loadlist(1:dimFlag))
        
        numSeriesNames = num2str(numSeries);
        numSeriesNameCell = cellstr(numSeriesNames);
        %Should I use strtrim()  ? -AAM 042116
        calib_GRBF_Tares = array2table(taresGRBF,'VariableNames',loadlist(1:dimFlag))
        
        calib_mean_GRBF_Resids_sqrd = array2table(resSquare2'./numpts,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Pcnt_Capacity_Max_Mag_Load_Resids = array2table(perGoop2,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Std_Dev_pcnt = array2table(stdDevPercentCapacity2,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Max_Load_Resids = array2table(maxTargets2,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Min_Load_Resids = array2table(minTargets2,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Ratio_Max_Mag_Load_Resid_and_Std_Dev = array2table(ratioGoop2,'VariableNames',loadlist(1:dimFlag))
        
        %% Prints the GRBF minmax
        calib_GRBF_minmaxband_per_capacity = array2table(theminmaxband2,'VariableNames',loadlist(1:dimFlag))
        %%
        
        %%%%%
        
        
    end
    
    
    
    %%%%%%%%
    
    if excel_FLAG == 1 && balCal_FLAG == 2;
        
    end
    
    %%%%%%%
    %        APPROX_AOX_GRBF_ws = wHist;
    %        APPROX_AOX_GRBF_coeffs  = cHist;
    %        APPROX_AOX_GRBF_Centers = centerIndexHist;
    
    filename = 'APPROX_AOX_GRBF_ws.csv';
    Z = wHist;
    %     xlRange = 'A1:JnumBasis';
    xlRange=char(strcat('A1:J',num2str(numBasis)));
    %xlswrite(filename,Z,xlRange)
    
    filename = 'APPROX_AOX_GRBF_coeffs.csv';
    Z = cHist;
    %     xlRange = 'A1:JnumBasis';
    xlRange=char(strcat('A1:J',num2str(numBasis)));
    %xlswrite(filename,Z,xlRange)
    
    filename = 'APPROX_AOX_GRBF_Centers.csv';
    Z = centerIndexHist;
    %     xlRange = 'A1:JnumBasis';
    xlRange=char(strcat('A1:J',num2str(numBasis)));
    %xlswrite(filename,Z,xlRange)
    %%%%%%%
    
    
end



%
% End GRBF Option
%

%end


%
%
% Start Validation Option
%
if balVal_FLAG == 1
    %%
    %
    %
    % Copyright ©2016 Andrew Meade and Ali Arya Mokhtarzadeh.  All Rights Reserved.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %
    
    excessVecvalid0 = excessVecvalid;
    
    % num of data points
    numptsvalid = length(seriesvalid);
    %
    
    dimFlagvalid = length(excessVecvalid(1,:));
    
    %find the average natural zeros (also called global zeros)
    globalZerosvalid = mean(natzerosvalid);
    
    %load capacities
    loadCapacitiesvalid(loadCapacitiesvalid == 0) = realmin;
    
    %find number of series; this will tell us the number of tares
    nseriesvalid = max(seriesvalid);
    
    %find zero points of each series and number of points in a series
    %localZerosAllPoints is the same as localZeroMatrix defined in the RBF
    %section - see if this can be cleaned up! is if/else necessary?  AAM042016
    prevVal = 1;
    for i=1:length(seriesvalid)
        if seriesvalid(i)==prevVal
            localZerosvalid(prevVal,:) = excessVecvalid0(i,:);
            localZerosAllPointsvalid(i,:) = localZerosvalid(prevVal,:);
            prevVal = prevVal+1;
        else
            counterval(prevVal) = i;
            localZerosAllPointsvalid(i,:) = localZerosvalid(prevVal-1,:);
        end
        globalZerosAllPointsvalid(i,:) = globalZerosvalid;
    end
    
    %add '1' to every value of counterval to find the indices of the local zeros
    indexLocalZerovalid = counterval(:)+1;
    
    %%%%%
    %    diaghistoryexcessVec0 = excessVecvalid0 - excessVec0; % AJM
    %%%%%%
    
    
    %% Subtract the Global Zeros from the Inputs and Local Zeros %%%%%%%%%%
    
    for k=1:dimFlagvalid
        %
        dainputsvalid(:,k) = excessVecvalid0(:,k)-globalZerosvalid(k);
        %
        dalzvalid(:,k) = localZerosAllPointsvalid(:,k)-globalZerosvalid(k);
        %
        dagzvalid(:,k) = 0;
        %
    end
    %%%%%%%%%%%%
    
    
    
    %find the number of points in each series and the number of series
    for j = 1:length(counterval)-1
        numPointsval(j) = counterval(j+1)-counterval(j);
        numSeriesval(j,:) = j;
    end
    
    
    
    %
    
    %%
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        %%% 5/16/18
    %Remember that  excessVec = excessVec0_complete - globalZerosAllPoints;
    excessVecvalidkeep = excessVecvalid0  - globalZerosAllPointsvalid;
    %%%
    
    %%
    %% Build the Algebraic Model
    %%
    
    
    lasttarevalid = seriesvalid(numptsvalid);
    
    %% Full Algebraic Model
    if model_FLAG == 1
        nterms = 2*dimFlag*(dimFlag+2);
    end
    
    %% Truncated Algebraic Model
    if model_FLAG == 2;
        nterms = dimFlag*(dimFlag+3)/2;
    end
    
    %% Linear Algebraic Model
    if model_FLAG == 3;
        nterms = dimFlag;
    end
    
    
    
    % Call the Algebraic Subroutine
    %
    
    comGZvalid = zeros(nterms+1,1);
    
    [comINvalid,comLZvalid,comGZvalid]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numptsvalid,seriesvalid,1, dainputsvalid,dalzvalid,dagzvalid);
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
    
    
    %%
    %% SOLVE FOR TARES BY TAKING THE MEAN
    for i=1:nseriesvalid
        zoop = zeros(length(excessVecvalidkeep(:,1)),dimFlag);
        kx=indexLocalZerovalid(i)-1;
        
        for m=indexLocalZerovalid(i):indexLocalZerovalid(i+1)-1
            zoop(m-kx,:) = checkitvalid(m,:);
        end
        
        zapvalid(i,:) = mean(zoop)*numptsvalid/(indexLocalZerovalid(i+1)-indexLocalZerovalid(i));
        %Bootstrap residuals
        zapvalid_ci(2*i-1:2*i,:)=bootci(5000,{meanResidual,zoop},'alpha',.01);
        %end bootstrap
    end
    %%
    %%
    for i=1:nseriesvalid
        for m=indexLocalZerovalid(i):indexLocalZerovalid(i+1)-1
            taretalvalid(m,:) = zapvalid(i,:);
        end
    end
    %%
    %%
    
    
    
    %RESIDUAL
    
    targetResvalid = targetMatrixvalid-aprxINminGZvalid+taretalvalid;
    
    %targetMatrixGlobalvalid = targetMatrixvalid + taretalvalid; % just for testing ajm 5_7_18
    
    % 0=b-Ax
    
    %find the sum of squares of the residual using the dot product
    
    resSquarevalid = dot(targetResvalid,targetResvalid)';
    %AAM note to self - in matlab, diag(A'*A) is the same as dot(A,A)'
    
    
    aprxINminGZvalidprime = targetMatrixvalid+taretalvalid;
    
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
    if hist_FLAG == 1 && balCal_FLAG == 2;
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
    
    
    
    %
    if print_FLAG == 1;
        %
        %% Full Algebraic Model
        if model_FLAG == 1;
            disp(' ');
            disp('%%%%%%%%%%%%%%%%%');
            disp(' ');
            disp('VALIDATION RESULTS: Full Algebraic Model');
        end
        
        %% Truncated Algebraic Model
        if model_FLAG == 2;
            disp(' ');
            disp('%%%%%%%%%%%%%%%%%');
            disp(' ');
            disp('VALIDATION RESULTS: Truncated Algebraic Model');
        end
        
        %% Linear Algebraic Model
        if model_FLAG == 3;
            disp(' ');
            disp('%%%%%%%%%%%%%%%%%');
            disp(' ');
            disp('VALIDATION RESULTS: Linear Algebraic Model');
        end
        
        
        %%%%%%%%
        
        disp('  ');
        disp('Validation data file read =');
        disp(out.savePathval);
        disp('  ');
        disp('Number of validation data points =');
        disp(numptsvalid);
        disp('  ');
        
        %%%%
        disp('ALG VALIDATION MODEL GLOBAL LOAD APPROXIMATION: VALID_AOX_GLOBAL_ALG_RESULT in Workspace');
        disp(' ');
        
        filename = 'VALID_AOX_GLOBAL_ALG_RESULT.csv';
        Z = aprxINminGZvalid;
        xlRange = 'matrixcolumnlabels(1)1:matrixcolumnlabels(dimFlag)numpts';
        %xlswrite(filename,Z,xlRange)
        %%%%
        
        
        alg_Tares_valid = array2table(zapvalid,'VariableNames',loadlist(1:dimFlag))
        
        mean_alg_Resids_sqrd_valid = array2table(resSquarevalid'./numptsvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Pcnt_Capacity_Max_Mag_Load_Resids_valid = array2table(perGoopvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Std_Dev_pcnt_valid = array2table(stdDevPercentCapacityvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Max_Load_Resids_valid = array2table(maxTargetsvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Min_Load_Resids_valid = array2table(minTargetsvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Ratio_Max_Mag_Load_Resid_and_Std_Dev_valid = array2table(ratioGoopvalid,'VariableNames',loadlist(1:dimFlag))
        
        %% Prints the minmaxband
        alg_per_minmaxband_valid = array2table(theminmaxbandvalid,'VariableNames',loadlist(1:dimFlag))
        %%
        
        %%%%
        
    end
    
    
    if res_FLAG == 1
        figure('Name','Algebraic Model Validation; Residuals of Load Versus Data Point Index','NumberTitle','off')
        plotResPages(seriesvalid, targetResvalid, loadCapacities, stdDevPercentCapacityvalid )
        %    hold off
    end
    
    
    %
    %
    %
    % Copyright ©2016 Andrew Meade and Ali Arya Mokhtarzadeh.  All Rights Reserved.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    RBF SECTION FOR VALIDATION     AJM 12/10/16                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %goal to use centers, width and coefficients to validate parameters against
    %independent data
    
    %%
    
    for k=1:dimFlagvalid
        for s=1:length(seriesvalid)
            targetResvalidX(s,k) = targetResvalid(s,k);
        end
    end
    targetRes2valid = targetResvalidX;
    aprxINminGZ2valid = aprxINminGZvalid;
    
    
    %% Subtract the Global Zeros from the Inputs and Local Zeros %%%%%%%%%%
    
    for k=1:dimFlagvalid
        %
        dainputsvalid(:,k) = excessVecvalid0(:,k)-globalZerosvalid(k);
        %
        dalzvalid(:,k) = localZerosAllPointsvalid(:,k)-globalZerosvalid(k);
        %
        dagzvalid(:,k) = 0;
        %
    end
    %%%%%%%%%%%%
    
    
    for k=1:dimFlag % ajm 6_8_18
        dainputscalib(:,k) = excessVec0(:,k) - globalZeros(k);
    end
    
    %%
    
    if balCal_FLAG == 2
        
        etaHistvalid = cell(numBasis,1);
        aprxINminGZ_Histvalid = cell(numBasis,1);
        tareHistvalid = cell(numBasis,1);
        
        for bleh=1:length(indexLocalZerovalid)-1
            for dleh = indexLocalZerovalid(bleh):indexLocalZerovalid(bleh+1)-1
                localZeroMatrixvalid(dleh,:) = localZerosvalid(bleh,:);
            end
        end
        
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
            
            
            
            
            % *************
            %% SOLVE FOR TARES BY TAKING THE MEAN
            
            for i=1:nseriesvalid
                zoop2 = zeros(length(excessVecvalid(:,1)),dimFlag);
                kx=indexLocalZerovalid(i)-1;
                
                for m=indexLocalZerovalid(i):indexLocalZerovalid(i+1)-1
                    zoop2(m-kx,:) = aprxINminGZ2valid(m,:)-targetMatrixvalid(m,:);
                end
                
                zapvalid2(i,:)=mean(zoop2)*numptsvalid/(indexLocalZerovalid(i+1)-indexLocalZerovalid(i));
                
            end
            
            %%
            %%
            for i=1:nseriesvalid
                for m=indexLocalZerovalid(i):indexLocalZerovalid(i+1)-1
                    taretalvalid2(m,:) = zapvalid2(i,:);
                end
            end
            %%
            %%
            
            taresGRBFvalid = zapvalid2;
            tareHistvalid{u} = taresGRBFvalid;
            
            
            targetRes2valid = targetMatrixvalid+taretalvalid2-aprxINminGZ2valid;      %0=b-Ax
            targetMatrixGlobalGRBFvalid = targetMatrixvalid+taretalvalid2;  % temp for ajm 6_7_18
            
            
            newRes2valid = targetRes2valid'*targetRes2valid;
            resSquare2valid = diag(newRes2valid);
            resSquareHistvalid(u,:) = resSquare2valid;
            
        end
        
        % *************
        
        
        for b=1:length(indexLocalZerovalid)-1
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
        
        
        %%%%%%
        
        
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
        if hist_FLAG == 1 && balCal_FLAG == 2;
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
        
        
        
        
        %%%%%%%%
        %% Number of bases
        %numBasis
        %%
        
        %
        if print_FLAG == 1;
            %
            disp(' ***** ');
            disp(' ');
            disp('Number of GRBFs =');
            disp(numBasis);
            
            %%%%
            disp(' ');
            disp('ALG+GRBF VALIDATION MODEL GLOBAL LOAD APPROXIMATION: Check VALID_AOX_GLOBAL_GRBF_RESULT.csv file');
            disp(' ');
            
            filename = 'VALID_AOX_GLOBAL_GRBF_RESULT.csv';
            Z = aprxINminGZ2valid;
            xlRange = 'matrixcolumnlabels(1)1:matrixcolumnlabels(dimFlag)numpts';
            %xlswrite(filename,Z,xlRange)
            %%%%
            
            
            %%%%
            twoSigmaGRBFvalid = standardDevvalid'.*2;
            GRBF_2Sigmavalid = array2table(twoSigmaGRBFvalid,'VariableNames',loadlist(1:dimFlag))
            
            numSeriesNames = num2str(numSeries);
            numSeriesNameCell = cellstr(numSeriesNames);
            %Should I use strtrim()  ? -AAM 042116
            GRBF_Taresvalid = array2table(taresGRBFvalid,'VariableNames',loadlist(1:dimFlag))
            
            
            mean_GRBF_Resids_sqrdvalid = array2table(resSquare2valid'./numptsvalid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Pcnt_Capacity_Max_Mag_Load_Resid_valid = array2table(perGoop2valid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Std_Dev_pcnt_valid = array2table(stdDevPercentCapacity2valid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Max_Load_Resids_valid = array2table(maxTargets2valid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Min_Load_Resids_valid = array2table(minTargets2valid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Ratio_Max_Mag_Load_Resid_and_Std_Dev_valid = array2table(ratioGoop2valid,'VariableNames',loadlist(1:dimFlag))
            
            %% Prints the GRBF minmax
            GRBF_minmaxband_per_capacity_valid = array2table(theminmaxband2valid,'VariableNames',loadlist(1:dimFlag))
            %%
            
            %%%%
            
            %
        end
        %
        %%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
        %%    hold off
        %end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if excel_FLAG == 1 && balCal_FLAG == 2;
            
        end
        
    end
    
    
    
    % End Validation Option
    %
end

%
%
% Start Approximation Option
if balApprox_FLAG == 1
    %%
    
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
    %     inputApprox_balCal = 'MK14C-ChkLds-Ames2011-Meade-8D_voltage.csv';
    %     excessVecapprox =         csvread(inputApprox_balCal,19,12,'M20..T143');
    load(out.savePathapp,'-mat');
    %
    
    % testing
    % nseriesapprox = nseriesvalid;
    nseriesapprox = 11;
    indexLocalZeroapprox(1) = 1;
    indexLocalZeroapprox(2) = 10;
    indexLocalZeroapprox(3) = 19;
    indexLocalZeroapprox(4) = 28;
    indexLocalZeroapprox(5) = 37;
    indexLocalZeroapprox(6) = 46;
    indexLocalZeroapprox(7) = 55;
    indexLocalZeroapprox(8) = 66;
    indexLocalZeroapprox(9) = 77 ;
    indexLocalZeroapprox(10) = 88;
    indexLocalZeroapprox(11) = 99;
    indexLocalZeroapprox(12) = 110;
    loadCapacitiesapprox(1) = 2500;
    loadCapacitiesapprox(2) = 2500;
    loadCapacitiesapprox(3) = 1250;
    loadCapacitiesapprox(4) = 1250;
    loadCapacitiesapprox(5) = 5000;
    loadCapacitiesapprox(6) = 700;
    % testing
    
    % num of data points
    numptsapprox = length(excessVecapprox);
    %
    
    
    dimFlagapprox = length(excessVecapprox(1,:));
    
    
    %find the average natural zeros (also called global zeros)
    globalZerosapprox = mean(natzerosapprox);
    
    
    %%% make an array out of the globalZerosapprox vector
    for i=1:numptsapprox
        globalZerosAllPointsapprox(i,:) = globalZerosapprox;
    end
    %%%
    
    
    
    %% Subtract the Global Zeros from the Inputs %%%%%%%%%%
    
    for k=1:dimFlagapprox
        dainputsapprox(:,k) = excessVecapprox(:,k)-globalZerosAllPointsapprox(:,k);
        dalzapprox(:,k) = globalZerosAllPointsapprox(:,k)-globalZerosAllPointsapprox(:,k);
    end
    
    %%%%%%%%%%%%
    
    
    
    
    
    %%
    %% Build the Algebraic Model
    %%
    
    %lasttareapprox = seriesapprox(numptsapprox);
    
    %% Full Algebraic Model
    if model_FLAG == 1
        nterms = 2*dimFlag*(dimFlag+2);
    end
    
    %% Truncated Algebraic Model
    if model_FLAG == 2;
        nterms = dimFlag*(dimFlag+3)/2;
    end
    
    %% Linear Algebraic Model
    if model_FLAG == 3;
        nterms = dimFlag;
    end
    
    
    
    % Call the Algebraic Subroutine
    %
    comGZapprox= zeros(nterms+1,1);
    
    for i=1:dimFlag
        biggeeapprox(:,i) = 0;
    end
    
    
    [comINapprox,comLZapprox,comGZapprox]=balCal_algEquations3(model_FLAG,nterms,dimFlag,numptsapprox,0,1,dainputsapprox,dalzapprox,biggeeapprox);
    
    %%
    %%
    
    %    xapprox = xcalib;
    
    for i=1:nterms+1
        xapprox(i,:) = xcalib(i,:);
    end
    
    
    %LOAD APPROXIMATION
    %define the approximation for inputs minus global zeros
    interceptsapprox = -(comGZapprox'*xapprox);
    aprxINapprox = (xapprox'*comINapprox)';        %to find ?? AJM111516
    %%
    %%
    
    for m=1:length(aprxINapprox)
        %%%% 3/23/18 Remove Intercepts %%%%%
        %        aprxINminGZapprox(m,:) = aprxINapprox(m,:)+interceptsapprox;
        aprxINminGZapprox(m,:) = aprxINapprox(m,:);
        %%%%%%%%%%%%%%%%%%
        
        stdevchecktestapprox(m,:) = aprxINminGZapprox(m,:);
        
    end
    %%
    %%
    
    
    %%
    %% SOLVE FOR TARES
    for i=1:nseriesapprox
        zoopapprox = zeros(length(excessVecapprox(:,1)),dimFlag);
        
        kx=indexLocalZeroapprox(i)-1;
        
        daemmlengthapprox = indexLocalZeroapprox(i+1) - indexLocalZeroapprox(i);
        
        stdevchecktestapprox = zeros(daemmlengthapprox,dimFlag);   %% ajm 7_18_18
        
        for m= indexLocalZeroapprox(i): indexLocalZeroapprox(i+1)-1
            zoopapprox(m-kx,:) = aprxINminGZapprox(m,:);
            stdevchecktestapprox(m-kx,:) = aprxINminGZapprox(m,:);
        end
        
        zapapprox(i,:) = mean(zoopapprox)*numptsapprox/(indexLocalZeroapprox(i+1)-indexLocalZeroapprox(i));
        
        zapstdevapprox(i,:) =  std(stdevchecktestapprox);  %%% ajm 7_17_18
        
    end
    %%
    %%
    
    for i=1:nseriesapprox
        
        for j= 1: dimFlag
            stdevfilterapprox(i,j) = 100.0*zapstdevapprox(i,j)/loadCapacitiesapprox(1,j); %% ajm 7_17_18
            
            
            if stdevfilterapprox(i,j) > 0.25
                zapapprox(i,j) = aprxINminGZapprox(indexLocalZeroapprox(i),j);
            end
            
        end
        
        
        for m=indexLocalZeroapprox(i):indexLocalZeroapprox(i+1)-1
            taretalapprox(m,:) = zapapprox(i,:);
        end
        
    end
    
    %%
    %%
    
    
    
    disp(' ');
    disp('%%%%%%%%%%%%%%%%%');
    disp('  ');
    disp('Approximation data file read =');
    disp(out.savePathapp);
    disp('  ');
    
    
    %%%%%%
    disp('  ');
    disp('ALG MODEL GLOBAL LOAD APPROXIMATION: Check APPROX_AOX_GLOBAL_ALG_RESULT.csv file');
    disp(' ');
    
    filename = 'APPROX_AOX_GLOBAL_ALG_RESULT.csv';
    Z = aprxINminGZapprox;
    xlRange = 'A1:JnumBasis';
    %xlswrite(filename,Z,xlRange)
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
            
            
            
            %UPDATE THE RESIDUAL
            
            %update the approximation
            
            aprxINminGZ2approx = aprxINminGZ2approx+rbfc_INminGZapprox;
            aprxINminGZ_Histapprox{u} = aprxINminGZ2approx;
            
            
            %%
            %% SOLVE FOR TARES BY TAKING THE MEAN
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
            %%
            %%
            
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
            
            %%
            %%
            
            
        end
        
        %%%%%%
        disp(' ');
        disp('ALG+GRBF MODEL GLOBAL LOAD APPROXIMATION: Check APPROX_AOX_GLOBAL_GRBF_RESULT.csv file');
        disp(' ');
        
        filename = 'APPROX_AOX_GLOBAL_GRBF_RESULT.csv';
        Z = aprxINminGZ2approx;
        xlRange = 'A1:JnumBasis';
        %xlswrite(filename,Z,xlRange)
        %%%%%%
        
        
        if excel_FLAG == 1 && balCal_FLAG == 2;
            
        end
        
    end
    
    
    
    
    %
    %End Approximation Option
end

disp('  ')
disp('Calculations Complete.')

% Tidy up the Workspace
%clearvars -except CALIB_AOX_GLOBAL_ALG_RESULT CALIB_AOX_GLOBAL_GRBF_RESULT VALID_AOX_GLOBAL_ALG_RESULT VALID_AOX_GLOBAL_GRBF_RESULT APPROX_AOX_GLOBAL_ALG_RESULT APPROX_AOX_GLOBAL_GRBF_RESULT APPROX_AOX_COEFF_MATRIX  APPROX_AOX_GRBF_Centers APPROX_AOX_GRBF_coeffs APPROX_AOX_GRBF_ws BALFIT_DATA_REDUCTION_MATRIX_IN_AMES_FORMAT  OUTLIER_ROWS OUTLIER_ROWS2 zeroed_series zeroed_targetMatrix zeroed_excessVec outlierIndices   % 5_11_18 ajm

%
% Copyright ©2016 Andrew Meade and Ali Arya Mokhtarzadeh.  All Rights Reserved.
%

%toc

