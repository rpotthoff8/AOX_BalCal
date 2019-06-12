% Copyright 2019 Andrew Meade, Ali Arya Mokhtarzadeh and Javier Villarreal.  All Rights Reserved.
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
% close all;
workspace;
fprintf('Copyright 2019 Andrew Meade, Ali Arya Mokhtarzadeh and Javier Villarreal.  All Rights Reserved.\n')
% The mean of the approximation residual (testmatrix minus local approximation) for each section is taken as the tare for that channel and section. The tare is subtracted from the global values to make the local loads. The accuracy of the validation, when compared to the known loads, is very sensitive to the value of the tares (which is unknown) and NOT the order of the calibration equations.
% Because of measurement noise in the voltage the APPROXIMATION tare is computed by post-processing. The average and stddev is taken of all channels per section. If the stddev is less than 0.25% of the capacity for any station the tare is equal to the average for that channel. If the stddev is greater than 0.25% then the approximation at the local zero is taken as the tare for that channel and section. The global approximation is left alone but the tare is subtracted from the values to make the local loads. Line 3133.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       USER INPUT SECTION
out = AOX_GUI;
if out.cancel == 1
    return
end
%TO SELECT Algebraic Model                             set FLAGS.balCal = 1;
%TO SELECT Algebraic and GRBF Model                    set FLAGS.balCal = 2;
FLAGS.balCal = out.grbf;
%DEFINE THE NUMBER OF BASIS FUNCTIONS
numBasis = out.basis;
%
%TO SELECT INDIRECT APPROACH                         set FLAGS.approach = 1;
FLAGS.approach = 0;%out.approach;
%
%SELECT ALGEBRAIC MODEL
%          set FLAGS.model = 1 (full), 2 (trunc), 3 (linear), or 4 (custom);
FLAGS.model = out.model;
%
%TO PRINT LOAD PERFORMANCE PARAMETERS                   set FLAGS.print = 1;
FLAGS.print = out.tables;
%
%TO SAVE DATA TO CSV                                    set FLAGS.excel = 1;
FLAGS.excel = out.excel;
%
%TO PRINT INPUT/OUTPUT CORRELATION PLOTS                 set FLAGS.corr = 1;
FLAGS.corr = out.corr;
%
%TO PRINT INPUT/RESIDUALS CORRELATION PLOTS           set FLAGS.rescorr = 1;
FLAGS.rescorr = out.rescorr;
%
%TO PRINT ORDER/RESIDUALS PLOTS                          set rest_FLAG = 1;
FLAGS.res = out.res;
%
%TO PRINT RESIDUAL HISTOGRAMS                            set FLAGS.hist = 1;
FLAGS.hist = out.hist;
%
%TO SELECT Validation of the Model                     set FLAGS.balVal = 1;
FLAGS.balVal = out.valid;
%
%TO SELECT Approximation from Cal Data              set FLAGS.balApprox = 1;
FLAGS.balApprox = out.approx;
%
%TO FLAG POTENTIAL OUTLIERS                            set FLAGS.balOut = 1;
FLAGS.balOut = out.outlier;
numSTD = out.numSTD;  %Number of standard deviations for outlier threshold.
%
%TO REMOVE POTENTIAL OUTLIERS                          set FLAGS.zeroed = 1;
FLAGS.zeroed = out.zeroed;
%
%TO USE LATIN HYPERCUBE SAMPLING set                          FLAGS.LHS = 1;
FLAGS.LHS = out.lhs;
numLHS = out.numLHS; %Number of times to iterate.
LHSp = out.LHSp; %Percent of data used to create sample.
%
%Uncertainty button outputs
numBoot=out.numBoot;
FLAGS.boot=out.bootFlag;
FLAGS.volt=out.voltFlag;
voltTrust=out.voltTrust;

FLAGS.anova = out.anova;
%                       END USER INPUT SECTION
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data and characterize series
load(out.savePathcal,'-mat');
series0 = series;
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

% Prints output vs. input and calculates correlations
if FLAGS.corr == 1
    figure('Name','Correlation plot','NumberTitle','off');
    correlationPlot(targetMatrix0, excessVec0, loadlist, voltagelist);
end

% Finds the average  of the natural zeros (called global zeros)
globalZeros = mean(natzeros);

% Subtracts global zeros from signal.
dainputs0 = excessVec0 - ones(numpts0,1)*globalZeros;

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
comIN0 = balCal_algEqns(FLAGS.model,dainputs0,series0);

if FLAGS.LHS == 0
    numLHS = 1;
end
lhs_check = zeros(length(excessVec0),1);
lhs_ind = 1:length(excessVec0(:,1));

if FLAGS.LHS == 1
    fprintf('\nNumber of LHS Iterations Selected: %i\n',numLHS)
end

fprintf('\nStarting Calculations\n')
% if FLAGS.approach == 1
%     fprintf('\n\nUsing the Indirect Approach for Calibration');
% else
%     fprintf('\n\nUsing the Direct Approach for Calibration');
% end

%%
for lhs = 1:numLHS

    % Creates an LHS sub-sample of the data.
    if FLAGS.LHS == 1
        sample = AOX_LHS(series0,excessVec0,LHSp);
        lhs_check(sample) = 1;
        lhs_ind(find(lhs_check-1)); % This line outputs which data points haven't been sampled yet
        pct_sampled = sum(lhs_check)/length(lhs_check); % This line outputs what percentage of points have been sampled
    else
        sample = [1:length(series0)]';
    end

    % Uses the sampling indices in "sample" to create the subsamples
    series = series0(sample);
    targetMatrix = targetMatrix0(sample,:);
    comIN = comIN0(sample,:);

    fprintf('\nWorking ...\n')

    %Calculate xcalib (coefficients)
    [xcalib, ANOVA]=calc_xcalib(comIN,targetMatrix,series,nterms,nseries0,dimFlag,FLAGS.model,customMatrix,FLAGS.anova);

    if FLAGS.LHS == 1
        x_all(:,:,lhs) = xcalib;
    end
end
if FLAGS.LHS == 1
    xcalib = mean(x_all,3);
    xcalib_std = std(x_all,[],3);
end

%%
% APPROXIMATION
% define the approximation for inputs minus global zeros
aprxIN = comIN0*xcalib;

% RESIDUAL
targetRes = targetMatrix0-aprxIN;

%%
% Identify Outliers After Filtering
% (Threshold approach) ajm 8/2/17
if FLAGS.balOut == 1

    %Identify outliers based on residuals
    [OUTLIER_ROWS,num_outliers,prcnt_outliers]=ID_outliers(targetRes,loadCapacities,numpts0,dimFlag,numSTD);

    % Use the reduced input and target files
    if FLAGS.zeroed == 1

        % Remove outlier rows for recalculation and all future calculations:
        numpts0 =  numpts0 - num_outliers;
        targetMatrix0(OUTLIER_ROWS,:) = [];
        excessVec0(OUTLIER_ROWS,:) = [];
        series0(OUTLIER_ROWS) = [];
        comIN0(OUTLIER_ROWS,:) = [];
        [~,s_1st0,~] = unique(series0);
        nseries0 = length(s_1st0);

        %Calculate xcalib (coefficients)
        [xcalib,ANOVA]=calc_xcalib(comIN0,targetMatrix0,series0,nterms,nseries0,dimFlag,FLAGS.model,customMatrix,FLAGS.anova);

        % APPROXIMATION
        % define the approximation for inputs minus global zeros
        aprxIN = comIN0*xcalib;

        % RESIDUAL
        targetRes = targetMatrix0-aprxIN;

    end
end

% Splits xcalib into Coefficients and Intercepts (which are negative Tares)
coeff = xcalib(1:nterms,:);
tares = -xcalib(nterms+1:end,:);
intercepts=-tares;
taretal=tares(series,:);
aprxINminGZ=aprxIN+taretal; %QUESTION: 29 MAR 2019: JRP

%Start uncertainty section
if FLAGS.boot==1
    %%start bootstrapfunction
    bootalpha=.05;
    f=@calc_xcalib;

    xcalib_ci=bootci(numBoot,{f,comIN0,targetMatrix0,series0,nterms,nseries0,dimFlag,FLAGS.model,customMatrix,0});
    
    %START EXTRA for testing uncert prop
    f2=@(input,input2,input3) calc_xcalib(input,input2,input3,nterms,nseries0,dimFlag,FLAGS.model,customMatrix,0);
    g=@(input,input2,input3) calc_yhat(input,input2,input3,nterms,nseries0,dimFlag,FLAGS.model,customMatrix,0,comIN0);
%     xcalib_ci2=bootci(numBoot,{f2,comIN0,targetMatrix0,series0},'type','norm');
    boot_yhat_ci=bootci(numBoot,{g,comIN0,targetMatrix0,series0}); %CHANGE BACK
%     boot_yhat_ci2=bootci(numBoot,{g,comIN0,targetMatrix0,series0},'type','norm');
    
    for i=1:size(aprxIN,1)
    for j=1:size(aprxIN,2)
        for k=1:2
            error(k,i,j)=abs(boot_yhat_ci(k,i,j)-aprxIN(i,j));
%             error2(k,i,j)=abs(boot_yhat_ci2(k,i,j)-aprxIN(i,j));
        end
        yhat_error(i,j)=max(error(:,i,j));
%         yhat_error2(i,j)=mean(error(:,i,j));
    end
    end
    %END extra for testing uncert prop
else
    xcalib_ci=zeros(2, size(xcalib,1),size(xcalib,2));
end
% END: bootstrap section

%ANOVA data for uncertainty
beta_CI_comb=zeros(size(xcalib,1),dimFlag);
y_hat_PI_comb=zeros(size(targetMatrix,1),size(targetMatrix,2));
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

[combined_uncert,tare_uncert, FL_uncert,xcalibCI_includeZero, xcalib_error,coeff_uncert_boot]=uncert_prop(xcalib,xcalib_ci,comIN0,dimFlag,uncert_comIN,s_1st0,nterms,targetMatrix0,series0,voltTrust,FLAGS.boot,FLAGS.volt);
[combined_uncert_anova,tare_uncert_anova, FL_uncert_anova,coeff_uncert_anova]=uncert_prop_anova(xcalib,beta_CI_comb,comIN,dimFlag,uncert_comIN,s_1st0,nterms,targetMatrix,series,voltTrust,FLAGS.anova,FLAGS.volt);
%end uncertainty section


% Temp for Tares AJM 4_20_19
filename = 'Tares_Double_Precision.csv';
dlmwrite(filename,tares,'precision','%.16f');
%%%%%

%  Creates Matrix for the volts to loads
%APPROX_AOX_COEFF_MATRIX = coeff;  AJM 4_19_19

xapprox = xcalib(1:nterms,:);

%xapprox = coeff;
if FLAGS.excel == 1
    filename = 'APPROX_AOX_COEFF_MATRIX.csv';
    dlmwrite(filename, [coeff;zeros(1,dimFlag)] ,'precision','%.16f'); % AJM 4_18_19
end
% Prints residual vs. input and calculates correlations
if FLAGS.rescorr == 1
    figure('Name','Residual correlation plot','NumberTitle','off');
    correlationPlot(excessVec0, targetRes, voltagelist, reslist);
end

% Calculates the Sum of Squares of the residual
resSquare = sum(targetRes.^2);

%find the sum of squares of the residual using the dot product
% resSquare = dot(targetRes,targetRes)';
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

%OUTPUT HISTOGRAM PLOTS
if FLAGS.hist == 1
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
if FLAGS.print == 1
    %
    %% END Direct or Indirect Approach to Calibration

    %% Identify the Possible Outliers
    if FLAGS.balOut == 1
        fprintf(' \n***** \n');
        fprintf(' ');
        fprintf('\nNumber of Outliers =');
        fprintf(string(num_outliers));
        fprintf('\nOutliers % of Data =');
        fprintf(string(prcnt_outliers));
    end

    % Recalculated Calibration with Reduced Matrices
    if FLAGS.zeroed == 1
        fprintf('\n ************************************************************************ \n');
        fprintf('\nFind the reduced data in zeroed_targetMatrix and zeroed_excessVec\n');
        fprintf('\n ************************************************************************ \n');
    end

    %%%%%%% 6_14_18 ajm
    calib_twoSigmaALGB = standardDev'.*2;
    calib_algebraic_2Sigma = array2table(calib_twoSigmaALGB,'VariableNames',loadlist(1:dimFlag))

    %Should I use strtrim()  ? -AAM 042116
    series_table = table([1:nseries0]','VariableNames',{'SERIES'});
    calib_algebraic_Tares = array2table(tares,'VariableNames',loadlist(1:dimFlag));
    calib_algebraic_Tares = [series_table, calib_algebraic_Tares]

    coefficientsALGB = [xcalib;intercepts];
    for coeffCount=1:size(coefficientsALGB(:,1),1)
        numCombin(coeffCount) = coeffCount;
    end
    numCombinName = cellstr(num2str(numCombin'));
    numCombinName(nterms+nseries0+1) = cellstr('Intercept');

    calib_mean_algebraic_Resids_sqrd = array2table((resSquare'./numpts0)','VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Pcnt_Capacity_Max_Mag_Load_Resids = array2table(perGoop,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Std_Dev_pcnt = array2table(stdDevPercentCapacity,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Max_Load_Resids = array2table(maxTargets,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Min_Load_Resids = array2table(minTargets,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Ratio_Max_Mag_Load_Resid_and_Std_Dev = array2table(ratioGoop,'VariableNames',loadlist(1:dimFlag))

    % Prints the minmaxband
    calib_alg_per_minmaxband = array2table(theminmaxband,'VariableNames',loadlist(1:dimFlag))
    %%%%%%%%%

end

if FLAGS.excel == 1
    % Output results to an excel file
    fprintf('\n  ');
    fprintf('\nALG CALIBRATION MODEL LOAD APPROXIMATION FILE: CALIB_AOX_ALG_RESULT.csv\n'); % AJM 4_19_19
    fprintf('\n ');

    filename = 'CALIB_AOX_ALG_RESULT.csv';
    dlmwrite(filename,aprxIN,'precision','%.16f');
end

if FLAGS.res == 1
    figure('Name','Algebraic Model Calibration; Residuals of Load Versus Data Point Index','NumberTitle','off')
    plotResPages(series0, targetRes, loadCapacities, stdDevPercentCapacity, loadlist)
    %    hold off
end

if FLAGS.balCal == 2
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
    localZerosAllPoints=tares(series0,:);
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

            coeffRBF(s) = dot(rbfINminGZ(:,s),targetRes2(:,s)) / dot(rbfINminGZ(:,s),rbfINminGZ(:,s));

            rbfc_INminLZ(:,s) = coeffRBF(s)*rbfINminLZ(:,s);
            rbfc_INminGZ(:,s) = coeffRBF(s)*rbfINminGZ(:,s);
            rbfc_LZminGZ(:,s) = coeffRBF(s)*rbfLZminGZ(:,s); %to find tares AAM042016
        end

        wHist(u,:) = w;
        cHist(u,:) = coeffRBF;
        centerIndexHist(u,:) = centerIndexLoop;
        etaHist{u} = eta;

        %update the approximation
        aprxINminGZ2 = aprxINminGZ2+rbfc_INminGZ;
        aprxINminGZ_Hist{u} = aprxINminGZ2;

        % SOLVE FOR TARES BY TAKING THE MEAN
        [taresAllPointsGRBF,taretalGRBFSTDDEV] = meantare(series0,aprxINminGZ2-targetMatrix0);

        taresGRBF = taresAllPointsGRBF(s_1st0,:);
        taresGRBFSTDEV = taretalGRBFSTDDEV(s_1st0,:);
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
    if FLAGS.res == 1
        figure('Name','GRBF + Algebraic Model Calibration; Residuals of Load Versus Data Point Index','NumberTitle','off')
        plotResPages(series0, targetRes2, loadCapacities, stdDevPercentCapacity2, loadlist)
        %    hold off
    end

    %%% Diagnostics %%%%% 5/16/18
    %
    %justtherbfs = aprxINminGZ2-aprxINminGZ;
    %
    %if FLAGS.res == 1
    %    figure('Name','Looking at GRBF Distribution in Calibration','NumberTitle','off')
    %    plotResPages(series0, justtherbfs, loadCapacities, stdDevPercentCapacity2)
    %    hold off
    %end

    %OUTPUT HISTOGRAM PLOTS
    if FLAGS.hist == 1
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

    if FLAGS.excel == 1
        fprintf('\n ***** \n');
        fprintf('\n  ');
        fprintf('\nALG+GRBF CALIBRATION MODEL LOAD APPROXIMATION: CALIB_AOX_GRBF_RESULT.csv\n'); %AJM 4_19_19

        filename = 'CALIB_AOX_GRBF_RESULT.csv';
        dlmwrite(filename,aprxINminGZ2,'precision','%.16f'); % Note that aprxINminGZ2 is actually aprxIN2

        filename = 'APPROX_AOX_GRBF_ws.csv';
        dlmwrite(filename,wHist,'precision','%.16f');

        filename = 'APPROX_AOX_GRBF_coeffs.csv';
        dlmwrite(filename,cHist,'precision','%.16f');

        filename = 'APPROX_AOX_GRBF_Centers.csv';
        dlmwrite(filename,centerIndexHist,'precision','%.16f');
    end


    if FLAGS.print == 1

        fprintf('\n ');
        %        fprintf('Number of GRBFs =');
        %        fprintf(string(numBasis));
        fprintf('\nNumber of GRBFs: %i\n',numBasis);
        fprintf('\n ');

        twoSigmaGRBF = standardDev'.*2;
        calib_GRBF_2Sigma = array2table(twoSigmaGRBF,'VariableNames',loadlist(1:dimFlag))

        %Should I use strtrim()  ? -AAM 042116
        taresGRBFactual = taresGRBF;
        series_table = table([1:nseries0]','VariableNames',{'SERIES'});
        calib_GRBF_Tares = array2table(taresGRBFactual,'VariableNames',loadlist(1:dimFlag));
        calib_GRBF_Tares = [series_table, calib_GRBF_Tares]

        calib_mean_GRBF_Resids_sqrd = array2table(resSquare2'./numpts0,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Pcnt_Capacity_Max_Mag_Load_Resids = array2table(perGoop2,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Std_Dev_pcnt = array2table(stdDevPercentCapacity2,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Max_Load_Resids = array2table(maxTargets2,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Min_Load_Resids = array2table(minTargets2,'VariableNames',loadlist(1:dimFlag))
        calib_GRBF_Ratio_Max_Mag_Load_Resid_and_Std_Dev = array2table(ratioGoop2,'VariableNames',loadlist(1:dimFlag))

        % Prints the GRBF minmax
        calib_GRBF_minmaxband_per_capacity = array2table(theminmaxband2,'VariableNames',loadlist(1:dimFlag))
    end

end

if FLAGS.balVal == 1
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
    [validSeries,s_1stV,~] = unique(seriesvalid);
    xvalid=[coeff;tares(validSeries,:)]; %QUESTION: 29 March 2019; JRP

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
    [localZerosvalid,localZerosAllPointsvalid] = localzeros(seriesvalid,excessVecvalid);
    globalZerosAllPointsvalid = ones(numptsvalid,1)*globalZerosvalid;

    % Subtract the Global Zeros from the Inputs and Local Zeros
    for k=1:dimFlagvalid
        dainputsvalid(:,k) = excessVecvalid(:,k)-globalZerosvalid(k);
        dalzvalid(:,k) = localZerosAllPointsvalid(:,k)-globalZerosvalid(k);
        dagzvalid(:,k) = 0;
    end

    %%% 5/16/18
    %Remember that  excessVec0 = excessVec0_complete - globalZerosAllPoints;
    excessVecvalidkeep = excessVecvalid  - globalZerosAllPointsvalid;
    %%%

    % Build the Algebraic Model
    lasttarevalid = seriesvalid(numptsvalid);
    % Full Algebraic Model
    if FLAGS.model == 1
        nterms = 2*dimFlag*(dimFlag+2);
    end
    % Truncated Algebraic Model
    if FLAGS.model == 2
        nterms = dimFlag*(dimFlag+3)/2;
    end
    % Linear Algebraic Model
    if FLAGS.model == 3
        nterms = dimFlag;
    end

    % Call the Algebraic Subroutine
    comGZvalid = zeros(nterms+1,1);
    [comINvalid,comLZvalid,comGZvalid]=balCal_algEquations3(FLAGS.model,nterms,dimFlag,numptsvalid,seriesvalid,nseriesvalid,dainputsvalid,dalzvalid,dagzvalid);

    comINminLZvalid = comINvalid-comLZvalid;

    %VALIDATION APPROXIMATION
    %define the approximation for inputs minus global zeros
    %    interceptsvalid = -(comGZvalid'*xvalid);  % ajm 5/17/18

    aprxINvalid = comINvalid'*xvalid;        %to find approximation AJM111516
    aprxLZvalid = comLZvalid'*xvalid;       %to find tares AAM042016

    aprxINminLZvalid = comINminLZvalid'*xvalid;

    for m=1:length(aprxINvalid(:,1))
        %%%%% 3/23/17 Zap intercepts %%%
        %        aprxINminGZvalid(m,:) = aprxINvalid(m,:)+interceptsvalid;
        aprxINminGZvalid(m,:) = aprxINvalid(m,:);
        %%%%%%

        checkitvalid(m,:) = aprxINminGZvalid(m,:)-targetMatrixvalid(m,:);
    end

    % SOLVE FOR TARES BY TAKING THE MEAN
    [taresAllPointsvalid,taretalstdvalid] = meantare(seriesvalid,checkitvalid);
    zapvalid     = taresAllPointsvalid(s_1stV,:);
    zapSTDEVvalid = taretalstdvalid(s_1stV,:);

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
    if FLAGS.hist == 1
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

    if FLAGS.print == 1
        %
        % Full Algebraic Model
        if FLAGS.model == 1
            fprintf('\n ');
            fprintf('\n%%%%%%%%%%%%%%%%%\n');
            fprintf('\n ');
            fprintf('\nVALIDATION RESULTS: Full Algebraic Model\n');
        end
        % Truncated Algebraic Model
        if FLAGS.model == 2
            fprintf('\n ');
            fprintf('\n%%%%%%%%%%%%%%%%%\n');
            fprintf('\n ');
            fprintf('\nVALIDATION RESULTS: Truncated Algebraic Model\n');
        end
        % Linear Algebraic Model
        if FLAGS.model == 3
            fprintf('\n ');
            fprintf('\n%%%%%%%%%%%%%%%%%\n');
            fprintf('\n ');
            fprintf('\nVALIDATION RESULTS: Linear Algebraic Model\n');
        end

        fprintf('\n  ');
        fprintf('\nValidation data file read: ');
        fprintf(fileNamevalid);
        fprintf('  ');
        fprintf('\nNumber of validation data points: %i\n',numptsvalid);
        fprintf('\n  ');

        series_table_valid = table([1:nseriesvalid]','VariableNames',{'SERIES'});
        alg_Tares_valid = array2table(zapvalid,'VariableNames',loadlist(1:dimFlag));
        alg_Tares_valid = [series_table_valid, alg_Tares_valid]

        mean_alg_Resids_sqrd_valid = array2table(resSquarevalid'./numptsvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Pcnt_Capacity_Max_Mag_Load_Resids_valid = array2table(perGoopvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Std_Dev_pcnt_valid = array2table(stdDevPercentCapacityvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Max_Load_Resids_valid = array2table(maxTargetsvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Min_Load_Resids_valid = array2table(minTargetsvalid,'VariableNames',loadlist(1:dimFlag))
        alg_Ratio_Max_Mag_Load_Resid_and_Std_Dev_valid = array2table(ratioGoopvalid,'VariableNames',loadlist(1:dimFlag))

        % Prints the minmaxband
        alg_per_minmaxband_valid = array2table(theminmaxbandvalid,'VariableNames',loadlist(1:dimFlag))
    end

    if FLAGS.excel == 1
        %%%%
        fprintf('\nALG VALIDATION MODEL GLOBAL LOAD APPROXIMATION: VALID_AOX_GLOBAL_ALG_RESULT in Workspace\n');
        fprintf('\n ');

        filename = 'VALID_AOX_GLOBAL_ALG_RESULT.csv';
        %        csvwrite(filename,aprxINminGZvalid)
        dlmwrite(filename,aprxINminGZvalid,'precision','%.16f');
    end

    if FLAGS.res == 1
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
        dainputsvalid(:,k) = excessVecvalid(:,k)-globalZerosvalid(k);
        dalzvalid(:,k) = localZerosAllPointsvalid(:,k)-globalZerosvalid(k);
        dagzvalid(:,k) = 0;
    end

    for k=1:dimFlag % ajm 6_8_18
        dainputscalib(:,k) = excessVec0(:,k) - globalZeros(k);
    end

    if FLAGS.balCal == 2
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
            [taresAllPointsvalid2,taretalstdvalid2] = meantare(seriesvalid,aprxINminGZ2valid-targetMatrixvalid);

            taresGRBFvalid = taresAllPointsvalid2(s_1st,:);
            taresGRBFSTDEVvalid = taretalstdvalid2(s_1st,:);
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
        if FLAGS.hist == 1
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
        if FLAGS.print == 1
            %
            fprintf('\n ***** \n');
            fprintf('\n ');
            fprintf('\nNumber of GRBFs: %i\n',numBasis);

            twoSigmaGRBFvalid = standardDevvalid'.*2;
            GRBF_2Sigmavalid = array2table(twoSigmaGRBFvalid,'VariableNames',loadlist(1:dimFlag))

            %Should I use strtrim()  ? -AAM 042116
            series_table_valid = table([1:nseriesvalid]','VariableNames',{'SERIES'});
            GRBF_Taresvalid = array2table(taresGRBFvalid,'VariableNames',loadlist(1:dimFlag));
            GRBF_Taresvalid = [series_table_valid, GRBF_Taresvalid]

            mean_GRBF_Resids_sqrdvalid = array2table(resSquare2valid'./numptsvalid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Pcnt_Capacity_Max_Mag_Load_Resid_valid = array2table(perGoop2valid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Std_Dev_pcnt_valid = array2table(stdDevPercentCapacity2valid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Max_Load_Resids_valid = array2table(maxTargets2valid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Min_Load_Resids_valid = array2table(minTargets2valid,'VariableNames',loadlist(1:dimFlag))
            GRBF_Ratio_Max_Mag_Load_Resid_and_Std_Dev_valid = array2table(ratioGoop2valid,'VariableNames',loadlist(1:dimFlag))

            % Prints the GRBF minmax
            GRBF_minmaxband_per_capacity_valid = array2table(theminmaxband2valid,'VariableNames',loadlist(1:dimFlag))
        end

        if FLAGS.res == 1
            figure('Name','GRBF + Algebraic Model Validation; Residuals of Load Versus Data Point Index','NumberTitle','off')
            plotResPages(seriesvalid, targetRes2valid, loadCapacities, stdDevPercentCapacity2valid, loadlist)
            %    hold off
        end

        % Diagnostics %%%%% 5/16/18
        %justtherbfsvalid = aprxINminGZ2valid-aprxINminGZvalid;
        %
        %%deltatherbfsvalid = justtherbfs - justtherbfsvalid;
        %
        %if FLAGS.res == 1
        %    figure('Name','Looking at GRBF Distribution in Validation','NumberTitle','off')
        %    plotResPages(seriesvalid, justtherbfsvalid, loadCapacities, stdDevPercentCapacity2valid)
        %    hold off
        %end

        if FLAGS.excel == 1
            fprintf('\n ');
            fprintf('\nALG+GRBF VALIDATION MODEL GLOBAL LOAD APPROXIMATION:  VALID_AOX_GLOBAL_GRBF_RESULT.csv\n');
            fprintf('\n ');

            filename = 'VALID_AOX_GLOBAL_GRBF_RESULT.csv';
            csvwrite(filename,aprxINminGZ2valid)
            dlmwrite(filename,aprxINminGZ2valid,'precision','%.16f');

        end
    end
end

if FLAGS.balApprox == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                        APPROXIMATION SECTION      AJM 6/29/17           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %DEFINE THE PRODUCTION CSV INPUT FILE AND SELECT THE RANGE OF DATA VALUES TO READ

    load(out.savePathapp,'-mat');

    % num of data points
    numptsapprox = length(excessVecapprox);

    %natural zeros (also called global zeros)
    globalZerosapprox = mean(natzerosapprox);

    % make an array out of the globalZerosapprox vector
    for i=1:numptsapprox
        globalZerosAllPointsapprox(i,:) = globalZerosapprox;
    end

    % Subtract the Global Zeros from the Inputs
    for k=1:dimFlag

        dainputsapprox(:,k) = excessVecapprox(:,k)-globalZerosAllPointsapprox(:,k);

        dalzapprox(:,k) = globalZerosAllPointsapprox(:,k)-globalZerosAllPointsapprox(:,k);

    end

    %% Build the Algebraic Model

    % Full Algebraic Model
    if FLAGS.model == 1
        nterms = 2*dimFlag*(dimFlag+2);
    end

    % Truncated Algebraic Model
    if FLAGS.model == 2
        nterms = dimFlag*(dimFlag+3)/2;
    end

    % Linear Algebraic Model
    if FLAGS.model == 3
        nterms = dimFlag;
    end

    % Call the Algebraic Subroutine
    comGZapprox= zeros(nterms,1);


    for i=1:dimFlag
        biggee(:,i) = 0;
    end

    [comINapprox,comLZapprox,comGZapprox]=balCal_algEquations3(FLAGS.model,nterms,dimFlag,numptsapprox,0,0,dainputsapprox,dalzapprox,biggee);

    %LOAD APPROXIMATION
    %define the approximation for inputs minus global zeros
    interceptsapprox = -(comGZapprox'*coeff);
    aprxINapprox = ( coeff'*comINapprox)';        %to find ?? AJM111516

    for m=1:length(aprxINapprox)
        aprxINminGZapprox(m,:) = aprxINapprox(m,:);
    end

    if FLAGS.excel == 1
        fprintf('\n ');
        fprintf('\n%%%%%%%%%%%%%%%%%\n');
        fprintf('\n ');
        fprintf('\nALG MODEL APPROXIMATION RESULTS: GLOBAL_ALG_APPROX.csv\n');

        filename = 'GLOBAL_ALG_APPROX.csv';
        csvwrite(filename,aprxINminGZapprox)
        dlmwrite(filename,aprxINminGZapprox,'precision','%.16f');
    else

        fprintf('\n ');
        fprintf('\n%%%%%%%%%%%%%%%%%\n');
        fprintf('\n ');
        fprintf('\nALG MODEL APPROXIMATION RESULTS: Check aprxINminGZapprox in Workspace\n');
    end


    %%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    RBF SECTION FOR APPROXIMATION     AJM 6/29/17                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %goal to use centers, width and coefficients to approxate parameters against
    %independent data

    aprxINminGZ2approx = aprxINminGZapprox;

    if FLAGS.balCal == 2

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

        if FLAGS.excel == 1
            fprintf('\n ');
            fprintf('\n%%%%%%%%%%%%%%%%%\n');
            fprintf('\n ');

            fprintf('\nALG + GRBF MODEL APPROXIMATION RESULTS: GLOBAL_ALG+GRBF_APPROX.csv\n');

            filename = 'GLOBAL_ALG+GRBF_APPROX.csv';

            csvwrite(filename,aprxINminGZ2approx)
            dlmwrite(filename,aprxINminGZ2approx,'precision','%.16f');
        else
            fprintf('\n ');
            fprintf('\nGRBF MODEL APPROXIMATION RESULTS: Check aprxINminGZ2approx in Workspace\n');
            fprintf('\n ');
        end

    end

    %End Approximation Option
end

fprintf('\n  ');
fprintf('\nCalculations Complete.\n');
fprintf('\n');
