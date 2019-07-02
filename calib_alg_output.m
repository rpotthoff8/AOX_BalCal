%Function creates all the outputs for the calibration, algebraic section
%This simplifies following the main code

function [] = calib_alg_output(FLAGS,coeff,dimFlag,excessVec0,targetRes,voltagelist,reslist,loadCapacities,nterms,nseries0,numpts0,ANOVA,loadlist,balfitcomIN,balfitxcalib,balfittargetMatrix,balfitANOVA,tares,tares_STDDEV,xcalib,intercepts,aprxIN,series0)
% Temp for Tares AJM 4_20_19
%filename = 'Tares_Double_Precision.csv';
%dlmwrite(filename,tares,'precision','%.16f');
%%%%%

%  Creates Matrix for the volts to loads
%APPROX_AOX_COEFF_MATRIX = coeff;  AJM 4_19_19

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

%OUTPUTS FOR ALGEBRAIC SECTION %SAME START
for k=1:length(targetRes(1,:))
    [goop(k),kstar(k)] = max(abs(targetRes(:,k)));
    goopVal(k) = abs(targetRes(kstar(k),k));
    xCent(k) = excessVec0(kstar(k),k);
    maxTargets(k) = max(targetRes(:,k));
    minTargets(k) = min(targetRes(:,k));
    tR2(k) = targetRes(:,k)'*targetRes(:,k);     % AJM 6_12_19
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
%SAME END

%%% ANOVA Stats AJM 6_12_19

if FLAGS.model ~= 4 && FLAGS.anova==1
    
    totalnum = nterms+nseries0;
    totalnumcoeffs = [1:totalnum];
    totalnumcoeffs2 = [2:totalnum+1];
    dsof = numpts0-nterms-1;
    
    loadstatlist = {'Load', 'Sum_Sqrs', 'PRESS_Stat', 'DOF', 'Mean_Sqrs', 'F_Value', 'P_Value', 'R_sq', 'Adj_R_sq', 'PRESS_R_sq'};
    
    regresslist = {'Term', 'Coeff_Value', 'CI_95cnt', 'T_Stat', 'P_Value', 'VIF_A', 'Signif'};
    
    
    for k=1:dimFlag
        
        RECOMM_ALG_EQN(:,k) = [1.0*ANOVA(k).sig([1:nterms])];
        
        manoa2(k,:) = [loadlist(k), tR2(1,k), ANOVA(k).PRESS, dsof, gee(1,k), ANOVA(k).F, ANOVA(k).p_F, ANOVA(k).R_sq, ANOVA(k).R_sq_adj, ANOVA(k).R_sq_p];
        
        ANOVA01(:,:) = [totalnumcoeffs; ANOVA(k).beta'; ANOVA(k).beta_CI'; ANOVA(k).T'; ANOVA(k).p_T'; ANOVA(k).VIF'; 1.0*ANOVA(k).sig']';
        
        ANOVA1_2(:,:) = [ANOVA01([1:nterms],:)];
        
        STAT_LOAD = array2table(manoa2(k,:),'VariableNames',loadstatlist(1:10));
        
        REGRESS_COEFFS = array2table(ANOVA1_2(:,:),'VariableNames',regresslist(1:7));
        
        filename = 'DIRECT_ANOVA_STATS.xlsx';
        writetable(STAT_LOAD,filename,'Sheet',k,'Range','A1');
        writetable(REGRESS_COEFFS,filename,'Sheet',k,'Range','A4');
        
    end
    
    filename = 'DIRECT_RECOMM_CustomEquationMatrix.csv';
    dlmwrite(filename,RECOMM_ALG_EQN,'precision','%.8f');
    
end

%%% ANOVA Stats AJM 6_8_19



%%% Balfit Stats and Regression Coeff Matrix AJM 5_31_19


balfitaprxIN = balfitcomIN*balfitxcalib;
balfittargetRes = balfittargetMatrix-balfitaprxIN;

for k=1:length(balfittargetRes(1,:))
    [balfitgoop(k),balfitkstar(k)] = max(abs(balfittargetRes(:,k)));
    balfitgoopVal(k) = abs(balfittargetRes(kstar(k),k));
    balfittR2(k) = balfittargetRes(:,k)'*balfittargetRes(:,k);     % AJM 6_12_19
end

balfitdavariance = var(balfittargetRes);
balfitgee = mean(balfittargetRes);
balfitstandardDev10 = std(balfittargetRes);
balfitstandardDev = balfitstandardDev10';


voltagestatlist = {'Voltage', 'Sum_Sqrs', 'PRESS_Stat', 'DOF', 'Mean_Sqrs', 'F_Value', 'P_Value', 'R_sq', 'Adj_R_sq', 'PRESS_R_sq'};

balfitregresslist = {'Term', 'Coeff_Value', 'CI_95cnt', 'T_Stat', 'P_Value', 'VIF_A', 'Signif'};

%balfitinterceptlist = ['Intercept', '0', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A'];
balfitinterceptlist = [1, 0, 0, 0, 0, 0, 0];


if FLAGS.model ~= 4 && FLAGS.anova==1
    
    for k=1:dimFlag
        
        BALFIT_RECOMM_ALG_EQN(:,k) = 1.0*balfitANOVA(k).sig;
        
        balfitANOVA01(:,:) = [totalnumcoeffs2; balfitANOVA(k).beta'; ANOVA(k).beta_CI'; balfitANOVA(k).T'; balfitANOVA(k).p_T'; balfitANOVA(k).VIF'; 1.0*balfitANOVA(k).sig']';
        
        balfitANOVA_intercept1(1,:) = balfitinterceptlist(1,:);
        
        balfitANOVA1([1:nterms+1],:) = [balfitANOVA_intercept1(1,:); balfitANOVA01([1:nterms],:)];
        
        toplayer2(k,:) = [voltagelist(k), balfittR2(1,k), balfitANOVA(k).PRESS, dsof, balfitgee(1,k), balfitANOVA(k).F, balfitANOVA(k).p_F, balfitANOVA(k).R_sq, balfitANOVA(k).R_sq_adj, balfitANOVA(k).R_sq_p];
        
        BALFIT_STAT_VOLTAGE_1 = array2table(toplayer2(k,:),'VariableNames',voltagestatlist(1:10));
        
        BALFIT_REGRESS_COEFFS_1 = array2table(balfitANOVA1([1:nterms],:),'VariableNames',balfitregresslist(1:7));
        
        filename = 'BALFIT_ANOVA_STATS.xlsx';
        writetable(BALFIT_STAT_VOLTAGE_1,filename,'Sheet',k,'Range','A1');
        writetable(BALFIT_REGRESS_COEFFS_1,filename,'Sheet',k,'Range','A4');
        
        
    end
    
    
    %filename = 'BALFIT_RECOMM_CustomEquationMatrixTemplate.csv';
    %dlmwrite(filename,BALFIT_RECOMM_ALG_EQN,'precision','%.8f');
    
end


%%% Balfit Stats and Matrix AJM 5_31_19

%SAME START
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
%END SAME

%START PRINT OUT PERFORMANCE INFORMATION TO THE SCREEN
if FLAGS.print == 1

    % Recalculated Calibration with Reduced Matrices
    if FLAGS.zeroed == 1
        fprintf('\n ************************************************************************ \n');
        fprintf('\nFind the reduced data in zeroed_targetMatrix and zeroed_excessVec\n');
        fprintf('\n ************************************************************************ \n');
    end
    
    fprintf('\n ');
    fprintf('\n%%%%%%%%%%%%%%%%%\n');
    fprintf('\n ');
    calib_twoSigmaALGB = standardDev'.*2;
    calib_algebraic_2Sigma = array2table(calib_twoSigmaALGB,'VariableNames',loadlist(1:dimFlag))
    
    %Should I use strtrim()  ? -AAM 042116
    %SAME START
    series_table = table([1:nseries0]','VariableNames',{'SERIES'});
    calib_algebraic_Tares = array2table(tares,'VariableNames',loadlist(1:dimFlag));
    calib_algebraic_Tares = [series_table, calib_algebraic_Tares]
    calib_algebraic_Tares_stdev = array2table(tares_STDDEV,'VariableNames',loadlist(1:dimFlag));
    calib_algebraic_Tares_stdev = [series_table, calib_algebraic_Tares_stdev]
    %SAME END
    
    coefficientsALGB = [xcalib;intercepts];
    for coeffCount=1:size(coefficientsALGB(:,1),1)
        numCombin(coeffCount) = coeffCount;
    end
    numCombinName = cellstr(num2str(numCombin'));
    numCombinName(nterms+nseries0+1) = cellstr('Intercept');
    
    %SAME START
    calib_mean_algebraic_Resids_sqrd = array2table((resSquare'./numpts0)','VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Pcnt_Capacity_Max_Mag_Load_Resids = array2table(perGoop,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Std_Dev_pcnt = array2table(stdDevPercentCapacity,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Max_Load_Resids = array2table(maxTargets,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Min_Load_Resids = array2table(minTargets,'VariableNames',loadlist(1:dimFlag))
    calib_algebraic_Ratio_Max_Mag_Load_Resid_and_Std_Dev = array2table(ratioGoop,'VariableNames',loadlist(1:dimFlag))
    
    % Prints the minmaxband
    calib_alg_per_minmaxband = array2table(theminmaxband,'VariableNames',loadlist(1:dimFlag))

    %SAME END
end

if FLAGS.excel == 1
    % Output results to an excel file
    fprintf('\n  ');
    fprintf('\nALG CALIBRATION MODEL LOAD APPROXIMATION FILE: CALIB_AOX_ALG_RESULT.csv\n'); % AJM 4_19_19
    fprintf('\n ');
    
    filename = 'CALIB_AOX_ALG_RESULT.csv';
    dlmwrite(filename,aprxIN,'precision','%.16f');
end

%SAME START
if FLAGS.res == 1
    figure('Name','Algebraic Model Calibration; Residuals of Load Versus Data Point Index','NumberTitle','off')
    plotResPages(series0, targetRes, loadCapacities, stdDevPercentCapacity, loadlist)
    %    hold off
end
%SAME END

end