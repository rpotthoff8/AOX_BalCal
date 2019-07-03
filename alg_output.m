%Function creates all the outputs for the calibration, algebraic section
%This simplifies following the main code

function [] = alg_output(section,FLAGS,targetRes,loadCapacities,fileName,numpts0,nseries0,tares,tares_STDDEV,loadlist,aprxIN,series0,excessVec0,dimFlag,coeff,voltagelist,reslist,nterms,ANOVA,balfitcomIN,balfitxcalib,balfittargetMatrix,balfitANOVA)

% Calculates the Sum of Squares of the residual
resSquare = sum(targetRes.^2);

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
%    theminmaxband = 100*(abs(maxTargets + minTargets)./loadCapacities); %QUESTION: JRP 2 July 19
theminmaxband = 100*(abs(maxTargets - minTargets)./loadCapacities);
%SAME END

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
%    theminmaxband = 100*(abs(maxTargets + minTargets)./loadCapacities); %QUESTION: JRP 2 July 19
theminmaxband = 100*(abs(maxTargets - minTargets)./loadCapacities);
%SAME END

%START PRINT OUT PERFORMANCE INFORMATION TO CSV
if FLAGS.print == 1
    empty_cells=cell(1,dimFlag+1);
    Header_cells=cell(8,dimFlag+1);
    Header_cells{1,1}=strcat(section, {' '},'Results');
    Header_cells{2,1}=strcat('Performed:',{' '},datestr( datetime(now,'ConvertFrom','datenum')));
    Header_cells{3,1}=strcat('Calibration Input File:',{' '},fileName);
    if FLAGS.balCal == 2
        Header_cells{4,1}='GRBF Addition Performed: TRUE';
    else
        Header_cells{4,1}='GRBF Addition Performed: FALSE';
    end
    if FLAGS.balOut == 1
        Header_cells{5,1}='Outliers Flagged: TRUE';
    else
        Header_cells{5,1}='Outliers Flagged: FALSE';
    end
    if FLAGS.zeroed == 1
        Header_cells{6,1}='Outliers Removed: TRUE';
    else
        Header_cells{6,1}='Outliers Removed: FALSE';
    end
    algebraic_models=[{'Full'},{'Truncated'},{'Linear'},{'Custom'}];
    Header_cells{7,1}=strcat('Algebraic Model Used:',{' '},algebraic_models(FLAGS.model));
    Header_cells{8,1}=strcat('Number of Datapoints:',{' '},string(numpts0));
    
    calib_output=[Header_cells;empty_cells];
    
    calib_twoSigmaALGB = standardDev'.*2;
    %     calib_algebraic_2Sigma = array2table(calib_twoSigmaALGB,'VariableNames',loadlist(1:dimFlag))
    output_name=cell(1,dimFlag+1);
    load_line=[cell(1),loadlist(1:dimFlag)];
    output_name{1}='Load Residual 2*(standard deviation)';
    calib_output=[calib_output;output_name;load_line;cell(1),num2cell(calib_twoSigmaALGB);empty_cells];
    
    %Should I use strtrim()  ? -AAM 042116
    %SAME START
    %     series_table = table([1:nseries0]','VariableNames',{'SERIES'});
    %     calib_algebraic_Tares = array2table(tares,'VariableNames',loadlist(1:dimFlag));
    %     calib_algebraic_Tares = [series_table, calib_algebraic_Tares]
    output_name{1}='Tares';
    calib_output=[calib_output;output_name;{'Series'},loadlist(1:dimFlag);num2cell([(1:nseries0)', tares]);empty_cells];
    
    %     calib_algebraic_Tares_stdev = array2table(tares_STDDEV,'VariableNames',loadlist(1:dimFlag));
    output_name{1}='Tares Standard Deviation';
    calib_output=[calib_output;output_name;{'Series'},loadlist(1:dimFlag);num2cell([[1:nseries0]', tares_STDDEV]);empty_cells];
    %     calib_algebraic_Tares_stdev = [series_table, calib_algebraic_Tares_stdev]
    %SAME END
    
    
    %SAME START
    %    calib_mean_algebraic_Resids_sqrd = array2table((resSquare'./numpts0)','VariableNames',loadlist(1:dimFlag))
    output_name{1}='Mean Load Residual Squared';
    calib_output=[calib_output;output_name;load_line;cell(1),num2cell((resSquare'./numpts0)');empty_cells];
    
    %    calib_algebraic_Pcnt_Capacity_Max_Mag_Load_Resids = array2table(perGoop,'VariableNames',loadlist(1:dimFlag))
    output_name{1}='Percent Load Capacity of Maximum Residual';
    calib_output=[calib_output;output_name;load_line;cell(1),num2cell(perGoop);empty_cells];
    
    %    calib_algebraic_Std_Dev_pcnt = array2table(stdDevPercentCapacity,'VariableNames',loadlist(1:dimFlag))
    output_name{1}='Percent Load Capacity of Residual Standard Deviation';
    calib_output=[calib_output;output_name;load_line;cell(1),num2cell(stdDevPercentCapacity);empty_cells];
    
    %    calib_algebraic_Max_Load_Resids = array2table(maxTargets,'VariableNames',loadlist(1:dimFlag))
    output_name{1}='Maximum Load Residual';
    calib_output=[calib_output;output_name;load_line;cell(1),num2cell(maxTargets);empty_cells];
    
    %    calib_algebraic_Min_Load_Resids = array2table(minTargets,'VariableNames',loadlist(1:dimFlag))
    output_name{1}='Minimum Load Residual';
    calib_output=[calib_output;output_name;load_line;cell(1),num2cell(minTargets);empty_cells];
    
    %    calib_algebraic_Ratio_Max_Mag_Load_Resid_and_Std_Dev = array2table(ratioGoop,'VariableNames',loadlist(1:dimFlag))
    output_name{1}='Ratio (Maximum Load Residual)/(Load Residual Standard Deviation)';
    calib_output=[calib_output;output_name;load_line;cell(1),num2cell(ratioGoop);empty_cells];
    
    % Prints the minmaxband
    %    calib_alg_per_minmaxband = array2table(theminmaxband,'VariableNames',loadlist(1:dimFlag))
    output_name{1}='Percent Load Capacity of Band Between Min and Max Residual';
    calib_output=[calib_output;output_name;load_line;cell(1),num2cell(theminmaxband);empty_cells];
    
    output_file=strcat(section,{' '},'Results.csv');
    writetable(cell2table(calib_output),char(output_file),'writevariablenames',0)
    %SAME END
end

%SAME START
if FLAGS.res == 1
    figure('Name',char(strcat(section,{' '},'Model; Residuals of Load Versus Data Point Index')),'NumberTitle','off','WindowState','maximized')
    plotResPages(series0, targetRes, loadCapacities, stdDevPercentCapacity, loadlist)
    %    hold off
end
%SAME END

%SAME START
%OUTPUT HISTOGRAM PLOTS
if FLAGS.hist == 1
    figure('Name',char(section),'NumberTitle','off','WindowState','maximized')
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

% Prints residual vs. input and calculates correlations
if FLAGS.rescorr == 1
    figure('Name',char(strcat(section,{' '},'Residual Correlation Plot')),'NumberTitle','off','WindowState','maximized');
    correlationPlot(excessVec0, targetRes, voltagelist, reslist);
end


if strcmp(section,{'Validation Algebraic'})==1
    if FLAGS.excel == 1
        %%%%
        fprintf('\nALG VALIDATION MODEL GLOBAL LOAD APPROXIMATION: VALID_AOX_GLOBAL_ALG_RESULT in Workspace\n');
        fprintf('\n ');
        
        filename = 'VALID_AOX_GLOBAL_ALG_RESULT.csv';
        %        csvwrite(filename,aprxINminGZvalid)
        dlmwrite(filename,aprxIN,'precision','%.16f');
    end
    
elseif strcmp(section,{'Calibration Algebraic'})==1
    
    %Prints coefficients to csv file
    if FLAGS.excel == 1 %ADD HEADER?
        filename = 'APPROX_AOX_COEFF_MATRIX.csv';
        dlmwrite(filename, [coeff;zeros(1,dimFlag)] ,'precision','%.16f'); % AJM 4_18_19
    end
    
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
    
    if FLAGS.excel == 1
        % Output results to an excel file
        fprintf('\n  ');
        fprintf('\nALG CALIBRATION MODEL LOAD APPROXIMATION FILE: CALIB_AOX_ALG_RESULT.csv\n'); % AJM 4_19_19
        fprintf('\n ');
        
        filename = 'CALIB_AOX_ALG_RESULT.csv';
        dlmwrite(filename,aprxIN,'precision','%.16f');
    end
    
end

end
