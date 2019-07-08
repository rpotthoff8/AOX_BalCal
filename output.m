%Function creates all the outputs for the calibration, algebraic section
%This simplifies following the main code

function [] = output(section,FLAGS,targetRes,loadCapacities,fileName,numpts,nseries0,tares,tares_STDDEV,loadlist,series0,excessVec0,dimFlag,voltagelist,reslist,numBasis,uniqueOut)
%Split uniqueOut structure into individual variables
names = fieldnames(uniqueOut);
for i=1:length(names)
    eval([names{i} '=uniqueOut.' names{i} ';' ]);
end

% Calculates the Sum of Squares of the residual
resSquare = sum(targetRes.^2);

%STATISTIC OUTPUTS %SAME START
for k=1:length(targetRes(1,:))
    [goop(k),kstar(k)] = max(abs(targetRes(:,k)));
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
twoSigma = standardDev'.*2;

%% START PRINT OUT PERFORMANCE INFORMATION TO CSV or command window
if FLAGS.print == 1 || FLAGS.disp==1
    %Initialize cell arrays
    empty_cells=cell(1,dimFlag+1);
    Header_cells=cell(9,dimFlag+1);
    output_name=cell(1,dimFlag+1);
    load_line=[cell(1),loadlist(1:dimFlag)];
    
    %Define Header section
    Header_cells{1,1}=char(strcat(section, {' '},'Results'));
    Header_cells{2,1}=char(strcat('Performed:',{' '},datestr( datetime(now,'ConvertFrom','datenum'))));
    Header_cells{3,1}=char(strcat(strtok(section),{' '}, 'Input File:',{' '},fileName));
    if FLAGS.balOut == 1
        Header_cells{4,1}='Calibration Outliers Flagged: TRUE';
    else
        Header_cells{4,1}='Calibration Outliers Flagged: FALSE';
    end
    if FLAGS.zeroed == 1
        Header_cells{5,1}='Calibration Outliers Removed: TRUE';
    else
        Header_cells{5,1}='Calibration Outliers Removed: FALSE';
    end
    algebraic_models=[{'FULL'},{'TRUNCATED'},{'LINEAR'},{'CUSTOM'}];
    Header_cells{6,1}=char(strcat('Algebraic Model Used:',{' '},algebraic_models(FLAGS.model)));
    Header_cells{7,1}=char(strcat('Number of Datapoints:',{' '},string(numpts)));
    if FLAGS.balCal == 2
        Header_cells{8,1}='GRBF Addition Performed: TRUE';
        Header_cells{9,1}=char(strcat('Number GRBFs:',{' '},string(numBasis)));
    else
        Header_cells{8,1}='GRBF Addition Performed: FALSE';
        Header_cells{9,1}='Number GRBFs: N/A';
    end
    csv_output=[Header_cells;empty_cells];
    %Command window printing;
    if FLAGS.disp==1
        fprintf('\n ********************************************************************* \n');
        for i=1:size(Header_cells,1)
            fprintf(Header_cells{i,:})
            fprintf('\n')
        end
        fprintf('\n')
    end
    
    %Statistics output section
    
    output_name{1}='Load Residual 2*(standard deviation)';
    section_out=[load_line;cell(1),num2cell(twoSigma)];
    csv_output=[csv_output;output_name;section_out;empty_cells];
    %Command window printing;
    if FLAGS.disp==1
        fprintf(output_name{:})
        fprintf('\n')
        disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end)))
    end
    
    output_name{1}='Tares';
    section_out=[{'Series'},loadlist(1:dimFlag);num2cell([(1:nseries0)', tares])];
    csv_output=[csv_output;output_name;section_out;empty_cells];
    %Command window printing;
    if FLAGS.disp==1
        fprintf(output_name{:})
        fprintf('\n')
        disp(cell2table(section_out(2:end,1:end),'VariableNames',section_out(1,1:end)))
    end
    
    output_name{1}='Tares Standard Deviation';
    section_out=[{'Series'},loadlist(1:dimFlag);num2cell([(1:nseries0)', tares_STDDEV])];
    csv_output=[csv_output;output_name;section_out;empty_cells];
    %Command window printing;
    if FLAGS.disp==1
        fprintf(output_name{:})
        fprintf('\n')
        disp(cell2table(section_out(2:end,1:end),'VariableNames',section_out(1,1:end)))
    end
    
    output_name{1}='Mean Load Residual Squared';
    section_out=[load_line;cell(1),num2cell((resSquare'./numpts)')];
    csv_output=[csv_output;output_name;section_out;empty_cells];
    %Command window printing;
    if FLAGS.disp==1
        fprintf(output_name{:})
        fprintf('\n')
        disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end)))
    end
    
    output_name{1}='Percent Load Capacity of Maximum Residual';
    section_out=[load_line;cell(1),num2cell(perGoop)];
    csv_output=[csv_output;output_name;section_out;empty_cells];
    %Command window printing;
    if FLAGS.disp==1
        fprintf(output_name{:})
        fprintf('\n')
        disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end)))
    end
    
    output_name{1}='Percent Load Capacity of Residual Standard Deviation';
    section_out=[load_line;cell(1),num2cell(stdDevPercentCapacity)];
    csv_output=[csv_output;output_name;section_out;empty_cells];
    %Command window printing;
    if FLAGS.disp==1
        fprintf(output_name{:})
        fprintf('\n')
        disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end)))
    end
    
    output_name{1}='Maximum Load Residual';
    section_out=[load_line;cell(1),num2cell(maxTargets)];
    csv_output=[csv_output;output_name;section_out;empty_cells];
    %Command window printing;
    if FLAGS.disp==1
        fprintf(output_name{:})
        fprintf('\n')
        disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end)))
    end
    
    output_name{1}='Minimum Load Residual';
    section_out=[load_line;cell(1),num2cell(minTargets)];
    csv_output=[csv_output;output_name;section_out;empty_cells];
    %Command window printing;
    if FLAGS.disp==1
        fprintf(output_name{:})
        fprintf('\n')
        disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end)))
    end
    
    output_name{1}='Ratio (Maximum Load Residual)/(Load Residual Standard Deviation)';
    section_out=[load_line;cell(1),num2cell(ratioGoop)];
    csv_output=[csv_output;output_name;section_out;empty_cells];
    %Command window printing;
    if FLAGS.disp==1
        fprintf(output_name{:})
        fprintf('\n')
        disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end)))
    end
    
    %Print Outlier Data:
    if strcmp(section,{'Calibration Algebraic'})==1 && FLAGS.balOut==1
        output_name{1}='Outlier Information:';
        outlier_sum=cell(3,dimFlag+1);
        outlier_sum{1,1}='Standard Deviation Cutoff';
        outlier_sum{1,2}=numSTD;
        outlier_sum{2,1}='Total # Outliers';
        outlier_sum{2,2}=num_outliers;
        outlier_sum{3,1}='Total % Outliers';
        outlier_sum{3,2}=prcnt_outliers;
        csv_output=[csv_output;output_name;outlier_sum];
        %Command window printing;
        if FLAGS.disp==1
            fprintf(output_name{:})
            fprintf('\n')
            for i=1:size(outlier_sum,1)
                fprintf(1,'%s: %d',outlier_sum{i,1:2})
                fprintf('\n')
            end
            fprintf('\n')
        end
        output_name{1}='Channel Specific Outlier Summary:';
        channel_N_out=zeros(1,dimFlag);
        for i=1:dimFlag
            channel_N_out(i)=sum(colOut==i);
        end
        channel_P_out=100*channel_N_out./numpts;
        outlier_channel_sum=[{'# Outliers'},num2cell(channel_N_out);{'% Outliers'},num2cell(channel_P_out)];
        section_out=[load_line;outlier_channel_sum];
        csv_output=[csv_output;output_name;section_out;empty_cells];
        %Command window printing;
        if FLAGS.disp==1
            fprintf(output_name{:})
            fprintf('\n')
            disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end),'RowNames',section_out(2:end,1)))
        end
        
        output_name{1}='Channel Specific Outlier Indices:';
        outlier_index=cell(max(channel_N_out),dimFlag+1);
        outlier_index(:,2:end)={'-'};
        for i=1:dimFlag
            outlier_index(1:channel_N_out(i),i+1)=cellstr(num2str(rowOut(colOut==i)));
        end
        section_out=[load_line;outlier_index];
        csv_output=[csv_output;output_name;section_out;empty_cells];
        %Command window printing;
        if FLAGS.disp==1
            fprintf(output_name{:})
            fprintf('\n')
            disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end)))
        end
    end
    
    %Write Statistics to xlsx file
    if FLAGS.print==1
        warning('off', 'MATLAB:xlswrite:AddSheet'); warning('off', 'MATLAB:DELETE:FileNotFound'); warning('off',  'MATLAB:DELETE:Permission')
        filename=strcat(strtok(section),{' '},'Statistics.xlsx');
        try
            if contains(section,'Algebraic')==1
                delete(char(filename))
            end
            %Print statisitics
            sheet=contains(section,'GRBF')+1; %Set sheet to print to: 1 for algebraic, 2 for GRBF
            writetable(cell2table(csv_output),char(filename),'writevariablenames',0,'Sheet',sheet)
            %Write filename to command window
            fprintf('\n'); fprintf(char(upper(strcat(section,{' '}, 'MODEL REPORT FILE:', {' '})))); fprintf(char(filename)); fprintf(', Sheet: '); fprintf(char(num2str(sheet))); fprintf('\n')
        catch ME
            fprintf('\nUNABLE TO PRINT PERFORMANCE PARAMETER XLSX FILE. ');
            if (strcmp(ME.identifier,'MATLAB:table:write:FileOpenInAnotherProcess'))
                fprintf('ENSURE "'); fprintf(char(filename));fprintf('" IS NOT OPEN AND TRY AGAIN')
            end
            fprintf('\n')
        end
        warning('on',  'MATLAB:DELETE:Permission'); warning('on', 'MATLAB:xlswrite:AddSheet'); warning('on', 'MATLAB:DELETE:FileNotFound')
    end
    
end

%% Residual vs datapoint plot
if FLAGS.res == 1
    figure('Name',char(strcat(section,{' '},'Model; Residuals of Load Versus Data Point Index')),'NumberTitle','off','WindowState','maximized')
    plotResPages(series0, targetRes, loadCapacities, stdDevPercentCapacity, loadlist)
end

%% OUTPUT HISTOGRAM PLOTS
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

%% Prints residual vs. input and calculates correlations
if FLAGS.rescorr == 1
    figure('Name',char(strcat(section,{' '},'Residual Correlation Plot')),'NumberTitle','off','WindowState','maximized');
    correlationPlot(excessVec0, targetRes, voltagelist, reslist);
end

%% Algebraic Validation Specific Outputs
if strcmp(section,{'Validation Algebraic'})==1
    if FLAGS.excel == 1
        filename = 'VALID_AOX_GLOBAL_ALG_RESULT.csv';
        input=aprxINminGZvalid;
        precision='%.16f';
        description='VALIDATION ALGEBRAIC MODEL GLOBAL LOAD APPROXIMATION';
        print_dlmwrite(filename,input,precision,description);
    end
end

%% Algebraic Calibration Specific Outputs
if strcmp(section,{'Calibration Algebraic'})==1
    
    %Prints coefficients to csv file
    if FLAGS.excel == 1
        filename = 'APPROX_AOX_COEFF_MATRIX.csv';
        input=[coeff;zeros(1,dimFlag)];
        precision='%.16f';
        description='CALIBRATION ALGEBRAIC MODEL COEFFICIENT MATRIX';
        print_dlmwrite(filename,input,precision,description);
    end
    
    %%% ANOVA Stats AJM 6_12_19
    if FLAGS.model ~= 4 && FLAGS.anova==1
        
        totalnum = nterms+nseries0;
        totalnumcoeffs = [1:totalnum];
        totalnumcoeffs2 = [2:totalnum+1];
        dsof = numpts-nterms-1;
        
        loadstatlist = {'Load', 'Sum_Sqrs', 'PRESS_Stat', 'DOF', 'Mean_Sqrs', 'F_Value', 'P_Value', 'R_sq', 'Adj_R_sq', 'PRESS_R_sq'};
        regresslist = {'Term', 'Coeff_Value', 'CI_95cnt', 'T_Stat', 'P_Value', 'VIF_A', 'Signif'};
        
        for k=1:dimFlag
            RECOMM_ALG_EQN(:,k) = [1.0*ANOVA(k).sig([1:nterms])];
            manoa2(k,:) = [loadlist(k), tR2(1,k), ANOVA(k).PRESS, dsof, gee(1,k), ANOVA(k).F, ANOVA(k).p_F, ANOVA(k).R_sq, ANOVA(k).R_sq_adj, ANOVA(k).R_sq_p];
            ANOVA01(:,:) = [totalnumcoeffs; ANOVA(k).beta'; ANOVA(k).beta_CI'; ANOVA(k).T'; ANOVA(k).p_T'; ANOVA(k).VIF'; 1.0*ANOVA(k).sig']';
            ANOVA1_2(:,:) = [ANOVA01([1:nterms],:)];
            STAT_LOAD = array2table(manoa2(k,:),'VariableNames',loadstatlist(1:10));
            REGRESS_COEFFS = array2table(ANOVA1_2(:,:),'VariableNames',regresslist(1:7));
        end
        
        warning('off', 'MATLAB:xlswrite:AddSheet')
        filename = 'DIRECT_ANOVA_STATS.xlsx';
        try
            for k=1:dimFlag
                writetable(STAT_LOAD,filename,'Sheet',k,'Range','A1');
                writetable(REGRESS_COEFFS,filename,'Sheet',k,'Range','A4');
            end
            fprintf('\nDIRECT METHOD ANOVA STATISTICS FILE: '); fprintf(filename); fprintf('\n ');
        catch ME
            fprintf('\nUNABLE TO PRINT DIRECT METHOD ANOVA STATISTICS FILE. ');
            if (strcmp(ME.identifier,'MATLAB:table:write:FileOpenInAnotherProcess'))
                fprintf('ENSURE "'); fprintf(char(filename)); fprintf('" IS NOT OPEN AND TRY AGAIN')
            end
            fprintf('\n')
        end
        warning('on', 'MATLAB:xlswrite:AddSheet')
        
        %Output recommended custom equation
        if FLAGS.Rec_Model==1
            filename = 'DIRECT_RECOMM_CustomEquationMatrix.csv';
            input=RECOMM_ALG_EQN;
            precision='%.8f';
            description='DIRECT METHOD ANOVA RECOMMENDED CUSTOM EQUATION MATRIX';
            print_dlmwrite(filename,input,precision,description);
        end
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
        end
        
        if FLAGS.BALFIT_ANOVA==1
            warning('off', 'MATLAB:xlswrite:AddSheet')
            filename = 'BALFIT_ANOVA_STATS.xlsx';
            try
                for k=1:dimFlag
                    writetable(BALFIT_STAT_VOLTAGE_1,filename,'Sheet',k,'Range','A1');
                    writetable(BALFIT_REGRESS_COEFFS_1,filename,'Sheet',k,'Range','A4');
                end
                fprintf('\nBALFIT ANOVA STATISTICS FILE: '); fprintf(filename); fprintf('\n');
                %filename = 'BALFIT_RECOMM_CustomEquationMatrixTemplate.csv';
                %dlmwrite(filename,BALFIT_RECOMM_ALG_EQN,'precision','%.8f');
            catch ME
                fprintf('\nUNABLE TO PRINT BALFIT ANOVA STATISTICS FILE. ');
                if (strcmp(ME.identifier,'MATLAB:table:write:FileOpenInAnotherProcess'))
                    fprintf('ENSURE "'); fprintf(char(filename)); fprintf('" IS NOT OPEN AND TRY AGAIN');
                end
                fprintf('\n')
            end
            warning('on', 'MATLAB:xlswrite:AddSheet')
        end
        
    end
    
    if FLAGS.anova==1
        filename = 'Testing_Sig.csv';
        input=balfitANOVA(1).sig;
        precision='%.16f';
        description='BALFIT TESTING SIG';
        print_dlmwrite(filename,input,precision,description);
    end
    
    if FLAGS.BALFIT_Matrix==1
        filename = 'BALFIT_DATA_REDUCTION_MATRIX_IN_AMES_FORMAT.csv';
        input=balfit_regress_matrix;
        precision='%.16f';
        description='BALFIT DATA REDUCTION MATRIX IN AMES FORMAT';
        print_dlmwrite(filename,input,precision,description);
        %%% Balfit Stats and Matrix AJM 5_31_19
    end
    
    
    if FLAGS.excel == 1
        %Output calibration load approximation
        filename = 'CALIB_AOX_ALG_RESULT.csv';
        input=aprxIN;
        precision='%.16f';
        description='CALIBRATION ALGEBRAIC MODEL LOAD APPROXIMATION';
        print_dlmwrite(filename,input,precision,description);
    end
end

%% GRBF Calibration Specific Outputs
if strcmp(section,{'Calibration GRBF'})==1
    if FLAGS.excel == 1
        %Output calibration load approximation
        filename = 'CALIB_AOX_GRBF_RESULT.csv';
        input=aprxINminGZ2;
        precision='%.16f';
        description='CALIBRATION ALGEBRAIC+GRBF MODEL LOAD APPROXIMATION';
        print_dlmwrite(filename,input,precision,description);
        
        %Output GRBF Widths
        filename = 'APPROX_AOX_GRBF_ws.csv';
        input=wHist;
        precision='%.16f';
        description='CALIBRATION GRBF WIDTHS';
        print_dlmwrite(filename,input,precision,description);
        
        %Output GRBF coefficients
        filename = 'APPROX_AOX_GRBF_coeffs.csv';
        input=cHist;
        precision='%.16f';
        description='CALIBRATION GRBF COEFFICIENTS';
        print_dlmwrite(filename,input,precision,description);
        
        %Output GRBF centers
        filename = 'APPROX_AOX_GRBF_Centers.csv';
        input=centerIndexHist;
        precision='%.16f';
        description='CALIBRATION GRBF CENTER INDICES';
        print_dlmwrite(filename,input,precision,description);
    end
end
%% GRBF Validation Specific Outputs
if strcmp(section,{'Validation GRBF'})==1
    if FLAGS.excel == 1
        %Output validation load approximation
        filename = 'VALID_AOX_GLOBAL_GRBF_RESULT.csv';
        input=aprxINminGZ2valid;
        precision='%.16f';
        description='VALIDATION ALGEBRAIC+GRBF MODEL LOAD APPROXIMATION';
        print_dlmwrite(filename,input,precision,description);
    end
end

end