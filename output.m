function [] = output(section,FLAGS,targetRes,fileName,numpts,nseries0,loadlist,series,excessVec0,voltdimFlag,loaddimFlag,voltagelist,reslist,numBasis,pointID,series2,output_location,REPORT_NO,algebraic_model,uniqueOut)
%Function creates all the outputs for the calibration and validation ALG
%and GRBF sections.  First common section runs for all sections.  Following
%sections run for specific section outputs
%This simplifies following the main code

%INPUTS:
%  section = current section of code. Expected values: 'Calibration Algebraic', 'Calibration GRBF', 'Validation Algebraic', 'Validation GRBF'
%  FLAGS = Structure containing flags from user inputs
%  targetRes = Matrix of residuals
%  loadCapacities = Vector of load capacities in each channel
%  fileName = Data Input filename
%  numpts = Number of datapoints (observations)
%  nseries0 = Number of series
%  tares = Matrix of calculated tares for each series
%  tares_STDDEV = Matrix of standard deviation of residuals used to calculate tares
%  loadlist = Labels for load variables
%  series = Series1 labels
%  excessVec0 = Matrix of measured voltages
%  voltdimFlag = Dimension of voltage data (# channels)
%  loaddimFlag = Dimension of load data (# channels)
%  voltagelist = Labels for voltage variables
%  reslist = Labels for residuals
%  numBasis = Number of RBFs to include in model
%  pointID = Point ID for each observation, from data input file
%  series2 = Series2 labels
%  output_location = Path for output files
%  REPORT_NO = Report number for run.
%  algebraic_model = Type of algebraic model used
%  uniqueOut = Structure containing unique variables to be output based on specific code section

%OUTPUTS:
%  []

%Split uniqueOut structure into individual variables
names = fieldnames(uniqueOut);
for i=1:length(names)
    eval([names{i} '=uniqueOut.' names{i} ';' ]);
end

% Calculates the Sum of Squares of the residual
resSquare = sum(targetRes.^2);

%STATISTIC OUTPUTS
for k=1:length(targetRes(1,:))
    [goop(k),kstar(k)] = max(abs(targetRes(:,k)));
    maxTargets(k) = max(targetRes(:,k));
    minTargets(k) = min(targetRes(:,k));
    tR2(k) = targetRes(:,k)'*targetRes(:,k);     % AJM 6_12_19
end
davariance = var(targetRes);
gee = mean(targetRes);
standardDev10 = std(targetRes);
standardDev = standardDev10';
ratioGoop = goop./standardDev';
ratioGoop(isnan(ratioGoop)) = realmin;
twoSigma = standardDev'.*2;
if FLAGS.mode==1
    perGoop = 100*(goop./loadCapacities);
    stdDevPercentCapacity = 100*(standardDev'./loadCapacities);
end

%% START PRINT OUT PERFORMANCE INFORMATION TO CSV or command window
if FLAGS.print == 1 || FLAGS.disp==1
    %Initialize cell arrays
    empty_cells=cell(1,loaddimFlag+1);
    Header_cells=cell(10,loaddimFlag+1);
    output_name=cell(1,loaddimFlag+1);
    load_line=[cell(1),loadlist(1:loaddimFlag)];
    
    %Define Header section
    Header_cells{1,1}=char(strcat(section, {' '},'Results'));
    if FLAGS.mode==1
        Header_cells{2,1}='Software Mode: FORCE BALANCE CALIBRATION';
    elseif FLAGS.mode==2
        Header_cells{2,1}='Software Mode: GENERAL FUNCTION APPROXIMATION';
    end
    Header_cells{3,1}=char(strcat('REPORT NO:',{' '},REPORT_NO));
    Header_cells{4,1}=char(strcat(strtok(section),{' '}, 'Input File:',{' '},fileName));
    if FLAGS.balOut == 1
        Header_cells{5,1}='Calibration ALG Outliers Flagged: TRUE';
    else
        Header_cells{5,1}='Calibration ALG Outliers Flagged: FALSE';
    end
    if FLAGS.zeroed == 1
        Header_cells{6,1}='Calibration ALG Outliers Removed: TRUE';
    else
        Header_cells{6,1}='Calibration ALG Outliers Removed: FALSE';
    end
    Header_cells{7,1}=char(strcat('Algebraic Model Used:',{' '},algebraic_model));
    Header_cells{8,1}=char(strcat('Number of Datapoints:',{' '},string(numpts)));
    if FLAGS.balCal == 2
        Header_cells{9,1}='GRBF Addition Performed: TRUE';
        Header_cells{10,1}=char(strcat('Number GRBFs:',{' '},string(numBasis)));
    else
        Header_cells{9,1}='GRBF Addition Performed: FALSE';
        Header_cells{10,1}='Number GRBFs: N/A';
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
      
    %Call outputs 'Load' for balance calibration, otherwise only 'Load'
    if FLAGS.mode==1
        outLabel='Load';
    else
        outLabel='Output';
    end
    
    if FLAGS.mode==1
        %Statistics output section
        output_name{1}='Percent Load Capacity of Residual Standard Deviation';
        section_out=[load_line;cell(1),num2cell(stdDevPercentCapacity)];
        csv_output=[csv_output;output_name;section_out;empty_cells];
        %Command window printing;
        if FLAGS.disp==1
            fprintf(output_name{:})
            fprintf('\n')
            disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end)))
        end
        
        output_name{1}='Tares';
        section_out=[{'Series'},loadlist(1:loaddimFlag);num2cell([(1:nseries0)', tares])];
        csv_output=[csv_output;output_name;section_out;empty_cells];
        %Command window printing;
        if FLAGS.disp==1
            fprintf(output_name{:})
            fprintf('\n')
            disp(cell2table(section_out(2:end,1:end),'VariableNames',section_out(1,1:end)))
        end
        
        output_name{1}='Tares Standard Deviation';
        section_out=[{'Series'},loadlist(1:loaddimFlag);num2cell([(1:nseries0)', tares_STDDEV])];
        csv_output=[csv_output;output_name;section_out;empty_cells];
        %Command window printing;
        if FLAGS.disp==1
            fprintf(output_name{:})
            fprintf('\n')
            disp(cell2table(section_out(2:end,1:end),'VariableNames',section_out(1,1:end)))
        end
    end
    
    output_name{1}=['Mean ' outLabel,' Residual Squared'];
    section_out=[load_line;cell(1),num2cell((resSquare'./numpts)')];
    csv_output=[csv_output;output_name;section_out;empty_cells];
    %Command window printing;
    if FLAGS.disp==1
        fprintf(output_name{:})
        fprintf('\n')
        disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end)))
    end
    
    if FLAGS.mode==1
        output_name{1}='Percent Load Capacity of Maximum Residual';
        section_out=[load_line;cell(1),num2cell(perGoop)];
        csv_output=[csv_output;output_name;section_out;empty_cells];
        %Command window printing;
        if FLAGS.disp==1
            fprintf(output_name{:})
            fprintf('\n')
            disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end)))
        end
    end
    
    output_name{1}=['Maximum ',outLabel,' Residual'];
    section_out=[load_line;cell(1),num2cell(maxTargets)];
    csv_output=[csv_output;output_name;section_out;empty_cells];
    %Command window printing;
    if FLAGS.disp==1
        fprintf(output_name{:})
        fprintf('\n')
        disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end)))
    end
    
    output_name{1}=['Minimum ',outLabel,' Residual'];
    section_out=[load_line;cell(1),num2cell(minTargets)];
    csv_output=[csv_output;output_name;section_out;empty_cells];
    %Command window printing;
    if FLAGS.disp==1
        fprintf(output_name{:})
        fprintf('\n')
        disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end)))
    end
    
    output_name{1}=[outLabel,' Residual 2*(standard deviation)'];
    section_out=[load_line;cell(1),num2cell(twoSigma)];
    csv_output=[csv_output;output_name;section_out;empty_cells];
    %Command window printing;
    if FLAGS.disp==1
        fprintf(output_name{:})
        fprintf('\n')
        disp(cell2table(section_out(2:end,2:end),'VariableNames',section_out(1,2:end)))
    end
    
    output_name{1}=['Ratio (Maximum ',outLabel,' Residual)/(',outLabel,' Residual Standard Deviation)'];
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
        outlier_sum=cell(3,loaddimFlag+1);
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
        channel_N_out=zeros(1,loaddimFlag);
        for i=1:loaddimFlag
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
        outlier_index=cell(max(channel_N_out),loaddimFlag+1);
        outlier_index(:,2:end)={'-'};
        for i=1:loaddimFlag
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
        if contains(section,'Calibration')
            filename='CALIB Report.xlsx';
        elseif contains(section,'Validation')
            filename='VALID Report.xlsx';
        end
        fullpath=fullfile(output_location,filename);
        try
            if contains(section,'Algebraic')==1
                delete(char(fullpath))
            end
            %Print statisitics
            sheet=contains(section,'GRBF')+1; %Set sheet to print to: 1 for algebraic, 2 for GRBF
            writetable(cell2table(csv_output),fullpath,'writevariablenames',0,'Sheet',sheet,'UseExcel', false)
            %Write filename to command window
            fprintf('\n'); fprintf(char(upper(strcat(section,{' '}, 'MODEL REPORT FILE:', {' '})))); fprintf(char(filename)); fprintf(', Sheet: '); fprintf(char(num2str(sheet))); fprintf('\n')
        catch ME
            fprintf('\nUNABLE TO PRINT PERFORMANCE PARAMETER XLSX FILE. ');
            if (strcmp(ME.identifier,'MATLAB:table:write:FileOpenInAnotherProcess'))
                fprintf('ENSURE "'); fprintf(char(filename));fprintf('" IS NOT OPEN AND TRY AGAIN')
            end
            fprintf('\n')
        end
        
        try %Rename excel sheets and delete extra sheets, only possible on PC
            [~,sheets]=xlsfinfo(fullpath);
            s = what;
            e = actxserver('Excel.Application'); % # open Activex server
            e.DisplayAlerts = false;
            e.Visible=false;
            ewb = e.Workbooks.Open(char(fullpath)); % # open file (enter full path!)
            
            if contains(section,'Algebraic')==1
                ewb.Worksheets.Item(sheet).Name = 'Algebraic Results'; % # rename 1st sheet
                %cycle through, deleting all sheets other than the 1st (Algebraic) sheet
                for i=sheet+1:max(size(sheets))
                    ewb.Sheets.Item(sheet+1).Delete;  %Delete 2nd sheet
                end
            elseif contains(section,'GRBF')==1
                ewb.Worksheets.Item(sheet).Name = 'GRBF Results'; % # rename 2nd sheet
            end
            ewb.Save % # save to the same file
            ewb.Close
            e.Quit
            delete(e);
        end
        warning('on',  'MATLAB:DELETE:Permission'); warning('on', 'MATLAB:xlswrite:AddSheet'); warning('on', 'MATLAB:DELETE:FileNotFound')
    end
    
end

%% Residual vs datapoint plot
if FLAGS.res == 1
    if FLAGS.mode==1
         figure('Name',char(strcat(section,{' '},'Model; Residuals of Load Versus Data Point Index')),'NumberTitle','off','WindowState','maximized')
        plotResPages(series, targetRes, loadlist, stdDevPercentCapacity, loadCapacities)
    else
        figure('Name',char(strcat(section,{' '},'Model; Residuals of Output Versus Data Point Index')),'NumberTitle','off','WindowState','maximized')
        plotResPages(series, targetRes, loadlist, standardDev)
    end
end

%% OUTPUT HISTOGRAM PLOTS
if FLAGS.hist == 1
    figure('Name',char(section),'NumberTitle','off','WindowState','maximized')
    for k0=1:length(targetRes(1,:))
        subplot(2,ceil(loaddimFlag/2),k0)
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


%% Algebraic Calibration Specific Outputs
if strcmp(section,{'Calibration Algebraic'})==1
    Term_Names=customMatrix_labels(loadlist,voltagelist,voltdimFlag,loaddimFlag,FLAGS.model,'voltages'); %Get label names for custom equation matrix
    %Prints coefficients to csv file
    if FLAGS.excel == 1
        filename = 'AOX_ALG_MODEL_COEFFICIENT_MATRIX.csv';
        description='CALIBRATION ALGEBRAIC MODEL COEFFICIENT MATRIX';
        print_coeff(filename,coeff,description,Term_Names,loadlist,output_location)
    end
    
    %%% ANOVA Stats AJM 6_12_19
    if FLAGS.anova==1
        
        RECOMM_ALG_EQN=anova_output(ANOVA,nterms,numpts,loaddimFlag,'CALIB ALG',output_location,Term_Names,loadlist,tR2,gee);
        
        %Output recommended custom equation
        if FLAGS.Rec_Model==1
            filename = 'DIRECT_RECOMM_EquationMatrix.csv';
            fullpath=fullfile(output_location,filename);
            [leftColumn,topRow]=customMatrix_labels(loadlist,voltagelist,voltdimFlag,loaddimFlag,FLAGS.model,'voltages'); %Get label names for custom equation matrix
            recTable=array2table(RECOMM_ALG_EQN,'VariableNames',topRow,'RowNames',leftColumn(:));
            description='DIRECT METHOD ANOVA RECOMMENDED EQUATION MATRIX';
            try
                writetable(recTable,fullpath,'WriteRowNames',true);
                fprintf('\n'); fprintf(description); fprintf(' FILE: '); fprintf(filename); fprintf('\n');
            catch ME
                fprintf('\nUNABLE TO PRINT '); fprintf('%s %s', upper(description),'FILE. ');
                if (strcmp(ME.identifier,'MATLAB:table:write:FileOpenInAnotherProcess'))
                    fprintf('ENSURE "'); fprintf(char(filename)); fprintf('" IS NOT OPEN AND TRY AGAIN')
                end
                fprintf('\n')
            end
        end
              
    end
    %%% ANOVA Stats AJM 6_8_19
    
    if FLAGS.calc_balfit==1
        %%% Balfit Stats and Regression Coeff Matrix AJM 5_31_19
        totalnum = size(balfitxcalib,1);
        totalnumcoeffs = [1:totalnum];
        totalnumcoeffs2 = [2:totalnum+1];
        dsof = numpts-nterms-1;
        
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
        balfitregresslist = {'Term_Index','Term_Name', 'Coeff_Value', 'CI_95cnt', 'T_Stat', 'P_Value', 'VIF_A', 'Signif'};
        
        %balfitinterceptlist = ['Intercept', '0', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A'];
        balfitinterceptlist = [1, 0, 0, 0, 0, 0, 0];
        
        BALFIT_STAT_VOLTAGE_1=cell(loaddimFlag,1);
        BALFIT_REGRESS_COEFFS_1=cell(loaddimFlag,1);
        Term_Names=customMatrix_labels(loadlist,voltagelist,voltdimFlag,loaddimFlag,FLAGS.model,'loads'); %Get label names for custom equation matrix
        if FLAGS.anova==1
            for k=1:loaddimFlag
                BALFIT_RECOMM_ALG_EQN(:,k) = 1.0*balfitANOVA(k).sig;
                balfitANOVA01(:,:) = [totalnumcoeffs2; balfitANOVA(k).beta'; ANOVA(k).beta_CI'; balfitANOVA(k).T'; balfitANOVA(k).p_T'; balfitANOVA(k).VIF'; 1.0*balfitANOVA(k).sig']';
                balfitANOVA_intercept1(1,:) = balfitinterceptlist(1,:);
                balfitANOVA1([1:nterms+1],:) = num2cell([balfitANOVA_intercept1(1,:); balfitANOVA01([1:nterms],:)]);
                toplayer2(k,:) = [voltagelist(k), balfittR2(1,k), balfitANOVA(k).PRESS, dsof, balfitgee(1,k), balfitANOVA(k).F, balfitANOVA(k).p_F, balfitANOVA(k).R_sq, balfitANOVA(k).R_sq_adj, balfitANOVA(k).R_sq_p];
                BALFIT_STAT_VOLTAGE_1{k} = array2table(toplayer2(k,:),'VariableNames',voltagestatlist(1:10));
                BALFIT_REGRESS_COEFFS_1{k} = cell2table([balfitANOVA1(1:nterms+1,1),[{'INTERCEPT'};Term_Names],balfitANOVA1(1:nterms+1,2:end)],'VariableNames',balfitregresslist);
            end
            
            if FLAGS.BALFIT_ANOVA==1
                warning('off', 'MATLAB:xlswrite:AddSheet'); warning('off', 'MATLAB:DELETE:FileNotFound'); warning('off',  'MATLAB:DELETE:Permission')
                filename = 'BALFIT_ANOVA_STATS.xlsx';
                fullpath=fullfile(output_location,filename);
                try
                    delete(char(fullpath))
                    for k=1:loaddimFlag
                        writetable(BALFIT_STAT_VOLTAGE_1{k},fullpath,'Sheet',k,'Range','A1');
                        writetable(BALFIT_REGRESS_COEFFS_1{k},fullpath,'Sheet',k,'Range','A4');
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
                warning('on',  'MATLAB:DELETE:Permission'); warning('on', 'MATLAB:xlswrite:AddSheet'); warning('on', 'MATLAB:DELETE:FileNotFound')
            end
            
        end
        
        if FLAGS.BALFIT_Matrix==1
            [leftColumn_coeff,voltRow]=customMatrix_labels(loadlist,voltagelist,voltdimFlag,loaddimFlag,FLAGS.model,'loads'); %Get label names for custom equation matrix
            Header_cells=cell(18,loaddimFlag);
            dash_row=cell(1,loaddimFlag);
            dash_row(:,:)={'- - - - - - - -'};
            dash_row{1,1}=repmat('- ',1,round(15*loaddimFlag/2));
            Header_cells(1,1)={'DATA_REDUCTION_MATRIX_IN_AMES_FORMAT'};
            Header_cells(2,1)={balance_type};
            Header_cells(3,1)={description};
            Header_cells(4,1)={'AOX_BalCal'};
            Header_cells(5,1)={REPORT_NO};
            Header_cells(6,1)={'Primary Load Iteration Method'};
            Header_cells(7,:)=loadlist(1:loaddimFlag);
            Header_cells(8,:)=loadunits(1:loaddimFlag);
            Header_cells(9,:)=num2cell(min(targetMatrix0));
            Header_cells(10,:)=num2cell(max(targetMatrix0));
            Header_cells(11,:)=num2cell(loadCapacities);
            Header_cells(12,:)=voltRow;
            Header_cells(13,:)=voltunits(1:loaddimFlag)';
            Header_cells(14,:)=num2cell(min(excessVec0));
            Header_cells(15,:)=num2cell(max(excessVec0));
            Header_cells(16,:)={'DEFINED'};
            Header_cells(17:18,:)=num2cell(balfit_regress_matrix(1:2,:));
            leftColumn_head=[{'FILE_TYPE'};{'BALANCE_NAME'};{'DESCRIPTION'};{'PREPARED_BY'};{'REPORT_NO'};{'ITERATION_METHOD'};{'LOAD_NAME'};{'LOAD_UNIT'};{'LOAD_MINIMUM'};{'LOAD_MAXIMUM'};{'LOAD_CAPACITY'};{'GAGE_OUT_NAME'};{'GAGE_OUT_UNIT'};{'GAGE_OUT_MINIMUM'};{'GAGE_OUT_MAXIMUM'};{'GAGE_SENSITIVITY'};{'NATURAL_ZERO'};{'INTERCEPT'}];
            leftColumn=[leftColumn_head;{'D0[TRANSPONSE (C1INV)]'};voltRow;{'D1[MATRIX IS NOT USED]'};leftColumn_coeff(1:loaddimFlag);'D2[TRANSPOSE(C1INVC2)]';leftColumn_coeff(loaddimFlag+1:end)];
            D0=balfit_regress_matrix(3:2+loaddimFlag,:);
            D1=balfit_regress_matrix(3+loaddimFlag:2+2*loaddimFlag,:);
            D2=balfit_regress_matrix(3+2*loaddimFlag:end,:);
            
            content=[Header_cells;dash_row;num2cell(D0);dash_row;num2cell(D1);dash_row;num2cell(D2)];
            balfit_matrix=[leftColumn,content];
            
            % Text file to output data
            filename = 'BALFIT_DATA_REDUCTION_MATRIX_IN_AMES_FORMAT.txt';
            fullpath=fullfile(output_location,filename);
            description='BALFIT DATA REDUCTION MATRIX IN AMES FORMAT';
            
            try
                % Open file for writing
                fid = fopen(fullpath, 'w');
                header_lines=[1:6];
                text_lines=[7,8,12,13,16];
                dash_lines=[19,19+loaddimFlag+1,19+2*loaddimFlag+2];
                
                for i=1:size(content,1) %Printing to text file
                    fprintf(fid,'%23s',leftColumn{i});
                    if ismember(i,header_lines) %Printing header cells
                        fprintf(fid, ' %s \r\n', content{i,1});
                    elseif ismember(i,text_lines) %Printing text only lines
                        for j=1:size(content,2)
                            fprintf(fid, ' %14s', content{i,j});
                        end
                        fprintf(fid,'\r\n');
                    elseif ismember(i,dash_lines) %Printing dashed break lines
                        dash_length=num2str(15*loaddimFlag);
                        dash_format=" %"+dash_length+"s";
                        for j=1:1
                            fprintf(fid, char(dash_format), content{i,1});
                        end
                        fprintf(fid,'\r\n');
                    else %Printing numerical results in scientific notation
                        for j=1:size(content,2)
                            A_str = sprintf('% 10.6e',content{i,j});
                            exp_portion=extractAfter(A_str,'e');
                            num_portion=extractBefore(A_str,'e');
                            if strlength(exp_portion)<4
                                zeros_needed=4-strlength(exp_portion);
                                zeros_add{1}=repmat('0',1,zeros_needed);
                                exp_portion=[exp_portion(1),zeros_add{1},exp_portion(2:end)];
                                A_str=[num_portion,'e',exp_portion];
                            end
                            fprintf(fid,' %s',A_str);
                        end
                        fprintf(fid,'\r\n');
                    end
                    
                end
                fclose(fid);
                %Write filename to command window
                fprintf('\n'); fprintf(description); fprintf(' FILE: '); fprintf(filename); fprintf('\n');
            catch
                fprintf('\nUNABLE TO PRINT '); fprintf('%s %s', upper(description),'FILE. ');
            end
        end
    end
    
    if FLAGS.excel == 1
        %Output calibration load approximation
        if FLAGS.mode==1
            filename = 'CALIB ALG Tare Corrected Load Approximation.csv';
            description='CALIBRATION ALGEBRAIC MODEL LOAD APPROXIMATION';
        else
            filename = 'CALIB ALG Output Approximation.csv';
            description='CALIBRATION ALGEBRAIC MODEL OUTPUT APPROXIMATION';
        end
        approxinput=aprxIN;
        print_approxcsv(filename,approxinput,description,pointID,series,series2,loadlist,output_location);
    end
    
    if FLAGS.calib_model_save==1
        %Output all calibration parameters
        filename = ['AOX_CALIBRATION_MODEL_',REPORT_NO,'.mat'];
        fullpath=fullfile(output_location,filename);
        description='CALIBRATION MODEL';
        try
            model=FLAGS.model;
            save(fullpath,'coeff','ANOVA','loadlist','model');
            fprintf('\n'); fprintf(description); fprintf(' FILE: '); fprintf(filename); fprintf('\n');
        catch ME
            fprintf('\nUNABLE TO SAVE '); fprintf('%s %s', upper(description),'FILE. ');
            fprintf('\n')
        end
    end
end

%% GRBF Calibration Specific Outputs
if strcmp(section,{'Calibration GRBF'})==1
    Term_Names=customMatrix_labels(loadlist,voltagelist,voltdimFlag,loaddimFlag,FLAGS.model,'voltages',numBasis); %Get label names for custom equation matrix
    
    if FLAGS.excel == 1
        %Prints coefficients to csv file
        filename = 'AOX_ALG-GRBF_MODEL_COEFFICIENT_MATRIX.csv';
        description='CALIBRATION GRBF MODEL ALGEBRAIC COEFFICIENT MATRIX';
        print_coeff(filename,coeff_algRBFmodel,description,Term_Names,loadlist,output_location)
        
        
        %Output calibration load approximation
        if FLAGS.mode==1
            filename = 'CALIB GRBF Tare Corrected Load Approximation.csv';
            description='CALIBRATION ALGEBRAIC+GRBF MODEL LOAD APPROXIMATION';
        else
            filename = 'CALIB GRBF Output Approximation.csv';
            description='CALIBRATION ALGEBRAIC+GRBF MODEL OUTPUT APPROXIMATION';
        end
        approxinput=aprxINminTARE2;
        print_approxcsv(filename,approxinput,description,pointID,series,series2,loadlist,output_location);
        
        %Output GRBF Widths
        filename = 'AOX_GRBF_Epsilon.csv';
        input=epsHist;
        precision='%.16f';
        description='CALIBRATION GRBF EPSILON';
        print_dlmwrite(filename,input,precision,description,output_location);
        
        %Output GRBF coefficients
        filename = 'AOX_GRBF_Coefficients.csv';
        input=coeff_algRBFmodel_RBF;
        precision='%.16f';
        description='CALIBRATION GRBF COEFFICIENTS';
        print_dlmwrite(filename,input,precision,description,output_location);
        
        %Output GRBF centers
        filename = 'AOX_GRBF_Centers.csv';
        input=centerIndexHist;
        precision='%.16f';
        description='CALIBRATION GRBF CENTER INDICES';
        print_dlmwrite(filename,input,precision,description,output_location);
        
        %Output GRBF h value
        filename = 'AOX_GRBF_h.csv';
        input=h_GRBF;
        precision='%.16f';
        description='CALIBRATION GRBF "h"  VALUE';
        print_dlmwrite(filename,input,precision,description,output_location);
    end
    
    if FLAGS.calib_model_save==1
        %Output all calibration parameters
        filename = ['AOX_CALIBRATION_MODEL_',REPORT_NO,'.mat'];
        fullpath=fullfile(output_location,filename);
        description='CALIBRATION MODEL';
        try
            model=FLAGS.model;
            save(fullpath,'coeff','ANOVA','loadlist','model','coeff_algRBFmodel','epsHist','center_daHist','h_GRBF','ANOVA_GRBF');
            fprintf('\n'); fprintf(description); fprintf(' FILE: '); fprintf(filename); fprintf('\n');
        catch ME
            fprintf('\nUNABLE TO SAVE '); fprintf('%s %s', upper(description),'FILE. ');
            fprintf('\n')
        end
    end
    
    %%% ANOVA Stats AJM 6_12_19
    if FLAGS.anova==1
        anova_output(ANOVA_GRBF,nterms,numpts,loaddimFlag,'CALIB GRBF',output_location,Term_Names,loadlist,tR2,gee);
    end
    %%% ANOVA Stats AJM 6_8_19
end

%% Algebraic Validation Specific Outputs
if strcmp(section,{'Validation Algebraic'})==1
    %OUTPUTING APPROXIMATION WITH PI FILE
    if FLAGS.approx_and_PI_print==1
        if FLAGS.mode==1
            section='VALID ALG Tare Corrected Load';
        else
            section='VALID ALG Output';
        end
        load_and_PI_file_output(aprxINminTAREvalid,loadPI_valid,pointID,series,series2,loadlist,output_location,section)
        
        %OUTPUTING APPROXIMATION FILE
    elseif FLAGS.excel == 1
        if FLAGS.mode==1
            filename = 'VALID ALG Tare Corrected Load Approximation.csv';
            description='VALIDATION ALGEBRAIC MODEL LOAD APPROXIMATION';
        else
            filename = 'VALID ALG Output Approximation.csv';
            description='VALIDATION ALGEBRAIC MODEL OUTPUT APPROXIMATION';            
        end
        approxinput=aprxINminTAREvalid;
        print_approxcsv(filename,approxinput,description,pointID,series,series2,loadlist,output_location);
    end
end

%% GRBF Validation Specific Outputs
if strcmp(section,{'Validation GRBF'})==1
    %OUTPUTING APPROXIMATION WITH PI FILE
    if FLAGS.approx_and_PI_print==1
        if FLAGS.mode==1
            section='VALID GRBF Tare Corrected Load';
        else
            section='VALID GRBF Output';
        end
        load_and_PI_file_output(aprxINminTARE2valid,loadPI_valid_GRBF,pointID,series,series2,loadlist,output_location,section)
        
        
    elseif FLAGS.excel == 1
        %Output validation load approximation
        if FLAGS.mode==1
            filename = 'VALID GRBF Tare Corrected Load Approximation.csv';
            description='VALIDATION ALGEBRAIC+GRBF MODEL LOAD APPROXIMATION';
        else
            filename = 'VALID GRBF Output Approximation.csv';
            description='VALIDATION ALGEBRAIC+GRBF MODEL OUTPUT APPROXIMATION';            
        end

        approxinput=aprxINminTARE2valid;
        print_approxcsv(filename,approxinput,description,pointID,series,series2,loadlist,output_location);
        
    end
end
end

function [RECOMM_ALG_EQN]=anova_output(ANOVA,nterms,numpts,loaddimFlag,section,output_location,Term_Names,loadlist,tR2,gee)
%Function creates all the outputs for ANOVA results

%INPUTS:
%  ANOVA = Structure containing results from ANOVA analysis
%  nterms = Number of terms (predictor variables) in regression model
%  nseries0 = Number of series
%  numpts = Number of datapoints (observations)
%  loaddimFlag = Dimension of load data (# channels)
%  section = current section of code. Expected values: 'Calibration Algebraic', 'Calibration GRBF', 'Validation Algebraic', 'Validation GRBF'
%  output_location = Path for output files
%  Term_Names =  Term labels for predictor variables
%  loadlist = Labels for load variables
%  tR2 = target residuals squared
%  gee = mean of target residuals in each channel

%OUTPUTS:
%  RECOMM_ALG_EQN = Matrix for recommended regression model.  1's and 0's for which predictor variables to include in the model

%Output ANOVA results

totalnum = size(ANOVA(1).sig,1);
totalnumcoeffs = [1:totalnum];
totalnumcoeffs2 = [2:totalnum+1];
dsof = numpts-nterms-1;

loadstatlist = {'Load', 'Sum_Sqrs', 'PRESS_Stat', 'DOF', 'Mean_Sqrs', 'F_Value', 'P_Value', 'R_sq', 'Adj_R_sq', 'PRESS_R_sq'};
regresslist = {'Term_Index','Term_Name', 'Coeff_Value', 'CI_95cnt', 'T_Stat', 'P_Value', 'VIF_A', 'Signif'};

STAT_LOAD=cell(loaddimFlag,1);
REGRESS_COEFFS=cell(loaddimFlag,1);
for k=1:loaddimFlag
    RECOMM_ALG_EQN(:,k) = [1.0*ANOVA(k).sig([1:nterms])];
    manoa2(k,:) = [loadlist(k), tR2(1,k), ANOVA(k).PRESS, dsof, gee(1,k), ANOVA(k).F, ANOVA(k).p_F, ANOVA(k).R_sq, ANOVA(k).R_sq_adj, ANOVA(k).R_sq_p];
    ANOVA01(:,:) = [totalnumcoeffs; ANOVA(k).beta'; ANOVA(k).beta_CI'; ANOVA(k).T'; ANOVA(k).p_T'; ANOVA(k).VIF'; 1.0*ANOVA(k).sig']';
    ANOVA1_2(:,:) = num2cell([ANOVA01([1:nterms],:)]);
    STAT_LOAD{k} = array2table(manoa2(k,:),'VariableNames',loadstatlist(1:10));
    REGRESS_COEFFS{k} = cell2table([ANOVA1_2(:,1),Term_Names,ANOVA1_2(:,2:end)],'VariableNames',regresslist);
end

warning('off', 'MATLAB:xlswrite:AddSheet'); warning('off', 'MATLAB:DELETE:FileNotFound'); warning('off',  'MATLAB:DELETE:Permission')
filename = [section,' ANOVA STATS.xlsx'];
fullpath=fullfile(output_location,filename);

try
    delete(char(fullpath))
    for k=1:loaddimFlag
        writetable(STAT_LOAD{k},fullpath,'Sheet',k,'Range','A1');
        writetable(REGRESS_COEFFS{k},fullpath,'Sheet',k,'Range','A4');
    end
    fprintf(['\n',section,' METHOD ANOVA STATISTICS FILE: ']); fprintf(filename); fprintf('\n ');
catch ME
    fprintf(['\nUNABLE TO PRINT ',section,' METHOD ANOVA STATISTICS FILE. ']);
    if (strcmp(ME.identifier,'MATLAB:table:write:FileOpenInAnotherProcess'))
        fprintf('ENSURE "'); fprintf(char(filename)); fprintf('" IS NOT OPEN AND TRY AGAIN')
    end
    fprintf('\n')
end
warning('on',  'MATLAB:DELETE:Permission'); warning('on', 'MATLAB:xlswrite:AddSheet'); warning('on', 'MATLAB:DELETE:FileNotFound')
end


function []=print_coeff(filename,coeff_input,description,termList,loadlist,output_location)
%Function prints coefficients to csv file.  Error handling
%included to catch if the file is open and unable to write

%INPUTS:
%  filename = Filename for file to write
%  coeff_input = Matrix of coefficients to write to file
%  description = Description of file for Command Window output
%  termList = List of term labels for predictor variables
%  loadlist = Labels for load variables
%  output_location = Path for location of output files

%OUTPUTS:
%  []

fullpath=fullfile(output_location,filename);
try
    top_row=[{'Term Name'},loadlist]; %Top label row
    full_out=[top_row; termList, num2cell(coeff_input)]; %full output
    writetable(cell2table(full_out),fullpath,'writevariablenames',0); %write to csv
    fprintf('\n'); fprintf(description); fprintf(' FILE: '); fprintf(filename); fprintf('\n');
catch ME
    fprintf('\nUNABLE TO PRINT '); fprintf('%s %s', upper(description),'FILE. ');
    if (strcmp(ME.identifier,'MATLAB:table:write:FileOpenInAnotherProcess')) || (strcmp(ME.identifier,'MATLAB:table:write:FileOpenError'))
        fprintf('ENSURE "'); fprintf(char(filename));fprintf('" IS NOT OPEN AND TRY AGAIN')
    end
    fprintf('\n')
end
end