%Function creates all the outputs for the calibration,  GRBF section
%This simplifies following the main code

function []= valid_GRBF_output(FLAGS,targetResvalid,targetRes2valid,loadCapacitiesvalid,loadlist,dimFlag,numBasis,validSeries,aprxINminGZ2valid,nseriesvalid,taresGRBFvalid,taresGRBFSTDEVvalid,resSquare2valid,numptsvalid)

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
    
    twoSigmaGRBFvalid = standardDev2valid'.*2;
    GRBF_2Sigmavalid = array2table(twoSigmaGRBFvalid,'VariableNames',loadlist(1:dimFlag))
    
    %Should I use strtrim()  ? -AAM 042116
    series_table_valid = table([1:nseriesvalid]','VariableNames',{'SERIES'});
    GRBF_Taresvalid = array2table(taresGRBFvalid,'VariableNames',loadlist(1:dimFlag));
    GRBF_Taresvalid = [series_table_valid, GRBF_Taresvalid]
    GRBF_TaresSTDEVvalid = array2table(taresGRBFSTDEVvalid,'VariableNames',loadlist(1:dimFlag));
    GRBF_TaresSTDEVvalid = [series_table_valid, GRBF_TaresSTDEVvalid]
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
    plotResPages(validSeries, targetRes2valid, loadCapacitiesvalid, stdDevPercentCapacity2valid, loadlist)
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