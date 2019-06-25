%Function creates all the outputs for the validation, algebraic section
%This simplifies following the main code

function []=valid_alg_output(targetResvalid,loadCapacitiesvalid,FLAGS,fileNamevalid,numptsvalid,nseriesvalid,zapvalid,zapSTDEVvalid,loadlist,resSquarevalid,aprxINminGZvalid,seriesvalid,excessVecvalidkeep,dimFlag)

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
    alg_Tares_stdev_valid = array2table(zapSTDEVvalid,'VariableNames',loadlist(1:dimFlag));
    alg_Tares_stdev_valid= [series_table_valid, alg_Tares_stdev_valid]
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
    plotResPages(seriesvalid, targetResvalid, loadCapacitiesvalid, stdDevPercentCapacityvalid, loadlist)
    %    hold off
end

end