%Function creates all the outputs for the calibration,  GRBF section
%This simplifies following the main code

function []=calib_GRBF_output(FLAGS,targetRes2,loadCapacities,series0,loadlist,aprxINminGZ2,wHist,cHist,centerIndexHist,numBasis,taresGRBF,taresGRBFSTDEV,dainputscalib,dimFlag,nseries0,resSquare2,numpts0)


for k2=1:length(targetRes2(1,:))
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
    
    twoSigmaGRBF = standardDev2'.*2;
    calib_GRBF_2Sigma = array2table(twoSigmaGRBF,'VariableNames',loadlist(1:dimFlag))
    
    %Should I use strtrim()  ? -AAM 042116
    taresGRBFactual = taresGRBF;
    series_table = table([1:nseries0]','VariableNames',{'SERIES'});
    calib_GRBF_Tares = array2table(taresGRBFactual,'VariableNames',loadlist(1:dimFlag));
    calib_GRBF_Tares = [series_table, calib_GRBF_Tares]
    calib_GRBF_STDEV_Tares = array2table(taresGRBFSTDEV,'VariableNames',loadlist(1:dimFlag));
    calib_GRBF_STDEV_Tares = [series_table, calib_GRBF_STDEV_Tares]
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