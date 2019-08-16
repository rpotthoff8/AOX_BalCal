function [OUTLIER_ROWS,num_outliers,prcnt_outliers,rowOut,colOut]=ID_outliers(targetRes,loadCapacities,numpts0,dimFlag,numSTD,FLAGS)
%function takes in the residual matrix, returns rows of outliers, number of
%outliers, and percent of datapoints that are outliers

% Use the modeled input for the rest of the calculations
targetRes_mean = mean(targetRes);

% Identify outliers. They are considered outliers if the residual
% is more than numSTD standard deviations from the mean.
targetRes_std = std(targetRes);

targetRes_norm = (targetRes-targetRes_mean)./targetRes_std;


[rowOut,colOut]=find(abs(targetRes_norm) > numSTD); %Find row and column indices of outliers

OUTLIER_ROWS = unique(rowOut);

num_outliers = length(OUTLIER_ROWS);
prcnt_outliers = 100.0*num_outliers/numpts0;

if FLAGS.print == 1
    %% Identify the Possible Outliers
    fprintf(' \n***** \n');
    fprintf('Number of Outliers =');
    fprintf(string(num_outliers));
    fprintf('\n%s','Outliers % of Data =');
    fprintf(string(prcnt_outliers));
    fprintf('\n ');
end

end