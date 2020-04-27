function [OUTLIER_ROWS,num_outliers,prcnt_outliers,rowOut,colOut]=ID_outliers(targetRes,numpts0,numSTD,FLAGS)
%function takes in the residual matrix, returns rows of outliers, number of
%outliers, and percent of datapoints that are outliers. Outliers are
%determined by the number of standard deviations away from the mean
%residual a point's residual is. It is possible (and likely) that a
%datapoint will be an outlier based based on the residual in 1 or more channels, but
%not all channels.  In these cases, the datapoints is still considered to
%be "bad" and thrown out as an outlier for calibration calculations in
%every channel.

%INPUTS:
%  targetRes = Matrix of residuals
%  numpts0 = Total number of datapoints (observations) in each channel
%  numSTD = User preference for number of standard deviations set as the outlier cutoff
%  FLAGS = Structure of flags determined by user inputs in GUI

%OUTPUTS:
%  OUTLIER_ROWS = Indices of Rows (observations) that have been flagged as outliers
%  num_outliers = Total number of outlier observations 
%  prcnt_outliers = Percentage of outliers in total number of observations
%  rowOut = Row (observation number) for outliers
%  colOut = Channel that denoted points as outliers

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