% Standalone AOX_approx program
% Copyright 2019 Andrew Meade, Ali Arya Mokhtarzadeh, Javier Villarreal and John Potthoff.  All Rights Reserved.
%
% Required files to run:
%   AOX_Approx_GUI.m
%   AOX_Approx_GUI.fig
%   AOX_approx_funct.m
%   balCal_algEqns.m
%   calc_alg_PI.m
%   print_approxcsv.m
%   place_GRBF.m

%%
%initialize the workspace
clc;
clearvars;
close all;
fprintf('Copyright 2019 Andrew Meade, Ali Arya Mokhtarzadeh, Javier Villarreal, and John Potthoff.  All Rights Reserved.\n')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       USER INPUT SECTION
%Run GUI
out = AOX_Approx_GUI;
if out.cancel == 1
    return
end

%DEFINE THE PRODUCTION CSV INPUT FILE AND SELECT THE RANGE OF DATA VALUES TO READ
load(out.savePathapp,'-mat');
%     Load Calibration Model
load(out.calib_model_path);


%Load Output Variables
FLAGS.balCal = out.grbf; %TO SELECT Algebraic Model: set FLAGS.balCal = 1;  %TO SELECT Algebraic and GRBF Model: set FLAGS.balCal = 2;
FLAGS.loadPI = out.loadPI;
FLAGS.excel = out.excel;    %TO SAVE DATA TO CSV: set FLAGS.excel = 1;
FLAGS.approx_and_PI_print=out.approx_and_PI_print;
FLAGS.PI_print=out.PI_print;
FLAGS.model=model;

PI_pct=out.PI_percent_confidence;

%Set Save Location for Files
if FLAGS.excel ==1 || FLAGS.approx_and_PI_print==1 || FLAGS.PI_print==1
    FLAGS.save_files=1;
else
    FLAGS.save_files=0;
end

REPORT_NO=datestr(now,'yyyy-mmdd-HHMMSS');
output_location=[out.output_location,filesep];
if out.subfolder_FLAG==1 && FLAGS.save_files==1
    try
        new_subpath=fullfile(output_location,['AOX_Approx_Results_',REPORT_NO]);
        mkdir(char(new_subpath));
        output_location=new_subpath;
    catch
        fprintf('Unable to create new subfolder. Saving results in: '); fprintf('%s',output_location); fprintf('\n');
    end
end

if FLAGS.balCal == 2 %If RBFs were placed, put parameters in structure
    GRBF.wHist=wHist;
    GRBF.cHist=cHist;
    GRBF.center_daHist=center_daHist;
    GRBF.h=h_GRBF;
else
    GRBF='GRBFS NOT PLACED';
end

%% Approximation Calculations Section
%Function that performs all ANOVA calculations and outputs
[aprxINminGZapprox,loadPI_approx]=AOX_approx_funct(coeff,natzerosapprox,excessVecapprox,FLAGS,seriesapprox,series2approx,pointIDapprox,loadlist,output_location,GRBF,ANOVA,PI_pct);


fprintf('\n  ');
fprintf('\nCalculations Complete.\n');
if FLAGS.save_files==1
    fprintf('%s',strcat('Check '," ",output_location,' for output files.'))
end
fprintf('\n');
if isdeployed % Optional, use if you want the non-deployed version to exit immediately
  input('Press enter to finish and close');
end
