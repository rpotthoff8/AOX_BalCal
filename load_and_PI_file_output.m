function []=load_and_PI_file_output(aprxINminGZ,loadPI,pointID,series1,series2,loadlist,output_location,section)
%Functio outputs a .xlsx file containing the load approximations and
%prediction intervals

%INPUTS:
%  aprxINminGZ = Load approximation matrix
%  loadPI = Matrix of Prediction interval (+/- value) for each load approximation
%  pointID = Point ID for each approximation, as determined by input file
%  series1 = Series1 labels for each point
%  series2 = Series2 labels for each point
%  loadlist = Labels for loads
%  output_location = Path for file output
%  section = Character array for current section, used for file labeling.  Expected possibilities: VALID ALG, VALID GRBF, APPROX ALG, APPROX GRBF

%OUTPUTS:
% []

warning('off', 'MATLAB:xlswrite:AddSheet'); warning('off', 'MATLAB:DELETE:FileNotFound'); warning('off',  'MATLAB:DELETE:Permission') %Surpress warnings
description=[upper(section), ' APPROXIMATION WITH PREDICTION INTERVAL']; %Text description for command window output
try %Check current section and determine filename accordingly
    if contains(section,'VALID')==1 %If in validation section
        filename=[section,' Tare Corrected Load Approximation w PI.xlsx'];
    elseif contains(section,'APPROX')==1 %If in approximation section
        filename=[section, ' Global Load Approximation w PI.xlsx'];
    end
    fullpath=fullfile(output_location,filename); %Full path for file output
    
    top_row=[{'Point ID','Series1','Series2'},loadlist]; %Top label row
    full_out{1}=[top_row; pointID, num2cell(num2cell(series1)), num2cell(series2),num2cell(aprxINminGZ)]; %full output
    full_out{2}=[top_row; pointID, num2cell(num2cell(series1)), num2cell(series2),num2cell(cellstr(string(aprxINminGZ)+' +/- '+string(loadPI)))]; %full output
    full_out{3}=[top_row; pointID, num2cell(num2cell(series1)), num2cell(series2),num2cell(loadPI)]; %full output
    
    for i=1:3
        writetable(cell2table(full_out{i}),fullpath,'writevariablenames',0,'Sheet',i,'UseExcel', false); %write to xlsx
    end
    fprintf('\n'); fprintf(description); fprintf(' FILE: '); fprintf(filename); fprintf('\n'); %Command Window output
    xlsx=1; %Tracking variable that xlsx write was successful
    
catch ME %Catch if unable to write to xlsx
    xlsx=0; %Tracking variable that xlsx write was not successful
    fprintf('\nUNABLE TO PRINT '); fprintf('%s %s', upper(description),'FILE. '); %Output to command window
    if (strcmp(ME.identifier,'MATLAB:table:write:FileOpenInAnotherProcess')) || (strcmp(ME.identifier,'MATLAB:table:write:FileOpenError')) %Check error
        fprintf('ENSURE "'); fprintf(char(filename));fprintf('" IS NOT OPEN AND TRY AGAIN. ') %Output for error cause
    end
    fprintf('\n ATTEMPTING TO PRINT APPROXIMATION ONLY FILE.'); %text output
    fprintf('\n')
    
    %If unable to print .xslx file containing load approximation and PI,
    %attempt to print only approximation .csv file
    
    %Determine filename
    if contains(section,'VALID')==1 %If in validation section
        filename=[section,' Tare Corrected Load Approximation.csv'];
    elseif contains(section,'APPROX')==1 %If in approximation section
        filename=[section, ' Global Load Approximation.csv'];
    end
    
    approxinput=aprxINminGZ; %Load approximation output
    description=[upper(section),' ALGEBRAIC MODEL GLOBAL LOAD APPROXIMATION']; %Description for command window output
    print_approxcsv(filename,approxinput,description,pointID,series1,series2,loadlist,output_location); %Write approximation output file
end


if xlsx==1 %If xlsx file was written

    try %Rename excel sheets and delete extra sheets, only possible on PC
        [~,sheets]=xlsfinfo(fullpath);
        s = what;
        e = actxserver('Excel.Application'); % # open Activex server
        e.DisplayAlerts = false;
        e.Visible=false;
        ewb = e.Workbooks.Open(char(fullpath)); % # open file (enter full path!)
        
        if max(size(sheets))>3
            %cycle through, deleting all sheets other than the 1st 3 sheets
            for i=4:max(size(sheets))
                ewb.Sheets.Item(i).Delete;  %Delete sheets
            end
        end
        
        ewb.Worksheets.Item(1).Name = 'Load Approximation'; %rename 1st sheet
        ewb.Worksheets.Item(2).Name = 'Load Approximation w_PI'; %rename 2nd sheet
        ewb.Worksheets.Item(3).Name = 'Prediction Interval'; %rename 3rd sheet
        
        ewb.Save % # save to the same file
        ewb.Close
        e.Quit
        delete(e);
    end
    
    
end
warning('on',  'MATLAB:DELETE:Permission'); warning('on', 'MATLAB:xlswrite:AddSheet'); warning('on', 'MATLAB:DELETE:FileNotFound') %Reset warning states
end