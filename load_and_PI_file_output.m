function []=load_and_PI_file_output(aprxINminGZ,loadPI,pointID,series1,series2,loadlist,output_location,section)
description=[upper(section), ' APPROXIMATION WITH PREDICTION INTERVAL'];
try
    filename=[upper(section),'_GLOBAL_ALG_APPROX_w_PI.xlsx'];
    fullpath=fullfile(output_location,filename);
    
    top_row=[{'Point ID','Series1','Series2'},loadlist]; %Top label row
    full_out{1}=[top_row; pointID, num2cell(num2cell(series1)), num2cell(series2),num2cell(aprxINminGZ)]; %full output
    full_out{2}=[top_row; pointID, num2cell(num2cell(series1)), num2cell(series2),num2cell(cellstr(string(aprxINminGZ)+' +/- '+string(loadPI)))]; %full output
    full_out{3}=[top_row; pointID, num2cell(num2cell(series1)), num2cell(series2),num2cell(loadPI)]; %full output
    
    for i=1:3
        writetable(cell2table(full_out{i}),fullpath,'writevariablenames',0,'Sheet',i,'UseExcel', false); %write to xlsx
    end
    fprintf('\n'); fprintf(description); fprintf(' FILE: '); fprintf(filename); fprintf('\n');
    xlsx=1;
    
catch ME
    xlsx=0;
    fprintf('\nUNABLE TO PRINT '); fprintf('%s %s', upper(description),'FILE. ');
    if (strcmp(ME.identifier,'MATLAB:table:write:FileOpenInAnotherProcess')) || (strcmp(ME.identifier,'MATLAB:table:write:FileOpenError'))
        fprintf('ENSURE "'); fprintf(char(filename));fprintf('" IS NOT OPEN AND TRY AGAIN. ')
    end
    fprintf('\n ATTEMPTING TO PRINT APPROXIMATION ONLY FILE.');
    fprintf('\n')
    
    filename = [upper(section),'_GLOBAL_ALG_APPROX.csv'];
    approxinput=aprxINminGZ;
    description=[upper(section),' ALGEBRAIC MODEL GLOBAL LOAD APPROXIMATION'];
    print_approxcsv(filename,approxinput,description,pointID,series1,series2,loadlist,output_location);
    
    
end


if xlsx==1
    warning('off', 'MATLAB:xlswrite:AddSheet'); warning('off', 'MATLAB:DELETE:FileNotFound'); warning('off',  'MATLAB:DELETE:Permission')
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
        
        ewb.Worksheets.Item(1).Name = 'GLOBAL_ALG_APPROX'; %rename 1st sheet
        ewb.Worksheets.Item(2).Name = 'GLOBAL_ALG_APPROX_w_PI'; %rename 2nd sheet
        ewb.Worksheets.Item(3).Name = 'ALG_PREDICTION_INTERVAL'; %rename 3rd sheet
        
        ewb.Save % # save to the same file
        ewb.Close
        e.Quit
        delete(e);
    end
    warning('on',  'MATLAB:DELETE:Permission'); warning('on', 'MATLAB:xlswrite:AddSheet'); warning('on', 'MATLAB:DELETE:FileNotFound')
    
end

end