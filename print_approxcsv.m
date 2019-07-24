%Function prints load approximation to csv file.  Error handling
%included to catch if the file is open and unable to write

function []=print_approxcsv(filename,approxinput,description,pointID,series1,series2,loadlist)

try
    top_row=[{'Point ID','Series1','Series2'},loadlist]; %Top label row
    full_out=[top_row; pointID, num2cell(num2cell(series1)), num2cell(series2),num2cell(approxinput)]; %full output
    writetable(cell2table(full_out),filename,'writevariablenames',0); %write to csv
    fprintf('\n'); fprintf(description); fprintf(' FILE: '); fprintf(filename); fprintf('\n');
catch ME
    fprintf('\nUNABLE TO PRINT '); fprintf('%s %s', upper(description),'FILE. ');
    if (strcmp(ME.identifier,'MATLAB:table:write:FileOpenInAnotherProcess')) || (strcmp(ME.identifier,'MATLAB:table:write:FileOpenError'))
        fprintf('ENSURE "'); fprintf(char(filename));fprintf('" IS NOT OPEN AND TRY AGAIN')
    end
    fprintf('\n')
end
end