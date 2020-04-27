function []=print_approxcsv(filename,approxinput,description,pointID,series1,series2,loadlist,output_location)
%Function prints load approximation to csv file.  Error handling
%included to catch if the file is open and unable to write

%INPUTS:
%  filename = Name for file to be saved
%  approxinput  =  Approximation matrix to be written to file
%  description = Description for Command Window output
%  pointID = Point ID labels for each datapoint, provided by input .csv
%  series1 = Series1 ID labels
%  series2 = Series2 ID labels
%  loadlist = Load label names for each channel
%  output_location = File output path

%OUTPUTS:
% []

fullpath=fullfile(output_location,filename); %Combine file name and output path for full path
try
    top_row=[{'Point ID','Series1','Series2'},loadlist]; %Top label row
    full_out=[top_row; pointID, num2cell(num2cell(series1)), num2cell(series2),num2cell(approxinput)]; %full output
    writetable(cell2table(full_out),fullpath,'writevariablenames',0); %write to csv
    fprintf('\n'); fprintf(description); fprintf(' FILE: '); fprintf(filename); fprintf('\n');
catch ME %Output error message if unable to output .csv
    fprintf('\nUNABLE TO PRINT '); fprintf('%s %s', upper(description),'FILE. ');
    if (strcmp(ME.identifier,'MATLAB:table:write:FileOpenInAnotherProcess')) || (strcmp(ME.identifier,'MATLAB:table:write:FileOpenError'))
        fprintf('ENSURE "'); fprintf(char(filename));fprintf('" IS NOT OPEN AND TRY AGAIN')
    end
    fprintf('\n')
end
end