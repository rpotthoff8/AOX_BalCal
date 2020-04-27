function []=print_dlmwrite(filename,input,precision,description,output_location)
%Function prints provided input to file using dlmwrite.  Error handling
%included to catch if the file is open and unable to write

%INPUTS:
%  filename = Name for file to be saved
%  input  =  input to be written to file
%  precision = Precision preference for dmlwrite
%  description = Description for Command Window output
%  output_location = File output path

%OUTPUTS:
% []

fullpath=fullfile(output_location,filename); %Combine file name and output path for full path
try
    dlmwrite(fullpath,input,'precision',precision); %Output file
    fprintf('\n'); fprintf(description); fprintf(' FILE: '); fprintf(filename); fprintf('\n');
catch ME %If unable to write file, display error message
    fprintf('\nUNABLE TO PRINT '); fprintf('%s %s', upper(description),'FILE. ');
    if (strcmp(ME.identifier,'MATLAB:dlmwrite:FileOpenFailure'))
        fprintf('ENSURE "'); fprintf(char(filename)); fprintf('" IS NOT OPEN AND TRY AGAIN')
    end
    fprintf('\n')
end

end