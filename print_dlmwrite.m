%Function prints provided input to file using dlmwrite.  Error handling
%included to catch if the file is open and unable to write

function []=print_dlmwrite(filename,input,precision,description,output_location)
fullpath=fullfile(output_location,filename);
try
    dlmwrite(fullpath,input,'precision',precision);
    fprintf('\n'); fprintf(description); fprintf(' FILE: '); fprintf(filename); fprintf('\n');
catch ME
    fprintf('\nUNABLE TO PRINT '); fprintf('%s %s', upper(description),'FILE. ');
    if (strcmp(ME.identifier,'MATLAB:dlmwrite:FileOpenFailure'))
        fprintf('ENSURE "'); fprintf(char(filename)); fprintf('" IS NOT OPEN AND TRY AGAIN')
    end
    fprintf('\n')
end

end