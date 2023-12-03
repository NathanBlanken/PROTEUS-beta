function [Parameters,status] = load_from_file()
% Load .mat file with simulation parameters into the GUI.

filter = '.mat';
defname = 'GUI_output_parameters.mat';
[filename,filepath] = uigetfile(filter,[],defname);

if filename == 0
    Parameters = 0;
    status = 'notloaded';
    return
    
else
    status = 'loaded';
end

try 
    Parameters = load(fullfile(filepath,filename));
catch
    Parameters = 0;
    status = 'invalid';
end