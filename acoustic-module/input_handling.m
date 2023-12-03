function [Acquisition, SimulationParameters] = input_handling(...
    Acquisition, SimulationParameters, inputCell)
% Assign the optional input arguments collected in inputCell to the structs
% Acquisition and SimulationParameters.
%
% inputCell{1}: continue with the same k-Wave medium (boolean)
% inputCell{2}: frame number to continue from (integer)
% inputCell{3}: GPU device number
%
% Nathan Blanken, University of Twente, 2023

% Parameters for continuation of an interrupted simulation:
Acquisition.Continue = false; % Reuse the same medium
Acquisition.StartFrame = 1;   % Continue the simulation from this frame

if ~isempty(inputCell)
    if isempty(inputCell{1})
        % Ignore input argument
    elseif islogical(inputCell{1})
        Acquisition.Continue = inputCell{1};        
    else
        error('varargin{1} must be a boolean value.')
    end
end

if size(inputCell,2) > 1
    if isempty(inputCell{2})
        % Ignore input argument
    elseif isnumeric(inputCell{2}) && ...
            (floor(inputCell{2}) == inputCell{2}) && (inputCell{2} > 0)
        Acquisition.StartFrame = inputCell{2};
    else
        error('varargin{2} must be a positive integer.')
    end
end

if size(inputCell,2) > 2
    if isempty(inputCell{3})
        % Ignore input argument
    elseif isnumeric(inputCell{3}) && (floor(inputCell{3}) == inputCell{3})
        SimulationParameters.DeviceNumber = inputCell{3};
    else
        error('varargin{2} must be an integer.')
    end
end

end