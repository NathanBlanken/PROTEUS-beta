function run_param = sim_setup(SimulationParameters)

% Full path to current function:
currentFile = mfilename('fullpath');

% Full path to acoustic module:
acousticModulePath = fileparts(currentFile);

% Full path to main directory:
startDirectory = fileparts(acousticModulePath);

addpath(acousticModulePath)
addpath(startDirectory)

% Get paths to the required folders:
PATHS = path_setup(startDirectory);

switch SimulationParameters.Solver
    
    case '3D'
        run_param.solver       = 'kspaceFirstOrder3D';
        run_param.DATA_CAST    = 'single'; % run locally and record movie
        run_param.record_movie = true;
    case '3DC'
        run_param.solver       = 'kspaceFirstOrder3DC';
        run_param.DATA_CAST    = 'single';
        run_param.record_movie = false;
        run_param.DATA_PATH    = PATHS.DataPath;
        run_param.BINARY_PATH  = PATHS.BinaryPath;
    case '3DG'
        run_param.solver       = 'kspaceFirstOrder3DG';
        run_param.DATA_CAST    = 'gpuArray-single';
        run_param.record_movie = false;
        run_param.DEVICE_NUM   = SimulationParameters.DeviceNumber;
        run_param.DATA_PATH    = PATHS.DataPath;
        run_param.BINARY_PATH  = PATHS.BinaryPath;
end


% Add toolbox paths:
addpath(PATHS.VoxelisationPath);
addpath(PATHS.kWavePath)

% Folder to save the simulation output:
run_param.savedir = PATHS.ResultsPath;

% Folder containing the vessel geometries:
run_param.GeometriesPath = PATHS.GeometriesPath;

% Folder containing the microbubble frames:
run_param.GroundTruthPath = PATHS.GroundTruthPath;

% Folder containing the microbubble simulation module:
run_param.MicrobubblePath = PATHS.MicrobubblePath;

run_param.N_interactions = SimulationParameters.NumberOfInteractions;

end