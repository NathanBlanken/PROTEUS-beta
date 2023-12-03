function main_pressure_field(filename, savename)
% =========================================================================
% SIMULATE PRESSURE MAPS
% Compute the pressure in the lateral-axial plane and the elevation-axial 
% plane. Store minimum and maximum pressure values.
%
% input:    filename:    the file containing the simulation settings
%           savename:    the folder to save the pressure maps data
%
% Alina Kuliesh,  Delft University of Technology
% Nathan Blanken, University of Twente
% 2022
% =========================================================================

load(filename,'Geometry', 'Medium', 'SimulationParameters', ...
    'Transducer', 'Transmit')

% simulation settings
run_param = sim_setup(SimulationParameters);

% Location of the geometry data:
Geometry.GeometriesPath = run_param.GeometriesPath;

savedir = [run_param.savedir filesep savename];

% Create save directory if nonexistent:
if ~isfolder(savedir)
    mkdir(savedir)
end
disp(['Pressure maps will be saved in: ' newline savedir '.' newline])

% Define the k-Wave grid:
[kgrid, Grid] = define_grid(SimulationParameters, Geometry);

run_param.PML = Grid.PML;

% Define the k-Wave medium:
[medium, vessel_grid] = define_medium(Grid, Medium, Geometry);

save([savedir '/medium.mat'],'medium','vessel_grid','Grid')

% record signals long enough for back and forth pass of the wave
run_param = compute_travel_times(run_param, ...
    Geometry,Medium,Transducer,Transmit);

% create the time array
kgrid.Nt = floor(run_param.tr(1) / kgrid.dt) + 1;

% Filter and resample transmit signal:
Transmit = preprocess_transmit(Transmit,Medium,kgrid);

% Distribute integration points at the transducer surface:
Transducer = get_transducer_integration_points(...
    Transducer, Transmit, Medium, Grid);

% define the transducer source
disp('Creating k-Wave sensor object for transducer.')
[sensor_transducer, sensor_weights] = define_sensor_transducer(...
    Transducer, Grid);

mask_idx_trans = find(logical(sensor_transducer.mask));

% Simulate the field with all elements on:
Transmit.SeqPulse = 'full';

disp('Creating k-Wave source object for transducer.')
source_transducer = define_source_transducer(Transducer, Transmit, ...
    Medium, Grid, transpose(sensor_weights), mask_idx_trans);

[~,iy] = min(abs(Grid.y));
[~,iz] = min(abs(Grid.z));

% DEFINE THE SENSOR FOR THE PRESSURE MAPS
% XY plane sensor:
sensor.mask = zeros(Grid.Nx,Grid.Ny,Grid.Nz,'logical');
sensor.mask(:,:,iz) = 1;
sensor_idx_xy = find(sensor.mask);

% XZ plane sensor:
sensor.mask = zeros(Grid.Nx,Grid.Ny,Grid.Nz,'logical');
sensor.mask(:,iy,:) = 1;
sensor_idx_xz = find(sensor.mask);

% Add XY plane to XZ plane sensor:
sensor.mask(:,:,iz) = 1;
sensor_idx_all = find(sensor.mask);

sensor.record={'p_max','p_min'};

sensor_data = run_simulation(run_param, kgrid, medium, ...
    source_transducer, sensor);

% Extract the pressure maps from the sensor data:

[~,sensor_data_xy_idx,~] = intersect(sensor_idx_all, sensor_idx_xy);
[~,sensor_data_xz_idx,~] = intersect(sensor_idx_all, sensor_idx_xz);

sensor_data_xy.p_max = sensor_data.p_max(sensor_data_xy_idx);
sensor_data_xz.p_max = sensor_data.p_max(sensor_data_xz_idx);

sensor_data_xy.p_min = sensor_data.p_min(sensor_data_xy_idx);
sensor_data_xz.p_min = sensor_data.p_min(sensor_data_xz_idx);

sensor_data_xy.p_max = reshape(sensor_data_xy.p_max,Grid.Nx,Grid.Ny);
sensor_data_xz.p_max = reshape(sensor_data_xz.p_max,Grid.Nx,Grid.Nz);

sensor_data_xy.p_min = reshape(sensor_data_xy.p_min,Grid.Nx,Grid.Ny);
sensor_data_xz.p_min = reshape(sensor_data_xz.p_min,Grid.Nx,Grid.Nz);

file_name = 'pressure_maps.mat';
save([savedir filesep file_name],'sensor_data_xy','sensor_data_xz','Grid')

end