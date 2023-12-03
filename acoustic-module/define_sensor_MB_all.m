function [sensor, MB_idx_all, max_mb] = define_sensor_MB_all(...
    Grid, folder, frame_finish, N_sequence, Geometry)
% =========================================================================
% DEFINE THE SENSOR: MBs and transducer record pressure
% input:    kgrid
%           grid_coordinates
%           folder 
%           frame
%           Geometry
% output:   sensor - mask of sensors
%           MB - struct with linear indexes of MBs in the sensor mask and 
%           indexes of recorded pressure lines (corresponding to MBs) in 
%           sensor_data
%           transducer - update of the transducer struct with same indexes
%           as MBs
% 
%
% =========================================================================

sensor.mask = zeros(Grid.Nx, Grid.Ny, Grid.Nz);
max_mb = 1;

for frame = 1 : frame_finish
       
    for pulse_seq_idx = 1:N_sequence
    
        MB = load_microbubbles(folder, frame, pulse_seq_idx, Geometry, frame_finish);

        % Put the microbubbles on the grid:
        [MB.points, ~, MB_idx, ~] = voxelize_media_points(MB.points, Grid);
        
        if size(MB.points, 1) > max_mb
            max_mb = size(MB.points, 1);
        end

        % Put sensor at the microbubbles
        mask_only = true;
        [sensor,~] = update_sensor(sensor, MB.points, MB_idx, ...
            Grid, mask_only);
    
    end

    
end

sensor.record={'p'};
MB_idx_all = find(sensor.mask == 1);

end