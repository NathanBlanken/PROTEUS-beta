%==========================================================================
% Backpropagating N streamlines to densely populate the inlet surface.
%
% The step tolerance parameter can be adjusted for maximum precision, but
% the voxelization of the simulation becomes the bottle neck for low
% values.
%==========================================================================

clear
clc
close all

% Use parallel computing for the streamline backpropagation:
useparfor = false;

%==========================================================================
% SETTINGS
%==========================================================================

% Folder containing the geometry data:
geometryFolder = ['..' filesep 'geometry_data'];
[filename, pathname] = uigetfile([geometryFolder filesep '*vtu.mat'], ...
    'Select vtu.mat file');

% Number of streamlines to propagate back to the inlet:
Nstreamlines = 10e3;

% Percentile to take for reference velocity:
velocityPercentile = 0.95;

% Maximum number of grid points the ODE solver can skip over in one 
% integration step:
stepTolerance = 5; % leniency parameter

% Maximum integration time:
frameRate = 250; % [Hz]
numberOfFrames = 5e3;
Tmax = (numberOfFrames-1)/frameRate;

skip = 25; % Skip this number of points for plotting vessel

% The maximum distance a streamline end point is allowed to be from any
% other end points is expressed as a fraction of the inlet diameter. This
% parameter is used to filter out isolated end points which will not be on
% the inlet surface. A lower factor should be chosen for a higher number of
% streamlines.
proximityFactor = 0.1;

%==========================================================================
% READ VESSEL DATA
%==========================================================================

% MATLAB file with VTU data of the flow simulation:
filepath = fullfile(pathname, filename);

[vtuStruct, Grid] = load_vessel_data(filepath);
load(filepath,'vtuProperties')

% Plot the vessel:
figure()
plot3(vtuStruct.points(1:skip:end,1),...
      vtuStruct.points(1:skip:end,2),...
      vtuStruct.points(1:skip:end,3),'.');
xlabel('X (m)')
ylabel('Z (m)')
zlabel('Z (m)')
title('Vessel mesh vertices')

%--------------------------------------------------------------------------
% ODE solver options
%--------------------------------------------------------------------------

maxStep = get_step_size(vtuStruct, velocityPercentile, stepTolerance);
options.MaxStep = maxStep;
save(fullfile(pathname, 'ode_options.mat'),'options') % Do not modify
options = odeset(options,'Events',@(t,y)exitVesselFcn(t,y,Grid));

vtuStruct.velocities = -vtuStruct.velocities; % Backpropagate

% Function handle to the ODE:
odefun = @(t,y) transpose(...
    get_velocity(transpose(y), Grid, vtuStruct.velocities));

%==========================================================================
% BACKPROPAGATE STREAMLINES TOWARDS THE INLET
%==========================================================================

points_start = zeros(Nstreamlines,3); % Array for streamline start points.
points_end   = zeros(Nstreamlines,3); % Array for streamline end points.
t_end        = zeros(Nstreamlines,1); % Array for streamline end times.

t1 = tic;

if useparfor == true
    
    %----------------------------------------------------------------------
    % PARALLEL COMPUTING OF STREAMLINES
    %----------------------------------------------------------------------
    
    point_start_cell = cell(1,Nstreamlines);
    points_end_cell  = cell(1,Nstreamlines);
    t_end_cell       = cell(1,Nstreamlines);

    parfor k = 1:Nstreamlines

        disp(['Computing streamline ' num2str(k) ' of ' ...
            num2str(Nstreamlines) '.'])

        % Position the bubble in the bulk of the vessel:
        startPosition = draw_start_position(1, vtuStruct);
        tspan = [0 Tmax];

        %------------------------------------------------------------------
        % COMPUTE STREAMLINE
        %------------------------------------------------------------------

        [t,positions] = ode23(odefun, tspan, startPosition(:),options);

        % Store the end point and end time of the streamline:
        point_start_cell{k} = startPosition(:);
        points_end_cell{k} = positions(end,:);
        t_end_cell{k} = t(end);

    end

    % Assign the values in the cells to the matrices:
    for k = 1:Nstreamlines
        points_start(k,:) = point_start_cell{k};
        points_end(k,:)   = points_end_cell{k};
        t_end(k) = t_end_cell{k};
    end
        
else
    
    %----------------------------------------------------------------------
    % SERIAL COMPUTING OF STREAMLINES
    %----------------------------------------------------------------------
    
    figure()

    for k = 1:Nstreamlines

        disp(['Computing streamline ' num2str(k) ' of ' ...
            num2str(Nstreamlines) '.'])

        % Position the bubble in the bulk of the vessel:
        startPosition = draw_start_position(1, vtuStruct);
        tspan = [0 Tmax];

        %------------------------------------------------------------------
        % COMPUTE STREAMLINE
        %------------------------------------------------------------------

        [t,positions] = ode23(odefun, tspan, startPosition(:),options);

        %------------------------------------------------------------------
        % PLOT STREAMLINE
        %------------------------------------------------------------------
        plot3(positions(:,1),positions(:,2),positions(:,3));
        xlabel('X (m)')
        ylabel('Y (m)')
        zlabel('Z (m)')
        title('Streamlines')
        hold on
        drawnow

        % Store the end point and end time of the streamline:
        points_start(k,:) = startPosition(:);
        points_end(k,:)   = positions(end,:);
        t_end(k) = t(end);

    end  

end

toc(t1)
hold off

%--------------------------------------------------------------------------
% PLOT STREAMLINE START POINTS
%--------------------------------------------------------------------------
figure()
plot3(points_start(:,1),points_start(:,2),points_start(:,3),'r.');
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
title('Streamline start points')
hold on

%==========================================================================
% FILTER INLET POINTS
%==========================================================================

%--------------------------------------------------------------------------
% Only keep end points of streamlines that have terminated before the
% maximum integration time. These streamlines are most likely to have
% reached the inlet. Streamlines will terminate just outside the inlet. 
% Move points back into the vessel by half a grid spacing.
%--------------------------------------------------------------------------
points = points_end;
points = points(t_end<Tmax,:);
points = points + 1/2*vtuStruct.cellsize.*vtuProperties.inletNormal;

%--------------------------------------------------------------------------
% Remove points not inside the vessel.
%--------------------------------------------------------------------------
vtuInd = get_vtu_indices(points,Grid);
points(vtuInd==0,:) = [];

%--------------------------------------------------------------------------
% Remove isolated points and remove points that are too far from the
% centroid of all other points.
%--------------------------------------------------------------------------
maxDistance = vtuProperties.inletDiameter*proximityFactor;
points = filter_points_proximity(points, maxDistance);
distance = vecnorm((points - mean(points)),2,2);
points(distance > 2*vtuProperties.inletDiameter,:) = [];

%--------------------------------------------------------------------------
% Determine the normal to the inlet plane. The normal must be aligned with
% one of the cartesian axes.
%--------------------------------------------------------------------------

[~, normal_axis] = min(std(points));

% Check if the determined inlet normal corresponds to the inlet normal in
% the vtu properties struct:
if isfield(vtuProperties,'inletNormal') && ...
        (normal_axis ~= find(vtuProperties.inletNormal))
    warning('Inlet normal discrepancy. Overwriting new inlet normal.')
end

vtuProperties.inletNormal = zeros(1,3);
vtuProperties.inletNormal(normal_axis) = 1;

save(filepath,'vtuProperties',"-append")

%--------------------------------------------------------------------------
% Remove points too far from the inlet plane
%--------------------------------------------------------------------------
% Round points to nearest grid point along the normal axis:
vtuInd = get_vtu_indices(points,Grid);
points(:,normal_axis) = vtuStruct.points(vtuInd,normal_axis);

% Estimate the coordinate of the inlet plane based on the most frequently
% occuring coordinate along the normal axis:
inlet_coordinate = mode(points(:,normal_axis));

% Compute the distance of each point to the inlet plane:
distance = abs(points(:,normal_axis) - inlet_coordinate);

% Remove points to far from the inlet plane:
points(distance>(vtuStruct.cellsize(normal_axis)/2),:) = [];

%--------------------------------------------------------------------------
% Remove points too far from the centroid of the inlet
%--------------------------------------------------------------------------

% Remove points that are in the infinite inlet plane but not in the inlet
% (subset of the inlet plane):
distance = vecnorm((points - mean(points)),2,2);
points(distance > vtuProperties.inletDiameter,:) = []; 
 

%==========================================================================
% PLOT AND SAVE RESULTS
%==========================================================================
figure();
plot3(points(:,1),points(:,2),points(:,3),'.');
xlabel('X (m)')
ylabel('Z (m)')
zlabel('Z (m)')
title('Streamline end points')

save(fullfile(pathname, 'backpropagation_points.mat'),'points');