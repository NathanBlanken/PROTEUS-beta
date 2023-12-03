%==========================================================================
% Backpropagating 2500 streamlines to densely populate the inlet surface.
%
% The step tolerance parameter is set to 5. It can be adjusted for maximum 
% precision, but the voxelization of the simulation becomes the bottle neck
% for low values.
%==========================================================================

clear
clc
close all

%==========================================================================
% SETTINGS
%==========================================================================

% Folder containing the geometry data:
geometryFolder = '../geometry_data/renal_tree';

% Number of streamlines to propagate back to the inlet:
Nstreamlines = 2500;

% Percentile to take for reference velocity:
velocityPercentile = 0.95;

% Maximum number of grid points the ODE solver can skip over in one 
% integration step:
stepTolerance = 5; % leniency parameter

% Maximum integration time:
frameRate = 500; % [Hz]
numberOfFrames = 1e4;
Tmax = (numberOfFrames-1)/frameRate;

skip = 100; % Skip this number of points for plotting vessel

%==========================================================================
% READ VESSEL DATA
%==========================================================================

% MATLAB file with VTU data of the flow simulation:
filename = [geometryFolder filesep 'vtu.mat'];

[vtuStruct, Grid] = load_vessel_data(filename);
load(filename,'vtuProperties')

% Axis perpendicular to the inlet:
normal_axis = find(vtuProperties.inletNormal);

% Plot the vessel:
figure()
plot3(vtuStruct.points(1:skip:end,1),...
      vtuStruct.points(1:skip:end,2),...
      vtuStruct.points(1:skip:end,3),'.');
xlabel('X (m)')
ylabel('Z (m)')
zlabel('Z (m)')

%--------------------------------------------------------------------------
% ODE solver options
%--------------------------------------------------------------------------

maxStep = get_step_size(vtuStruct, velocityPercentile, stepTolerance);
options.MaxStep = maxStep;
save([geometryFolder filesep 'ode_options.mat'],'options') % Do not modify
options = odeset(options,'Events',@(t,y)exitVesselFcn(t,y,Grid));

vtuStruct.velocities = -vtuStruct.velocities; % Backpropagate

% Function handle to the ODE:
odefun = @(t,y) transpose(...
    get_velocity(transpose(y), Grid, vtuStruct.velocities));

%==========================================================================
% BACKPROPAGATE STREAMLINES TOWARDS THE INLET
%==========================================================================

points = zeros(Nstreamlines,3); % Array for holding streamline end points.
t_end  = zeros(Nstreamlines,1); % Array for holding streamline end times.

t1 = tic;
figure()

for k = 1:Nstreamlines
    
    disp(['Computing streamline ' num2str(k) ' of ' ...
        num2str(Nstreamlines) '.'])
       
    % Position the bubble in the bulk of the vessel:
    startPosition = draw_start_position(1, vtuStruct);
    
    tspan = [0 Tmax];
    
    %----------------------------------------------------------------------
    % COMPUTE STREAMLINE
    %----------------------------------------------------------------------

    [t,positions] = ode23(odefun, tspan, startPosition(:),options);
    
    %----------------------------------------------------------------------
    % PLOT STREAMLINE
    %----------------------------------------------------------------------
    plot3(positions(:,1),positions(:,2),positions(:,3));
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')
    hold on
    drawnow
    
    % Store the end point and end time of the streamline:
    points(k,:) = positions(end,:);
    t_end(k) = t(end);
              
end

toc(t1)
hold off


%==========================================================================
% FILTER INLET POINTS
%==========================================================================

%--------------------------------------------------------------------------
% Only keep end points of streamlines that have terminated before the
% maximum integration time. These streamlines are most likely to have
% reached the inlet. Streamlines will terminate just outside the inlet. 
% Move points back into the vessel by half a grid spacing.
%--------------------------------------------------------------------------
points = points(t_end<Tmax,:);
points = points + 1/2*vtuStruct.cellsize.*vtuProperties.inletNormal;

%--------------------------------------------------------------------------
% Remove points not inside the vessel. Round points to nearest grid point
% along the normal axis.
%--------------------------------------------------------------------------
vtuInd = get_vtu_indices(points,Grid);
points(vtuInd==0,:) = [];
vtuInd(vtuInd==0)   = [];
points(:,normal_axis) = vtuStruct.points(vtuInd,normal_axis);

%--------------------------------------------------------------------------
% Remove points too far from the inlet plane
%--------------------------------------------------------------------------

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

streamlineProperties

save([geometryFolder filesep 'backpropagation_points.mat'],'points');