%Read a VTU file and convert it to MATLAB format using the
%vtkToolbox written by Steffen Schuler, Karlsruhe Institute of Technology.
%(https://github.com/KIT-IBT/vtkToolbox).
%
% - filename: the name of the VTU file
% - savename: the name of the .mat file to store the output struct
%
% Nathan Blanken, University of Twente, 2023

addpath('vtkToolbox/MATLAB')

geometryFolder = '/renal_tree';

% Source VTU file:
filename = [geometryFolder filesep 'renal_tree.vtu'];

% Destination .MAT file (do not modify):
savename = [geometryFolder filesep 'vtu.mat'];

% Read VTU data
verbose = true; % Display progress messages
vtuStruct = readVTK(filename, verbose);

% Units of the vtu file:
vtuProperties.lengthUnit    = 1e-6; % [m];
vtuProperties.velocityUnit  = 1;    % [m/s]

% Normal to the inlet pointing inwards (currently supported only along one
% of the cartesian axes):
vtuProperties.inletNormal = [0 0 1];
vtuProperties.inletDiameter = 0.6e-3; % Maximum diameter of the inlet [m]

% Field of cellData that contains the velocities:
vtuProperties.velocityField = 'velocity_phy';

disp('Saving output struct...')
save(savename,'vtuStruct','vtuProperties','-v7.3')

rmpath('vtkToolbox/MATLAB')

