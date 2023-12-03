function [outStruct, Grid] = load_vessel_data(filename)
%LOAD_VESSEL_DATA Read VTK data represented as a MATLAB struct (as defined
%in vtkToolbox/README.md, https://github.com/KIT-IBT/vtkToolbox). Remove
%zero-velocity cells, convert to standard units, and set up a cartesian
%grid.
%
% INPUT:
% - filename:  the name of the VTU file
% - vtuUnit:   length unit of the VTU file in meters
%
% OUTPUT:
% - outstruct: the VTU file as a MATLAB struct with fields:
%              - points: the cell centres.
%              - velocities: the cell velocities.
%
% - Grid:      grid properties struct with fields:
%              - X,  Y,  Z:   array of grid coordinates in each dimension
%              - dX, dY, dZ:  grid spacing in each dimension
%              - NX, NY, NZ:  grid size in each dimension
%              - vtu_indices: sparse array with the mapping of each grid
%                             point to the corresponding entry in the list
%                             of points in the vtu file 
%                             ((NX x NY x NZ) x 3 array)
%
% Nathan Blanken, University of Twente, 2023

% Read the VTK data represented as a MATLAB struct:
disp('Loading VTK struct...')
load(filename, 'vtuStruct','vtuProperties');

% VTU file units for length and velocity:
xUnit = vtuProperties.lengthUnit;
vUnit = vtuProperties.velocityUnit;

% Get the velocity data:
outStruct.velocities = vtuStruct.cellData.(vtuProperties.velocityField);

% Removing cells with zero velocity:
disp('Removing cells with zero velocity norm...')
zeroVelocityIdx = logical(vecnorm(outStruct.velocities,2,2));
vtuStruct.cells = vtuStruct.cells(zeroVelocityIdx,:);
outStruct.velocities = outStruct.velocities(zeroVelocityIdx,:);

% Get the first vertex of each cell:
cellVertices = vtuStruct.points(vtuStruct.cells(:,1),:);

% Get unique grid coordinates:
disp('Retrieving cartesian grid...')
X = unique(cellVertices(:,1)); dX = mean(diff(X)); NX = length(X);
Y = unique(cellVertices(:,2)); dY = mean(diff(Y)); NY = length(Y);
Z = unique(cellVertices(:,3)); dZ = mean(diff(Z)); NZ = length(Z);

% Translate cell vertices to cell centers:
X = X + dX/2; Y = Y + dY/2; Z = Z + dZ/2;
cellsize = [dX dY dZ];
cellCenters = cellVertices + cellsize/2;

% Get the grid subscripts for the cells:
I = round((cellCenters(:,1) - X(1))/dX) + 1;
J = round((cellCenters(:,2) - Y(1))/dY) + 1;
K = round((cellCenters(:,3) - Z(1))/dZ) + 1;

% Convert the subscripts to a linear index:
ind = sub2ind([NX NY NZ], I, J, K);

% Create a matrix of size [NX*NY*NZ, 3] that holds the index in the vtu 
% velocity list for each grid point:
disp('Creating sparse matrix...')
Grid.vtu_indices = sparse(NX*NY*NZ,1);
Grid.vtu_indices = sparse(ind,ones(size(ind)),1:length(ind),NX*NY*NZ,1);

% Convert grid properties to standard units:
Grid.X = transpose(X)*xUnit; Grid.dX = dX*xUnit; Grid.NX = NX;
Grid.Y = transpose(Y)*xUnit; Grid.dY = dY*xUnit; Grid.NY = NY;
Grid.Z = transpose(Z)*xUnit; Grid.dZ = dZ*xUnit; Grid.NZ = NZ;

outStruct.velocities = outStruct.velocities*vUnit;
outStruct.points     = cellCenters*xUnit;
outStruct.cellsize   = cellsize*xUnit;

end