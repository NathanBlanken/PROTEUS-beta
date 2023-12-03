function [medium, vessel_grid] = define_medium(Grid, Medium, Geometry)
% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% input:    kgrid
%           Medium
%           Geometry
% output:   medium - struct with speed of sound, density, B/A, alpha
%           coefficient maps. Maps consist of masks which correspond to 
%           different tissue types
%           vessel - struct with mask of the vessel branch geometry
% 
% vessel mask has properties of the blood, the rest of the medium has 
% tissue average properties.
% =========================================================================

% Grid dimensions
Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz;

%==========================================================================
% Medium inhomegeneity
%==========================================================================

% Normal distribution truncated at +/-cutoff standard deviations:
cutoff = Medium.InhomogeneityCutoff;
pd = truncate(makedist('Normal'),-cutoff,cutoff);

% Make the speed of sound and density inhomogeneous to generate a linear
% scatterer background:
background_map   = 1 + Medium.Inhomogeneity*...
                   random(pd,[Nx, Ny, Nz]);
               
%==========================================================================
% Background tissue maps
%==========================================================================
medium.sound_speed = Medium.SpeedOfSound * background_map;  %[m/s]
medium.density     = Medium.Density      * background_map;  %[kg/m^3]

medium.BonA        = Medium.BonA         * ones(Nx,Ny,Nz);

% power law absorption prefactor [dB/(MHz^y cm)]:
medium.alpha_coeff = Medium.AttenuationA * ones(Nx,Ny,Nz);
medium.alpha_power = Medium.AttenuationB;

%==========================================================================
% Dispersion
%==========================================================================
% The dispersion term in k-Wave is derived via the Kramers-Kronig
% relations. There exists a singularity for alpha_power = 1 (see page 33 of
% the k-Wave manual, Manual Version 1.1,  for further details). Do not
% simulate dispersion for powers close to 1.
if abs(medium.alpha_power - 1) < 0.1
    medium.alpha_mode =  'no_dispersion';
end

%==========================================================================
% Vessel
%==========================================================================

% Location of the STL file:
Geometry.STLfile = [Geometry.GeometriesPath filesep Geometry.Folder ...
    filesep Geometry.STLfile];

stl_coordinates = get_stl_coordinates(Grid,Geometry);

if Geometry.EmbedVessel
    % convert stl mesh to voxels
    vessel_grid = VOXELISE(stl_coordinates.X, stl_coordinates.Y, ...
        stl_coordinates.Z, Geometry.STLfile);

    % Rotate the vessel grid: 
    vessel_grid = rotate_3D_array(vessel_grid, Geometry.Rotation);
else
    vessel_grid = zeros(Nx, Ny, Nz,'logical');
end

% Assign the tissue properties of the vessel:
medium.sound_speed(vessel_grid) = Medium.Vessel.SpeedOfSound;
medium.density(vessel_grid)     = Medium.Vessel.Density;
medium.BonA(vessel_grid)        = Medium.Vessel.BonA; 
    
end


%==========================================================================
% FUNCTIONS
%==========================================================================

function [stl_coordinates] = get_stl_coordinates(Grid,Geometry)
% Transform the coordinates of the simulation grid to the coordinate system
% of the original STL file.
%
% INPUT:
% - grid_coordinates: the coordinates of the simulation grid
% - Geometry:         struct containing the orientation and position of the
%                     vessel tree.
% OUTPUT:
% - stl_coordinates:  coordinates of the simulation grid in the coordinate
%                     system of the STL file

BB = Geometry.BoundingBox;          % Bounding box of the STL file

% Centre the grid:
X = Grid.x - Geometry.Center(1); % [m]
Y = Grid.y - Geometry.Center(2); % [m]
Z = Grid.z - Geometry.Center(3); % [m]

% Rotate the grid:
R = transpose(Geometry.Rotation);
[X,Y,Z] = rotate_coordinate_vectors(X,Y,Z,R);

% Translate to the centre of the STL bounding box and convert to the units
% of the STL file:
stl_coordinates.X = (X + BB.Center(1))/Geometry.STLunit;
stl_coordinates.Y = (Y + BB.Center(2))/Geometry.STLunit;
stl_coordinates.Z = (Z + BB.Center(3))/Geometry.STLunit;

end


function [X,Y,Z] = rotate_coordinate_vectors(X,Y,Z,R)
% INPUT: 
% - X,Y,Z: coordinate vectors (row vectors) spanning a cartesian grid
% - R:     rotation matrix
%
% OUTPUT
% - X,Y,Z: coordinate vectors (row vectors) of the cartesian grid after
%          application of the rotation matrix R

% Matrix with extreme points (corners) of grid coordinates as column 
% vectors:
V  = [min(X) min(X) min(X) min(X) max(X) max(X) max(X) max(X); ...
      max(Y) max(Y) min(Y) min(Y) max(Y) max(Y) min(Y) min(Y); ...;
      max(Z) min(Z) min(Z) max(Z) max(Z) min(Z) min(Z) max(Z)];
 
% Column vector with coordinate vector length in each dimension:
N = [length(X); length(Y); length(Z)];

% Column vector with coordinate increment in each dimension:
dx = (max(X) - min(X))/(length(X)-1);
dy = (max(Y) - min(Y))/(length(Y)-1);
dz = (max(Z) - min(Z))/(length(Z)-1);
d  = [dx; dy; dz];

% Rotate:
V = R*V;
N = R*N;
d = R*d;

NX = abs(N(1)); NY = abs(N(2)); NZ = abs(N(3));
dx = abs(d(1)); dy = abs(d(2)); dz = abs(d(3));

% Coordinate vectors of rotated grid:
X = (0:(NX-1))*dx + min(V(1,:));
Y = (0:(NY-1))*dy + min(V(2,:));
Z = (0:(NZ-1))*dz + min(V(3,:));

end