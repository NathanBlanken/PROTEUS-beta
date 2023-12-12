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

if Geometry.EmbedVessel
    % Read STL data:
    V = READ_stl(Geometry.STLfile);

    % Translate and rotate STL vertex coordinates:
    V = V*Geometry.STLunit;
    V = V - transpose(Geometry.BoundingBox.Center);
    V = rotate_stl(V,Geometry.Rotation);
    V = V + transpose(Geometry.Center);

    % convert STL mesh to voxels:
    vessel_grid = VOXELISE(Grid.x, Grid.y, Grid.z, V);
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

function meshXYZ = rotate_stl(meshXYZ,R)
%ROTATE_STL rotates STL coordinates meshXYZ with rotation matrix R.
%
% The N-by-3-by-3 array meshXYZ is defined as the output of READ_stl.m in
% the mesh voxelisation toolbox by Adam H. Aitkenhead:
%  1 row for each facet
%  3 cols for the x,y,z coordinates
%  3 pages for the three vertices
%
% The 3-by-3 rotation matrix is defined such that B = R*A, where A is a
% column vector and B is the rotated column vector.

for i = 1:3
    V = transpose(meshXYZ(:,:,i)); % 3-by-N matrix of vertex coordinates
    V = R*V;                       % Rotate the vertex coordinates
    meshXYZ(:,:,i) = transpose(V); % Assign result
end

end