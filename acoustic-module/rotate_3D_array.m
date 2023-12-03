function A = rotate_3D_array(A,R)
% The function rotate_3D_array rotates a 3D array A with a rotation
% described by the rotation matrix R. R must be a rotation that can be
% expressed as a product of 90 degree rotations.
%
% INPUT:
% - A: 3D array to be rotated
% - R: rotation matrix
%
% OUTPUT:
% - B: rotated array
%
% Nathan Blanken, University of Twente, 2022

% Check dimensions of rotation matrix:
if ~(size(R,1)==3 && size(R,2)==3 && ismatrix(R))
    error('Rotation matrix must have dimensions 3x3.')
end

% Find the dimensions that need to be flipped:
[flip_dim,~] = find(transpose(R)==-1);
flip_dim = transpose(flip_dim);

% Execute the required flips:
for n = flip_dim
    A = flip(A,n);
end

% Find the permutation order:
[perm_order,~] = find(transpose(abs(R))==1);
perm_order = transpose(perm_order);

% Execute the permutation:
A = permute(A, perm_order);


end