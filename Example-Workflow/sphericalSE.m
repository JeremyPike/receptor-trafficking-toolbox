function SE3D = sphericalSE(SESize, voxelSizeX, voxelSizeZ)
% SPHERICALSE creates a spherical structural element with specified radius.
% Accounts for variation in axial and lateral spacing
%
% INPUT SEsize: radius of structural element in physical distance units 
%               (eg microns)           
%       voxelSizeX, voxelSizeZ: Voxel spacing for axial and lateral
%                               dimensions. Should be in the same units as 
%                               SEsize  
%
% OUTPUT SE3D: Spherical structural element
 
% created by: Jeremy Pike
% DATE: 19-Oct-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    [xx,yy,zz] = ndgrid(-ceil(SESize/voxelSizeX):ceil(SESize/voxelSizeX),-ceil(SESize/voxelSizeX):ceil(SESize/voxelSizeX),-ceil(SESize/voxelSizeZ):ceil(SESize/voxelSizeZ));
    SE3D = sqrt((voxelSizeX*xx).^2 + (voxelSizeX*yy).^2 + (voxelSizeZ*zz).^2) <= SESize;
end

