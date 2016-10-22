
function dataBS = backgroundSubtract3D(data, SEsize, voxelSizeX, voxelSizeZ)
 
% BACKGROUNDSUBTRACT3D using a spherical structural element to perform 3D 
% rolling ball background subtraction
%
% INPUT data:5D matrix containing the image data. Format - (X, Y, Z, T, C) 
%       SEsize: radius of structural element in physical distance units 
%               (eg microns)           
%       voxelSizeX, voxelSizeZ: Voxel spacing for axial and lateral
%                               dimensions. Should be in the same units as 
%                               SEsize  
%
% OUTPUT dataBS: 5D matrix containing the background subtracted image data. 
 
% created by: Jeremy Pike
% DATE: 19-Oct-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
disp('Performing background subtraction...')
 
% Find number of time-points and channels in data
[~,~,~,numTimePoints,numChannels]=size(data);
 
% Create spherical structural element with given radius
SE = sphericalSE(SEsize, voxelSizeX, voxelSizeZ);      
 
% Create empty matrix for background subtracted data
dataBS = zeros(size(data));
 
%Loop through time-points
for t=1:numTimePoints
    %Loop though channels
    for c=1:numChannels
        %Pull data for given time-point and channel
        dataTemp=data(:,:,:,t,c);
        %Use morphological opening to calculate the background
        background=imopen(dataTemp,SE);
        % Subtract the background
        dataBSTemp = dataTemp-background;
        %Set any voxels less than zero to zero
        dataBSTemp(dataBSTemp<0)=0;
        % Store result in dataBS
        dataBS(:,:,:,t,c)=dataBSTemp;
        
    end
end
           
 
          
