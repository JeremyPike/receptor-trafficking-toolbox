function thresholds = otsuThresholds(data, ROI)
 
% otsuThresholds performs Otsu thresholding on each channel of the input
% data using only voxels within the specified ROI using Matlab’s multithresh
% function
%
% INPUT data: 5D matrix containing the image data. Format - (X, Y, Z, T, C)           
%       ROI: 4D binary region of interest. Format - (X, Y, Z, T)
%
% OUTPUT thresholds: vector containing threshold values
%
% created by: Jeremy Pike
% DATE: 19-Oct-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Find number of channels
numChannels = size(data,5);
 
% Create empty vector for threshold values
thresholds=zeros(numChannels,1);
 
disp('Computing Otsu thresholds...')
 
% Reshape ROI as vector
ROILine=ROI(:);
for c=1:numChannels
    % Pull out channel data
    channelData=data(:,:,:,:,c);
    % Reshape channel data as vector
    channelDataLine=channelData(:);
    % Delete all pixels outside ROI
    channelDataLine(ROILine==0)=[];
    % Calculate channel threshold (across all time-points) using the Otsu
    % approach
    thresholds(c,1)=multithresh(channelDataLine,1);
end
 
    
end
