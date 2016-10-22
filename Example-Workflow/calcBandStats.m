function stats = calcBandStats(data, ROI, thresholds, meta, params)
 
% CALCBANDSTATS splits the ROI into banded volumes based on distance from
% the ROI edge. Each band has a fixed width in microns. The integrated 
% signal and colocalization (co-occurrence and correlation) is then 
% calculated for each band. For more details see the manuscript by Pike et
% al.
%
% INPUT data: 5D matrix containing the image data. Format - (X, Y, Z, T, C)           
%       ROI: 4D binary region of interest. Format - (X, Y, Z, T)
%       thresholds: vector containing threshold values
%       meta: Metadata containing voxel dimensions
%       bandWidth: Width of each band in microns
%
% OUTPUT stats: struct containing the calculated statistics
%
% REMARKS: Requires bwdistsc.m to perform the 3D distance transform. This is
%          available at  https://uk.mathworks.com/matlabcentral/
%          fileexchange/15455-3d-euclidean-distance-transform-for-variable-
%          data-aspect-ratio and described by Y. Mishchenko (2015) Signal, 
%          Image and Video Processing, 9(1), 19-27.
%
% created by: Jeremy Pike
% DATE: 19-Oct-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
disp('Calculating band based statistics...') 
 
 
% Find number of time-points and channels
[~, ~, numSlices, numTimePoints, numChannels] = size(data);
 
% Check there are two channels
if numChannels ~= 2
    disp('Please provide 2 channel data in the specified 5D format')
    return  
end
 
% Loop though time-points
for t=1:numTimePoints 
   
    %Pull ROI for given time-point
    ROITP = ROI(:,:,:,t);
    %Calculate the 3D weighted distance transform using the ROI
    distTrans = bwdistsc(~ROITP,[meta.voxelSizeX, meta.voxelSizeY, meta.voxelSizeZ]);
    % maximum distance from segmentation edge
    maxDist = max(distTrans(:));
    % Calculate number of bands for given cell (round down)
    numBands=floor(maxDist/params.bandWidth);
    % Calculate the total volume (in voxels) of the ROI
    totalVol = sum(ROITP(:));
    %Pull data for given time-point and analysis channel
    dataTemp = data(:,:,:,t,params.analysisChannel);
    % Zero pixels outside ROI
    dataTemp(ROITP == 0) = 0;
    %Calcualte total signal for this channel
    totalSignal = sum(dataTemp(:));
    % Calculate overlap of both channels after thresholding
    dataC1 = data(:,:,:,t,1);
    dataC2 = data(:,:,:,t,2);
    overlap = (dataC1 >= thresholds(1,1)) .* (dataC2 >= thresholds(2,1));
    % Calculate co-occurrence map for specified channel
    cooccurrenceMap =  dataTemp .* overlap;
    %Calcualte total signal from the map
    totalMap = sum(cooccurrenceMap(:));
    % Loop through bands
    for n = 1:numBands
        
       % Define distance bounds for band
       lowerBound = (n - 1) * params.bandWidth;
       upperBound = n * params.bandWidth;
       % Create ROI for the band
       ROIBand = ones(size(ROITP));
       ROIBand(distTrans > upperBound) = 0;
       ROIBand(distTrans <= lowerBound) = 0;
       % Calculate band volume
       bandVol = sum(ROIBand(:));
       
       %Pull data for given time-point and analysis channel
       dataTemp = data(:,:,:,t,params.analysisChannel);
       % Zero pixels outside ROI
       dataTemp(ROIBand == 0) = 0;
       % get co-occurrence map 
       mapTemp = cooccurrenceMap;
       % Zero pixels outside ROI
       mapTemp(ROIBand == 0) = 0;
       
       % Calculate percentage of signal for this band
       stats.signalPercentage{t,1}(n,1) = sum(dataTemp(:)) / totalSignal * 100;
       % Calculate percentage of signal after volume correction
       % The expected value (fractional volume) of the band is subtracted
       stats.signalPercentageVolumeCor{t,1}(n,1) = stats.signalPercentage{t,1}(n,1) - bandVol / totalVol * 100;
       
       % Calculate percentage of co co-occurrence map for this band
       stats.cooccurrenceMapPercentage{t,1}(n,1) = sum(mapTemp(:)) / totalMap * 100;
       % Calculate percentage of co-occurrence map after volume correction
       stats.cooccurrenceMapPercentageVolumeCor{t,1}(n,1) = stats.cooccurrenceMapPercentage{t,1}(n,1) - bandVol / totalVol * 100;
       
       % Calculate colocalization measures for band
       bandColocStats = calcColocStats(data(:,:,:,t,:), ROIBand , thresholds, false);
       % Store in stats
       stats.M1{t,1}(n,1) = bandColocStats.M1;
       stats.M2{t,1}(n,1) = bandColocStats.M2;
       stats.Pearson{t,1}(n,1) = bandColocStats.Pearson;
       stats.M1diff{t,1}(n,1) = bandColocStats.M1diff;
       stats.M2diff{t,1}(n,1) = bandColocStats.M2diff;
 
    end
    
end
 
end
 

