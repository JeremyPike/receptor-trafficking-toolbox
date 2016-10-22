function exampleVisualisations(data, ROI, thresholds, meta, timePoint, analysisChannel)
 
% EXAMPLEVISULAISATIONS creates 2 example visualisation figures using the 
% results of the workflow
%
% INPUT data: 5D matrix containing the image data. Format - (X, Y, Z, T, C)           
%       ROI: 4D binary region of interest. Format - (X, Y, Z, T)
%       thresholds: vector containing threshold values
%       meta: Metadata containing voxel dimensions
%       timepoint: Time-point to be visualised
%       analysisChannel: Channel for weighting of co-occurrence map
%
% REMARKS The 3D visualisation requires the vol3d function by Joe Conti and
%         Oliver Woodford available at, http://uk.mathworks.com/
%         matlabcentral/fileexchange/22940-vol3d-v2
 
%
% created by: Jeremy Pike
% DATE: 21-Oct-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
disp('Creating example visualisations...') 
 
% Retrieve necessary data dimensions
[rows, cols, numSlices, ~, numChannels] = size(data);
% Check there are two channels
if numChannels ~= 2
    disp('Please provide 2 channel data in the specified 5D format')
    return  
end
% Create 3D Cartesian grid for 3D visualisation
[X,Y,Z] = meshgrid(1:cols,1:rows,1:numSlices);
% Scale to physical distance using metadata
X = X * meta.voxelSizeX;
Y = Y * meta.voxelSizeY;
Z = Z * meta.voxelSizeZ;
 
% Calculate overlap of both channels after thresholding
dataC1 = data(:,:,:,timePoint,1);
dataC2 = data(:,:,:,timePoint,2);
overlap = (dataC1 >= thresholds(1,1)) .* (dataC2 >= thresholds(2,1));
% Calculate co-occurrence map weighted with the specified channel
cooccurrenceMap =  data(:,:,:,timePoint,analysisChannel) .* overlap .* ROI(:,:,:,timePoint);
%Normalise between zero and one
cooccurrenceMap = (cooccurrenceMap-min(cooccurrenceMap(:)))/(max(cooccurrenceMap(:))-min(cooccurrenceMap(:)));
 
%% 3D Visualisation of the co-occurrence map and cellular segmetation
figure('units','normalized','outerposition',[0 0 1 1])
% 3D rendering of co-occurrence map using vol3d
% http://uk.mathworks.com/matlabcentral/fileexchange/22940-vol3d-v2
vol3d('cdata',(cooccurrenceMap),'XData', [0 cols* meta.voxelSizeX],'YData', [0 (rows)*meta.voxelSizeY],'ZData', [0 (numSlices)*meta.voxelSizeZ] );
% change colour map
colormap(hot); 
% change transparency of map
alphamap('rampup'); 
% make patch object for cellular segmentation
p = patch(isosurface(X,Y,Z,ROI(:,:,:,timePoint),0)); 
% transparent faces and black edges for patch
set(p,'FaceColor','green','EdgeColor','none','FaceAlpha',0.1);
% set axis using data dimension
axis([0 cols 0 rows 0 numSlices]) 
axis equal 
% set view angle
view([-24 15])
% set background colour to black
set(gcf, 'Color', 'k');
% figure box on
box on
% axis not visible
axis off
 
%% Maximal projections of both channels and co-occurrence map
 
% Calculate maxiumal projection
dataC1MP=max(dataC1,[],3);
dataC2MP=max(dataC2,[],3);
cooccurrenceMapMP=max(cooccurrenceMap,[],3);
 
figure('units','normalized','outerposition',[0 0 1 1])
% Figure title
suptitle('Maximal projections')
% Subplot for channel 1
subplot(1,3,1)
imshow(dataC1MP,[])
title('Channel 1')
% Subplot for channel 2
subplot(1,3,2)
imshow(dataC2MP,[])
title('Channel 2')
% Subplot for co-occurrence map
ax3 = subplot(1,3,3);
imshow(cooccurrenceMapMP,[])
title(['Co-occurrence map (channel ' num2str(analysisChannel) ' weighted)'])
% Use a false colour look-up-table (Jet)
colormap(ax3, jet)
end
