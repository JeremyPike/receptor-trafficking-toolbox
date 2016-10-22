function segmentation = cellularSegment(data, meta, params)
 
% CELLULARSEGMENT performs cellular segmentation using a 3D implementation of
% DRLSE as described by Li et al. IEEE transactions on image processing 
% 19.12 (2010): 3243-3254.
%
% INPUT data: 4D matrix containing the image data. Format - (X, Y, Z, T)           
%       meta: Metadata containing voxel dimensions
%       params: Segmentation params as specified in exampleWorkflow.m
%
% OUTPUT segmentation: 4D binary matrix containing the segmentation result. 
%               Format - (X, Y, Z, T)
%
% REMARKS: A K-means based segmentation estimate is used to initialize the 
% level set function for the DRLSE
%
% created by: Jeremy Pike
% DATE: 15-Oct-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%% Calculate initial segmentation estimate using a K-means approach
 
disp('Computing initial segmentation estimate...')
 
[rows, cols, numSlices, numTimePoints] = size(data);
 
%Create spherical structual element
SE = sphericalSE(params.SEradius, meta.voxelSizeX, meta.voxelSizeZ);
 
%Create empty matrix for estimate
initialEstimate = zeros(size(data));
 
%Loop through time-points
for t = 1:numTimePoints
  
    %Data for given time-point
    dataTP = data(:,:,:,t);
    %Perform K-means clustering of data using specified number of clusters
    [kMeansInd, clust] = kmeans(dataTP(:), params.numClusters);
    % Combine the three clusters with the highest-means
    kMeansComb = zeros(size(kMeansInd));
    [~, orderInd]=sort(clust,'descend');
    for i=1:length(clust)-1;
        kMeansComb(kMeansInd==orderInd(i))=1;       
    end   
 
    % Reshape the vector to the size of the original data
    estimateTP = reshape(kMeansComb, [rows, cols, numSlices]); 
    
    %Clear all but the largest connected component (volume)
    lab = bwlabeln(estimateTP);
    stats = regionprops(lab, 'Area');
    centroids = cat(1, stats.Area);
    [~, ind]=max(centroids);
    estimateTP(lab ~= ind)=0;
    
    % Dilate, fill and then erode with spherical structural element
    estimateTP = imdilate(estimateTP, SE);   
    for z=1:numSlices
        estimateTP(:,:,z) = imfill(estimateTP(:,:,z));
    end
    estimateTP = imerode(estimateTP, SE);
    
    % Store in 4D matrix
    initialEstimate(:,:,:,t) = estimateTP;   
    
    disp(['Time-point ' num2str(t) ' of '  num2str(numTimePoints) ' complete'])
end
 
%% 3D DRLSE (level set segmentation)
 
disp('Computing level set segmentation...')
 
%Establish ratio between axial and lateral dimensions
deltaXY = 1;
deltaZ = meta.voxelSizeZ/ meta.voxelSizeX;
 
% Calculate the maximum weighting of distance regularization term using 
% the Courant–Friedrichs–Lewy (CFL) condition (see Appendix A. of Li et al.
% IEEE transactions on image processing 19.12 (2010): 3243-3254)
lh = 1/(deltaXY^2)+1/(deltaXY^2)+1/(deltaZ^2);
rh = 1/(2*params.timeStep);
 
%Subtract a small value from the bound;
mu = rh/lh - 0.05;
 
%Create empty matrix for segmentation result
segmentation = zeros(size(data));
 
%Loop though all time-points
for t=1:numTimePoints
    
    %Data for given time-point
    dataTP = data(:,:,:,t);
    % Calculate first order derivatives of data (scaled)
    [Ix,Iy,Iz]=gradient(dataTP, deltaXY, deltaXY, deltaZ);    
    %Calculate absolute value of gradient
    grad=sqrt(Ix.^2+Iy.^2+Iz.^2);    
    % Calculate edge indicator function g
    g=1./(1+grad.^params.gradValue);    
    %Scale between 0 and 1
    g=(g-min(g(:)))/(max(g(:))-min(g(:)));   
    % Calculate first order derivatives of g
    [vx, vy, vz]=gradient(g);   
    % Scale between -1 and +1
    vx=scaleNegPos(vx);
    vy=scaleNegPos(vy);
    vz=scaleNegPos(vz);   
    % Scale with ratio between axial and lateral voxel dimensions
    vx=vx/deltaXY;
    vy=vy/deltaXY;
    vz=vz/deltaZ;
    % Initialize phi as binary step function using initial segmentation
    % estimate
    phi=params.c0*ones(rows, cols, numSlices);
    phi(initialEstimate(:,:,:,t)==1)=-params.c0;
    
   %Variable for storing change volume contained with zero level set at the
   % previous iteration 
   volumePrev=0;
   % Evolve the level set function with a maximum number of iterations
   for n=1:params.iterMax 
        
        % Update phi using 3D edge based DRLSE
        phi=drlse_edge_3D(phi, g,vx,vy,vz, params.lambda, mu, params.alfa, params.epsilon, params.timeStep, deltaXY,deltaZ);
        
        %Construct binary matrix containing the current segmentation
        binary=zeros(rows, cols, numSlices);
        binary(phi<0)=1;   
        %Calculate volume of current segmentation
        volume=sum(binary(:));
        % Calculate change in volume from previous iteration
        volumeDiff=abs(volumePrev-volume);
        
        % If greater than minimum number of iterations and percentage
        % change in volume less than stopping condition for two sequential
        % iterations stop the evolution
        if n>params.iterMin && volumeDiff<(volume*params.stopPercent/100) && volumeDiffPrev<(volume*params.stopPercent/100) 
            break
        end
        % Update previous volume and previous volume difference
        volumePrev=volume;
        volumeDiffPrev=volumeDiff;
   end
   
   % Store final segmentation result for time-point
   
     %Clear all but the largest connected component (volume)
    lab = bwlabeln(binary);
    stats = regionprops(lab, 'Area');
    centroids = cat(1, stats.Area);
    [~, ind]=max(centroids);
    binary(lab ~= ind)=0;
    %Fill holes in each slice
    for z=1:numSlices
        binary(:,:,z) = imfill(binary(:,:,z));
    end
   segmentation(:,:,:,t)=binary;
   
   disp(['Time-point ' num2str(t) ' of '  num2str(numTimePoints) ' complete'])
end
 
 
end
 
 
function vScale=scaleNegPos(v)
%Scale vx between -1 and 1. Scale positive and negative components
%separately
    vpos=v;
    vpos(v<0)=0;   
    vpos=(vpos-min(vpos(:)))/(max(vpos(:))-min(vpos(:)));
    vneg=v;
    vneg(v>0)=0;
    vneg=-vneg; 
    vneg=(vneg-min(vneg(:)))/(max(vneg(:))-min(vneg(:))); 
    vScale=vpos-vneg;
end
