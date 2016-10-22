function stats = calcColocStats(data, ROI, thresholds, displayOn)
 
% CALCCOLOCSTATS calculates the Manders (M1 and M2) and Pearson 
% coefficients in 3D for the specified ROI using the given threshold values.
% M1diff and M2diff are also calculated.
%
% INPUT data: 5D matrix containing the image data. Format - (X, Y, Z, T, C)           
%       ROI: 4D binary region of interest. Format - (X, Y, Z, T)
%       thresholds: vector containing threshold values
%       displayOn: Boolean, set to true to display message at start
%
% OUTPUT stats: struct containing the calculated colocalization statistics
%
% REMARKS: M1diff, M2diff are the difference between the Manders
%          coefficients and their expected value. For more details see 
%          Mcdonald et al. Journal of microscopy 252.3 (2013): 295-302.
%
%          The Pearson is coefficient is calculated using only voxels above 
%          the threshold values for both channels. For more details see 
%          Adler et al. Cytometry Part A 77.8 (2010): 733-742.
%
% created by: Jeremy Pike
% DATE: 19-Oct-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% If true display a message to let the user know this function has been
% called
if displayOn  
   disp('Calculating colocalization statistics...') 
end
 
% Find number of time-points and channels
[~, ~, ~, numTimePoints, numChannels] = size(data);
 
% Check there are two channels
if numChannels ~= 2
    disp('Please provide 2 channel data in the specified 5D format')
    return  
end
 
% Loop though time-points
for t=1:numTimePoints 
    %Pull data for given time-point
    dataTP = data(:,:,:,t,:);
    % Calculate M1, M2, M1diff and M2diff
    [stats.M1(t,1), stats.M2(t,1), stats.M1diff(t,1), stats.M2diff(t,1)] = manders(dataTP, ROI(:,:,:,t), thresholds);
    % Calculate the Pearson coefficient
    stats.Pearson(t,1) = pearson(dataTP, ROI(:,:,:,t), thresholds);   
end
 
end
 
function [M1, M2, M1diff, M2diff] = manders(data,ROI,thresholds)
 
% MANDERS calculates the Manders (M1 and M2) coefficients 
%
% INPUT data: 5D matrix containing the 2 channel image data. 
%             Format - (X, Y, Z, T, C)           
%       ROI: 4D binary region of interest. Format - (X, Y, Z, T)
%       thresholds: vector containing threshold values
%
% OUTPUT M1, M2: Manders coefficients
%        M1diff, M2diff: The Manders coefficients after subtracting the
%                        expected value assuming random signal distribution
%                        within the ROI. For more details see Mcdonald et
%                        al. Journal of microscopy 252.3(2013): 295-302.
%
% created by: Jeremy Pike
% DATE: 19-Oct-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Extract the data for each channel
dataC1=data(:,:,:,:,1);
dataC2=data(:,:,:,:,2);
% Reshape ROI as a vector
ROILine=ROI(:);
% Calculate the total volume (in voxels) of the ROI
totalVol = sum(ROILine);
% Reshape channel data as vectors
dataC1Line=dataC1(:); 
dataC2Line=dataC2(:);  
% Delete voxels outside ROI
dataC1Line(ROILine==0)=[];
dataC2Line(ROILine==0)=[];
% Set voxels less than corresponding threshold to zero
dataC1Line(dataC1Line<thresholds(1,1))=0; 
dataC2Line(dataC2Line<thresholds(2,1))=0; 
% Calculate the volume (in voxels) of isolated signal in each channel
volumeC1=sum(dataC1Line>0);
volumeC2=sum(dataC2Line>0);
% Find the voxels in channel 1 which colocalize with signal in channel 2
dataC1Overlap=dataC1Line;
dataC1Overlap(dataC2Line==0)=[];
% Find the voxels in channel 2 which colocalize with signal in channel 1
dataC2Overlap=dataC2Line;
dataC2Overlap(dataC1Line==0)=[];
 
% If there is no signal in channel 1 display a warning
if sum(dataC1Line)==0
   disp('Warning: cannot calculate valid coefficient, returning M1=0')
   % Set M1 to zero
   M1=0;
else
% Else calculate the first Manders coefficient    
    M1=sum(dataC1Overlap)/sum(dataC1Line);
end
% If there is no signal in channel 2 display a warning
if sum(dataC2Line)==0
   disp('Warning: cannot calculate valid coefficient, returning M2=0')
   % Set M2 to zero
   M2=0;
else
% Else calculate the second Manders coefficient
    M2=sum(dataC2Overlap)/sum(dataC2Line);
end
 
% Calculate M1diff by subtracting the fractional volume of channel 2
M1diff = M1 - volumeC2/totalVol;
% Calculate M2diff by subtracting the fractional volume of channel 1
M2diff = M2 - volumeC1/totalVol;
end
 
function R = pearson(data, ROI, thresholds)
 
% PEARSON calculates the Pearson coefficient using only voxels within the
% ROI with intensity above the specified thresholds for each channel
%
% INPUT data: 5D matrix containing the 2 channel image data. 
%             Format - (X, Y, Z, T, C)           
%       ROI: 4D binary region of interest. Format - (X, Y, Z, T)
%       thresholds: vector containing threshold values
%
% OUTPUT R: The Pearson coefficient
 
% REMARKS: For further justification of this approach see Adler et al. 
%          Cytometry Part A 77.8 (2010): 733-742.
%
% created by: Jeremy Pike
% DATE: 19-Oct-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
 
% Extract the data for each channel
dataC1=data(:,:,:,:,1);
dataC2=data(:,:,:,:,2);
% Reshape ROI as a vector
ROILine=ROI(:);
% Reshape channel data as vectors
dataC1Line=dataC1(:);  
dataC2Line=dataC2(:);
% Remove voxels with channel 1 intensity below channel 1 threshold from ROI
ROILine(dataC1Line<thresholds(1,1))=0;  
% Remove voxels with channel 2 intensity below channel 2 threshold from ROI
ROILine(dataC2Line<thresholds(2,1))=0;
% Isolate voxels in this new ROI
dataC1Line(ROILine==0)=[];
dataC2Line(ROILine==0)=[];
 
% If there is less than 2 isolated voxels display a warning
if isempty(dataC1Line)==1 ||  size(dataC1Line,1)==1 
    disp('Warning: cannot calculate valid coefficient, returning R=0')
    %Set Pearson coefficient to zero
    R = 0;
else
% Calculate the Pearson coefficient across isolated voxels
    meanC1=mean(dataC1Line);
    meanC2=mean(dataC2Line);
    R = (sum((dataC1Line-meanC1).*(dataC2Line-meanC2)))/(sqrt(sum((dataC1Line-meanC1).^2)*sum((dataC2Line-meanC2).^2)));
end
 
end


