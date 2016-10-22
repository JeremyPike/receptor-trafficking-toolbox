
function denoisedData = denoise(data, params)
 
% DENOISE denoises data using the PureDenoise plugin which implements the 
% PURE-LET scheme for removing of Poisson noise introduced by
% F. Luisier et al., Signal Processing, vol. 90, no. 2, 415-427, Feb 2010
%
% INPUT data: 5D matrix containing the raw image data
%             Format - (X, Y, Z, T, C)
%       params: denoising parameters as specified in exampleWorkflow.m
%
% OUTPUT data: 5D matrix containing the denoised data
%
% REMARKS: 
% - Requires ImageJ (https://imagej.nih.gov/ij/) with the PureDenoise plugin
%   installed (available at http://bigwww.epfl.ch/algorithms/denoise/)
% - Uses MIJ to run ImageJ within Matlab as described by Daniel Sage et al.
%   ImageJ User & Developer Conference,24-26 Oct 2012
%   (available at http://bigwww.epfl.ch/sage/soft/mij/)
%
% created by: Jeremy Pike
% DATE: 16-Oct-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
disp('De-noising data...');
 
% Get current working directory
directory=pwd;
% Add Mij.jar to the java class path
javaaddpath(params.MIJpath);
% Add ij.jar to the java class path
javaaddpath([params.IJpath '\ij.jar']);
% Start ImageJ using MIJ
MIJ.start(params.IJpath)
%Check the working directory hasn’t moved
cd(directory)
%Create parameter string for PureDenoise plugin
 paramStr=['parameters=''' num2str(params.numFrames) ' ' num2str(params.numSpinCycles) ''' estimation=''Auto Global'' '];
 
% Retrieve number of time-points and channels
[~, ~, ~, numTimePoints, numChannels]=size(data);
% Create empty matrix for denoised data
denoisedData=zeros(size(data));
 
% Loop through time-points
for t=1:numTimePoints
  
   % Loop through channels
    for c=1:numChannels
 
        % Pass data for time-point and channel to ImageJ
        MIJ.createImage(squeeze(data(:,:,:,t,c)));
        
        % Run pureDenoise plugin
        MIJ.run('PureDenoise ...', paramStr);
        
        % Retrieve denoised data
        denoisedData(:,:,:,t,c)=MIJ.getCurrentImage;
        
        % Close open windows in ImageJ
        MIJ.run('Close All');
    end
     disp(['Time-point ' num2str(t) ' of '  num2str(numTimePoints) ' complete'])
end
 
% Exit MIJ
MIJ.exit

