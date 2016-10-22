
%EXAMPLEWORKFLOW implements the key workflows for the quantification of 
% sub-cellular receptor distribution, and colocalization with endosomes.
% The example confocal time-lapse movie consists of a HeLa cell expressing
% EGFR-EGFP and rab5-mRFP. The cell was stimulated with EGF immediately
% prior to imaging.
%
% REMARKS: Please refer to README.txt and the manuscript by Pike et al. for
%          more details. Ensure you have all the necessary dependencies
%          correctly configured (see REEDME.txt) and that the ImageJ and
%          MIJ paths are correctly specified (see below).
%
% created by: Jeremy Pike
% DATE: 16-Oct-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%% User specified parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
close all; clear all; clc
 
% Name of example data file
filePath = 'exampleTimelapseAquisition.nd2';
 
%%  Denoising parameters
% Path of ImageJ installation
denoiseParams.IJpath = 'C:\ImageJ';
% Location of mij.jar
denoiseParams.MIJpath = 'C:\Program Files\MATLAB\R2014a\java\mij.jar';
% Number of frames
denoiseParams.numFrames = 3;
% Number of spin-cycles
denoiseParams.numSpinCycles = 4;
 
%% Segmentation parameters
% Number of clusters for K-means based segmentation estimate
segmentationParams.numClusters = 4;
% Radius in microns for morphological operations on segementation estimate
segmentationParams.SEradius = 1.5;
 
% For details concerning the following DRLSE parameters see the paper by 
% Li et al. IEEE transactions on image processing 19.12 (2010): 3243-3254
% weight of the weighted length term
segmentationParams.lambda = 8;
% weight of the weighted area term
segmentationParams.alfa = 0;
% width of Dirac Delta function
segmentationParams.epsilon = 1.5;
% time step for update weighting
segmentationParams.timeStep = 1;
% power, p, for edge indicator function
segmentationParams.gradValue = 2;
% Size of binary step for initialization of level set function
segmentationParams.c0 = 2;    
% Maximum number of iterations for level set evolution
segmentationParams.iterMax = 50;
% Minimum number of iterations for level set evolution
segmentationParams.iterMin = 5;
% Mimnimum volume percentage change for stopping condition
segmentationParams.stopPercent = 0.01;
 
%% Background subtraction parameters
% Radius in microns of the spherical structural element used for 
% background subtraction
backgroundSubtractParams.SEradius = 1;
 
%% Band-based analysis parameters
% Band width in microns
bandParams.bandWidth = 0.5;
% Channel for calculation of signal distribution
bandParams.analysisChannel = 1;
 
%% Visualisation and plotting parameters
% Maximum distance from segmentation edge for plotting (in microns)
visParams.maxDistance = 5;
% Time-point for visualisation
visParams.timePoint = 2;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Workflow implementation
 
% Load the example file using the Bio-Formats Library
[rawData, meta] = loadRawData(filePath, 1);
% Denoise with PURE-LET scheme (within IJ using MIJ)
 denoisedData = denoise(rawData, denoiseParams);
% Segment cell using 3D DRLSE 
cellSegmentation = cellularSegment(denoisedData(:,:,:,:,1), meta, segmentationParams);
% Calculate background subtracted data using a 3D rolling ball approach
dataBackgroundSubtracted = backgroundSubtract3D(denoisedData, backgroundSubtractParams.SEradius, meta.voxelSizeX, meta.voxelSizeZ);
% Find threshold values for single isolation (Otsu approach)
thresholdsOtsu = otsuThresholds(dataBackgroundSubtracted, cellSegmentation);
% Calculate bulk colocalization statistics for the entire cell
bulkColocStats = calcColocStats(dataBackgroundSubtracted, cellSegmentation, thresholdsOtsu, true);
% Calculate band based statistics
bandStats = calcBandStats(dataBackgroundSubtracted, cellSegmentation, thresholdsOtsu, meta, bandParams);
 
%% Plot results and visualisation
 
% Construct some example figures using the results of the workflow
exampleFigures(bulkColocStats, bandStats, visParams.maxDistance, bandParams.bandWidth)
% Create some example visualisations
exampleVisualisations(dataBackgroundSubtracted, cellSegmentation, thresholdsOtsu, meta, visParams.timePoint, bandParams.analysisChannel)

