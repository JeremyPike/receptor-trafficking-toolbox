
function exampleFigures(bulkColocStats, bandStats, maxDepth, bandWidth)
 
% EXAMPLEFIGURES creates 4 example figures using the results of the
% workflow
%
% INPUT bulkColocStats: struct containing the calculated colocalization 
%                       measures for the whole cell 
%       bandStats: struct containing the calculated band based measures
%       maxDepth: maximum distance from cell edge for plotting
%       bandWidth: Width of each band in microns
%
%
% created by: Jeremy Pike
% DATE: 21-Oct-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
disp('Creating example figures...') 
 
% Calculate number of bands which should be plotted
numBands = floor(maxDepth/bandWidth);
% Acquisition times as vector
time = [0 10 20 30];
% Acquisition times as text
timeLabels = {'0 minutes', '10 minutes', '20 minutes', '30 minutes'};
numTimePoints = length(time);
% Maximum edge distance for each band
bandDist = (1:numBands)*bandWidth;
 
%% Figure for the bulk colocalization measures for the cell
figure('units','normalized','outerposition',[0.25 0.25 0.5 0.5])
% Title
suptitle('Colocalization measures for whole cell')
% M1 plot
subplot(1,3,1)
plot(time, bulkColocStats.M1,'.-') 
xlabel('Time (minutes)')
ylabel('M1')
% M2 plot
subplot(1,3,2)
plot(time, bulkColocStats.M2,'.-') 
xlabel('Time (minutes)')
ylabel('M2')
% Pearson plot
subplot(1,3,3)
plot(time, bulkColocStats.Pearson,'.-') 
xlabel('Time (minutes)')
ylabel('Pearson')
 
%% Figure for the percentage of receptor distribution in each banded volume
figure('units','normalized','outerposition',[0.25 0.25 0.5 0.5])
% title
suptitle('Analysis channel signal distribution')
% Subplot for measurements without volume correction
subplot(1,2,1)
hold all
% Plot all time-points on single axis
for t = 1:numTimePoints
    plot(bandDist, bandStats.signalPercentage{t,1}(1:numBands),'.-') 
end
xlabel('Edge distance (\mu m)')
ylabel('% Signal')
% Subplot for measurements with volume correction
subplot(1,2,2)
hold all
% Plot all time-points on single axis
for t = 1:numTimePoints
    plot(bandDist, bandStats.signalPercentageVolumeCor{t,1}(1:numBands),'.-')   
end
xlabel('Edge distance (\mu m)')
ylabel('% Signal (volume corrected)')
legend(timeLabels)
 
%% Figure for the percentage of the co-occurrence map in each banded volume
figure('units','normalized','outerposition',[0.25 0.25 0.5 0.5])
% title
suptitle('Co-occurrence map signal distribution')
% Subplot for measurements without volume correction
subplot(1,2,1)
hold all
% Plot all time-points on single axis
for t = 1:numTimePoints
    plot(bandDist, bandStats.cooccurrenceMapPercentage{t,1}(1:numBands),'.-') 
end
xlabel('Edge distance (\mu m)')
ylabel('% Signal')
% Subplot for measurements with volume correction
subplot(1,2,2)
hold all
% Plot all time-points on single axis
for t = 1:numTimePoints
    plot(bandDist, bandStats.cooccurrenceMapPercentageVolumeCor{t,1}(1:numBands),'.-')   
end
xlabel('Edge distance (\mu m)')
ylabel('% Signal (volume corrected)')
legend(timeLabels)
 
%% Figure showing colocalization measures for each banded volume
figure('units','normalized','outerposition',[0.25 0.25 0.5 0.5])
% Title
suptitle('Colocalization measures for each band')
% M1 subplot
subplot(1,3,1)
hold all
% Plot all time-points on single axis
for t = 1:numTimePoints
    plot(bandDist, bandStats.M1{t,1}(1:numBands),'.-') 
end
xlabel('Edge distance (\mu m)')
ylabel('M1')
% M2 subplot
subplot(1,3,2)
hold all
% Plot all time-points on single axis
for t = 1:numTimePoints
    plot(bandDist, bandStats.M2{t,1}(1:numBands),'.-')   
end
xlabel('Edge distance (\mu m)')
ylabel('% M2')
% Pearson subplot
subplot(1,3,3)
hold all
% Plot all time-points on single axis
for t = 1:numTimePoints
    plot(bandDist, bandStats.Pearson{t,1}(1:numBands),'.-')   
end
xlabel('Edge distance (\mu m)')
ylabel('% Pearson')
legend(timeLabels)
