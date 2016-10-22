function [T1, T2] = costesThresh(dataC1, dataC2, ROI, updateStep)
 
%costesThresh implements Costes thresholding as described by
% Costes, Sylvain V., et al. Biophysical journal 86.6 (2004): 3993-4003.
%
% INPUT dataC1: Matrix containing data for channel 1
%       dataC2: Matrix containing data for channel 2
%       ROI: Boolean matrix containing analysis ROI. 
%       updateStep: Update for each for iteration down regression line
%
% OUTPUT T1: threshold for channel 1
%        T2: threshold for channel 2
%
% REMARKS: Requires linortfit2 (http://uk.mathworks.com/matlabcentral/
% fileexchange/16800-orthogonal-linear-regression)
%
% created by: Jeremy Pike
% DATE: 15-Oct-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Convert data and ROI to vectors
dataC1=dataC1(:);
dataC2=dataC2(:);
ROI=ROI(:);
 
% Delete voxels not within the ROI
dataC1(ROI==0)=[];
dataC2(ROI==0)=[];
 
% Find maximum values (within ROI) for both channels
maxC1=max(dataC1);
maxC2=max(dataC2);
   
% Use orthogonal linear regression to find line of best fit for joint
% histogram 
fit=linortfit2(dataC1, dataC2);
 
% if |fit(1)|<=1 then move threshold down C1 and calculate C2 as 
% dataC2 = fit(1)*dataC1 + fit(2)
if abs(fit(1))<=1
    [T1, T2] = costesOpt(dataC1, dataC2, fit(1), fit(2), maxC1, updateStep);      
else
% if |fit(1)|>1 then move threshold down C2 and calculate C1 as 
% dataC1 = dataC2/fit(1) - fit(2)/fit(1)        
    [T2, T1] = costesOpt(dataC2, dataC1, 1/fit(1), -fit(2)/fit(1), maxC2, updateStep);
    
end
 
   
end
 
function [T1, T2] = costesOpt(dataC1, dataC2, a, b, start, updateStep)
 
    % Counting variable for number of iterations
    count = 0;
    % Set starting value for T1
    T1 = start;
    % Find minum values for both channels
    minC1=min(dataC1);
    minC2=min(dataC2);
    
    while true
      
        %Increase iteration
        count=count+1;
        %Update threshold
        T1=T1-updateStep;
        %Calculate T2 using line of best fit
        T2=a*T1+b;
        %Store current thresholds in vector containing all values
        T1All(count,1)=T1;
        T2All(count,1)=T2;
        %Calculate Pearson coefficient using all voxels below either
        %threshold
        pearsonCoefficient(count,1)=pearsonCostes(dataC1, dataC2, T1 ,T2);
        
        %If threshold values are less than the minimum of either channels
        if (T1<=minC1) || (T2<minC2)
            %Find the minimum Pearson coefficient (with highest T1)
            pInd = find(pearsonCoefficient==min(pearsonCoefficient),1, 'last' );
            %Use the corresponding threshold values
            T1=T1All(pInd);
            T2=T2All(pInd);           
           %Stop iterations
            break
        end
 
        % If the Pearson coefficient <=0     
        if pearsonCoefficient(count,1) <=0
            %Use the previous threshold values
            T1=T1+updateStep;
            T2=a*T1+b;
            %Stop iterations
            break
        end     
        
    end
    
end
 
function R = pearsonCostes(dataC1, dataC2, T1 ,T2)
 
%Create mask to identify all pixels bellow either threshold
mask = zeros(size(dataC1));
mask(dataC1 < T1) = 1;
mask(dataC2 < T2) = 1;
 
%Delete all other pixels
dataC1(mask==0)=[];
dataC2(mask==0)=[];
 
% Calculate Pearson coefficient 
R =(sum((dataC1-mean(dataC1)).*(dataC2-mean(dataC2))))/(sqrt(sum((dataC1-mean(dataC1)).^2)*sum((dataC2-mean(dataC2)).^2)));
 
end
 
 
 

