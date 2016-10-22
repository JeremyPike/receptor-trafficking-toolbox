
function phi = drlse_edge_3D(phi, g, vx, vy, vz, lambda, mu, alfa, epsilon, timestep, deltaXY, deltaZ)
 
% DRLSE_EDGE_3D is a 3D implementation of edge based Distance Regularized 
% Level Set Evolution (DRLSE) as described in 2D by Li et al. IEEE 
% transactions on image processing 19.12 (2010): 3243-3254.
%
% INPUT phi: 3D matrix containing level set function to be updated by DRLSE           
%       g: edge indicator function
%       vx, vy, vz: lateral and axial components of the gradient of g  
%       lambda: weight of the weighted length term
%       mu: weight of distance regularization term
%       alfa: weight of the weighted area term
%       epsilon: width of Dirac Delta function
%       timestep: time step for update weighting
%       delataXY, deltaZ: Scaling for lateral and axial dimensions
%
% OUTPUT phi: 3D matrix containing level set function after single
%             iteration of DRLSE
%
% REMARKS: This code is a modified version of the 2D implementation available
%          at http://www.imagecomputing.org/~cmli/
%
% created by: Original 2D implementation by Chunming Li, 3D modification by
%             Jeremy Pike
% DATE: 19-Oct-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Enforce Neumann boundary conditions on level set function, phi
phi=NeumannBoundCond(phi);
 
% Calculate first order dervivatives of phi
[phi_x,phi_y,phi_z]=gradient(phi,deltaXY,deltaXY,deltaZ);
 
% Gradient of phi
s=sqrt(phi_x.^2 + phi_y.^2+ phi_z.^2);
smallNumber=1e-10;  
 
%Normal vector
Nx=phi_x./(s+smallNumber); 
Ny=phi_y./(s+smallNumber);
Nz=phi_z./(s+smallNumber);
 
%Calculate curvature
curvature=div(Nx,Ny,Nz,deltaXY,deltaZ);
 
% Compute the distance regularization term in equation with the double-well 
% potential.
distRegTerm=distReg_p2(phi,phi_x,phi_y,phi_z,s,deltaXY,deltaZ);  
 
% Dirac delta function
diracPhi=Dirac(phi,epsilon);
 
% Compute area term
areaTerm=diracPhi.*g; % balloon/pressure force
 
%Compute edge term
edgeTerm=diracPhi.*(vx.*Nx+vy.*Ny+vz.*Nz) + diracPhi.*g.*curvature;
 
% Update phi
phi=phi + timestep*(mu*distRegTerm + lambda*edgeTerm + alfa*areaTerm);
 
 
 
function f = distReg_p2(phi,phi_x,phi_y,phi_z,s,deltaXY,deltaZ)
% compute the distance regularization term with the double-well potential p2 in equation (16)
 
a=(s>=0) & (s<=1);
b=(s>1);
ps=a.*sin(2*pi*s)/(2*pi)+b.*(s-1);  % compute first order derivative of the double-well potential p2 in equation (16)
dps=((ps~=0).*ps+(ps==0))./((s~=0).*s+(s==0));  % compute d_p(s)=p'(s)/s in equation (10). As s-->0, we have d_p(s)-->1 according to equation (18)
f = div(dps.*phi_x - phi_x, dps.*phi_y - phi_y, dps.*phi_z - phi_z,deltaXY,deltaZ) + 6*del2(phi,deltaXY,deltaXY,deltaZ);  
 
function f = div(nx,ny,nz,deltaXY,deltaZ)
%Divergence
[nxx,~,~]=gradient(nx,deltaXY,deltaXY,deltaZ);  
[~,nyy,~]=gradient(ny,deltaXY,deltaXY,deltaZ);
[~,~,nzz]=gradient(nz,deltaXY,deltaXY,deltaZ);
f=nxx+nyy+nzz;
 
function f = Dirac(x, sigma)
%Dirac-delta function
f=(1/2/sigma)*(1+cos(pi*x/sigma));
b = (x<=sigma) & (x>=-sigma);
f = f.*b;
 
function g = NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[rows, cols , numSlices] = size(f);
g = f;
g([1 rows],[1 cols],[1 numSlices]) = g([3 rows-2],[3 cols-2],[3 numSlices-2]);  
g([1 rows],2:end-1) = g([3 rows-2],2:end-1);          
g(2:end-1,[1 cols],[1 numSlices]) = g(2:end-1,[3 cols-2],[3 numSlices-2]);  
g([1 rows],2:end-1,[1 numSlices]) = g([3 rows-2],2:end-1,[3 numSlices-2]);      
g([1 rows],[1 cols],2:numSlices-1) = g([3 rows-2], [3 cols-2], 2:numSlices-1);
