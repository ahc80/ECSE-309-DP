%% ECSE 309 Example Torque Calculation due to Lorentz Force
% Scripts are provided -as is- and it is up to the student to use and
% modify them as needed to complete assigned work

clc
close all
clear all 

currentPath = pwd;                          % store current path

%% Generate Computation Space
% Build a workspace of point of interest in planar arrangement below the
% plane of the magnets with 32 sensors along x and 32 sensors along y.
% The plane of sensors is repeated 32 times along z
% Sensors lies in a equally spaced grid with step 0.001 m (1 mm).
xres = 32;
yres = 32;
zres = 32;
spacing = 0.001;    %spacing in meters between computed points
[SensorPosMatrix,WorkspaceDim] = buildWorkspace([xres yres zres],spacing,2);

%% Generate Magnets
% Positions of the magnets: the first three coordinates (per line) are the
% positions, while the last three coordinates are the dipole magnetic
% moment unit vectors, e.g. the orientation of each magnet (in this case 
% they are pointing downward, perpedicular to the sensor plane)
MagPos = [0   0    0   0 0 -1;
          0   0  0.003 0 0 -1;
          0   0  0.006 0 0 -1]'; 
      
%MagPos = [0.01  0.00 0.000 0 1 0]';    %example for a single magnet   

MM = size(MagPos,2); 
D = 0.012*ones(MM,1); % Diameter of the magnets
L = 0.003*ones(MM,1); % Length of the magnets
M = 1.25/(4*pi*1e-7)*ones(MM,1); % Magnetization of the magnets [A/m]

%% Generate Current Carrying Wire
% a wire in 3D space may be defined using a function over the computation
% space or a wire can be specified with piecewise coordinates
% ---- snap the wire to the nearest simulation points ----
% ---- break the wire where it intersects magnets ---- 

% if you generate a wire manually make sure to put points the same spacing
% apart as in the sensor grid above
% Otherwise you can specify the wire vertices and use the function to
% generate a piecewise linear wire

I = 50;    %current in wire (A)

%%
%Testing out helix wire

% Parameters for the helix
R = 0.008; % Radius of the spiral, just outside the magnets' diameter
c = 0.0015; % Pitch of the helix, can be the same as magnet thickness for tight spiral
t = linspace(0, 4*pi*3, 100); % Parameter t, adjust 2*pi*3 for the number of turns

% Helix equations
x = R * cos(t);
y = R * sin(t);
z = c * t / (2*pi); % Normalize by 2*pi to make 'c' represent the pitch directly

% Generating the spiral wire points
wirePoints = [x; y; z]';
%%

% now at each point on the wire, generate a vector of current components
currVec = 0.*wirePoints;
for k = 1:length(wirePoints)-1
    uvec = wirePoints(k+1,:)-wirePoints(k,:);
    mag = sqrt(uvec(1)^2 + uvec(2)^2 + uvec(3)^2);
    currVec(k,:) = uvec ./ mag;
end

currVec(end,:) = currVec(end-1,:);
currVec = currVec .* I;   %scale to the magnitude of current at each point

%% Run Simulation
disp('Computation of the magnetic field with GenerateReadings.m :')
tic
Readings = GenerateReadings(MagPos',SensorPosMatrix,M,D,L);
toc

figure(1)
scatter3(SensorPosMatrix(:,1),SensorPosMatrix(:,2),SensorPosMatrix(:,3),0.5,[0.9 0.9 0.9],'.')
hold on
plot3(wirePoints(:,1), wirePoints(:,2), wirePoints(:,3),'--k')

for i = 1:MM
    hold on
    drawCylindricalMagnet(L(i),D(i),MagPos(:,i),'texture','axial')
    % hold on
    % drawVec(MagPos(1:3,i),MagPos(4:6,i)*0.005)
end
hold on
showFrame([0 0 0],[1 0 0]/500,[0 1 0]/500,[0 0 1]/500,...
                'sphereradius',0.0002,...
                'arrowheadlength',0.05,...
                'arrowheadwidth',0.05)

xcomp = Readings(:,1);
ycomp = Readings(:,2);
zcomp = Readings(:,3);

% q = quiver3(SensorPosMatrix(:,1),SensorPosMatrix(:,2),SensorPosMatrix(:,3),...
%     xcomp,ycomp,zcomp,...
%     'r','LineWidth',1);
% q.ShowArrowHead = 'on';
% %q.Marker = '.';
% axis equal
% hold off

colormap turbo
quiverC3D(SensorPosMatrix(:,1),SensorPosMatrix(:,2),SensorPosMatrix(:,3),...
    xcomp,ycomp,zcomp,...
    1,5000,'scaleMode','auto')
axis equal
%clim("auto");
hold off
% axis([-80 80 -80 80 -80 80]/1000)
setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera')
grid off
axis off

%% Estimate Wire Torque
% because the wire is only in one plane only 2D integration is needed 
% extract the magnetic field at points along the wire
% this is much easier if the wire is on the same grid that the field is
% calculated on

[~,I] = pdist2(SensorPosMatrix,wirePoints,'euclidean','Smallest',1); % the indices in I are the rows in SensorPosMatrix to evaluate B on
fieldVec = Readings(I,:);
forceVec = cross(currVec,fieldVec);

fToty = cumtrapz(forceVec(:,2)).*(spacing);    %total force in N along y axis
%fToty = cumtrapz([0:spacing:0.002*(length(forceVec)-1)],forceVec(:,2));

% to determine torque we need a torque moment arm at each location of the
% force vector
% assume the pivot point is at x,y = 0,0 and the pivot is around z
tMomVec = wirePoints;
torqueVec = cross(tMomVec,forceVec);
torqueZ = torqueVec(:,3);   %Z component of torque vector is the torque along the rotation axis for each piece of wire analyzed
tTotz = cumtrapz(torqueVec(:,3)).*(spacing);  %cumulative torque along wire
zTorque = abs(tTotz(end));                    %total torque moment along Z axis

f2 = figure(2);
scatter3(SensorPosMatrix(:,1),SensorPosMatrix(:,2),SensorPosMatrix(:,3),0.5,[0.9 0.9 0.9],'.')
hold on
plot3(wirePoints(:,1), wirePoints(:,2), wirePoints(:,3),'-k','LineWidth',1)

vecType = 'field';  %options: force, torque, field, moment

scaling = 'log';
xref = wirePoints(:,1);
yref = wirePoints(:,2);
zref = wirePoints(:,3);

if strcmp(vecType,'force') 
    xcomp = forceVec(:,1);
    ycomp = forceVec(:,2);
    zcomp = forceVec(:,3);
elseif strcmp(vecType,'torque') 
    xcomp = torqueVec(:,1);
    ycomp = torqueVec(:,2);
    zcomp = torqueVec(:,3);

    xref = 0 .* xref;
    yref = 0 .* yref;
    zref = 0 .* zref;
    %scaling = 'none';   %allows vectors to be seen

elseif strcmp(vecType,'field') 
    xcomp = fieldVec(:,1);
    ycomp = fieldVec(:,2);
    zcomp = fieldVec(:,3);
elseif strcmp(vecType,'moment') 
    xcomp = tMomVec(:,1);
    ycomp = tMomVec(:,2);
    zcomp = tMomVec(:,3);
    
    xref = 0 .* xref;
    yref = 0 .* yref;
    zref = 0 .* zref;
    scaling = 'none';   %allows vectors to be seen

else
    error('select a vector field to plot')
end

for i = 1:MM    
    %drawCylindricalMagnet(L(i),D(i),MagPos(:,i),'texture','normal')
    % hold on
    drawVec(MagPos(1:3,i),MagPos(4:6,i)*0.005)  %vectors of magnetic moments
end
% frame axes: red = x, green = y, blue = z
showFrame([0 0 0],[1 0 0]/500,[0 1 0]/500,[0 0 1]/500,...
                'sphereradius',0.0002,...
                'arrowheadlength',0.05,...
                'arrowheadwidth',0.05)

colormap turbo
quiverC3D(xref,yref,zref,xcomp,ycomp,zcomp,...
    0.00001,5000,'scaleMode',scaling,'LineWidth', 3)
%quiver3(xref,yref,zref,xcomp,ycomp,zcomp)
axis equal
%clim("auto");
hold off
% axis([-80 80 -80 80 -80 80]/1000)
setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera')
grid off
axis off