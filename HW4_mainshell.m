% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This is the main shell script for ASEN 5218 Homework #4. It facilitates
% running the problems in the homework. To check this homework, the only 
% thing that is required is to run this script.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear
close all

% Question 1: Build ray tracing code
%% Set up variables to make this easily modular
f = 20; %m focal length
ecc = 1; %eccentricity
p0 = [0,0,0]; %vertex of the mirror
paxis = [-1;0;0]; %principle axis of the mirror
d = 10; %m aperture diameter
al = [-2*f,0,0]; %aperture location (center of aperture)

%% define the grid of rays to simulate
n = 5; %number of gridpoints per side
appgrid = zeros(n,n,3);
appgrid(:,:,1) = al(1,1); %all the x coordinates are at the aperture (2*f)
linvec = linspace(-d/2,d/2,n);
for i = 1:1:n
    appgrid(:,i,2) = linvec; %change x values
    appgrid(i,:,3) = linvec; %change y values
end

% seemygrid(appgrid)

%% Define the direction of the incoming ray at each gridpoint
rayunitvecs = raydirection(45,appgrid,paxis);

% seemyrays(rayunitvecs,appgrid)

%% Calculate M and P for each gridpoint
% M is the same for each gridpoint
M = eye(3) - ecc^2*(paxis*paxis');

% Question 2: Countour plot of optical path difference

% Question 3: Coma - the plot of optical path difference as a function of a

% Question 4: Effect of the reference sphere

%% functions
function [] = seemygrid(gridpoint3layer)
% reshape the three layer matrix into things that are plottable
x = reshape(gridpoint3layer(:,:,1),1,size(gridpoint3layer,1)*size(gridpoint3layer,2));
y = reshape(gridpoint3layer(:,:,2),1,size(gridpoint3layer,1)*size(gridpoint3layer,2));
z = reshape(gridpoint3layer(:,:,3),1,size(gridpoint3layer,1)*size(gridpoint3layer,2));

% make a plot
plot3(x,y,z,'*');
xlabel('x');ylabel('y');zlabel('z')
end

function [] = seemyrays(rayvecs,gridpoint3layer)
% reshape the three layer matrix into things that are plottable
x = reshape(gridpoint3layer(:,:,1),1,size(gridpoint3layer,1)*size(gridpoint3layer,2));
y = reshape(gridpoint3layer(:,:,2),1,size(gridpoint3layer,1)*size(gridpoint3layer,2));
z = reshape(gridpoint3layer(:,:,3),1,size(gridpoint3layer,1)*size(gridpoint3layer,2));

figure; hold on
for i = 1:1:size(rayvecs,2)
    % get the gridpoint for the ray vector
    X = [x(i), x(i) + rayvecs(1,i)];
    Y = [y(i), y(i) + rayvecs(2,i)];
    Z = [z(i), z(i) + rayvecs(3,9)];
    
    % plot the vector
    plot3(X,Y,Z,'-b*')
end
xlabel('x');ylabel('y');zlabel('z')
axis equal
end

function [rayunitvecs] = raydirection(tiltangle,gridpoint3layer,paxis)
% set aside space for all the unit vectors
rayunitvecs = zeros(3,size(gridpoint3layer,1)*size(gridpoint3layer,2));

if tiltangle == 0
    % establish the nominal direction of the rays (antiparallel to principal axis
    for i = 1:1:size(rayunitvecs,2)
        rayunitvecs(:,i) = -paxis;
    end
else
    % Tilt the unit vectors by the specified angle
    for i = 1:1:size(rayunitvecs,2)
        rayunitvecs(2,i) = sind(tiltangle);
        rayunitvecs(1,i) = cosd(tiltangle);
    end
end
end