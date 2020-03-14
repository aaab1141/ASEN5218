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
n = 50; %number of gridpoints per side
appgrid = zeros(n,n,3);
appgrid(:,:,1) = al(1,1); %all the x coordinates are at the aperture (2*f)
linvec = linspace(-d/2,d/2,n);
for i = 1:1:n
    appgrid(:,i,2) = linvec; %change x values
    appgrid(i,:,3) = linvec; %change y values
end

% seemygrid(appgrid)

%% Define the direction of the incoming ray at each gridpoint
rayunitvecs = raydirection(10e-2,appgrid,paxis);


% Question 2: Countour plot of optical path difference

% Question 3: Coma - the plot of optical path difference as a function of a

% Question 4: Effect of the reference sphere

function [] = seemygrid(threelayer)
% reshape the three layer matrix into things that are plottable
x = reshape(threelayer(:,:,1),1,size(threelayer,1)*size(threelayer,2));
y = reshape(threelayer(:,:,2),1,size(threelayer,1)*size(threelayer,2));
z = reshape(threelayer(:,:,3),1,size(threelayer,1)*size(threelayer,2));

% make a plot
plot3(x,y,z,'*');
xlabel('x');ylabel('y');zlabel('z')
end

function [rayunitvecs] = raydirection(tiltangle,gridpoint3layer,psi)
% set aside space for all the unit vectors
rayunitvecs = zeros(3,size(gridpoint3layer,1)*size(gridpoint3layer,2));

if tiltangle == 0
    % establish the nominal direction of the rays (antiparallel to principal axis
    for i = 1:1:size(rayunitvecs,2)
        rayunitvecs(:,i) = -psi;
    end
else
    % Tilt the unit vectors by the specified angle
    for i = 1:1:size(rayunitvecs,2)
        rayunitvecs(2,i) = sind(tiltangle);
        rayunitvecs(1,i) = cosd(tiltangle);
    end
end
end