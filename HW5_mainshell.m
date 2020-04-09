% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% This is the main shell script for ASEN 5218 Homework #4. It facilitates
% running the problems in the homework. To check this homework, the only
% thing that is required is to run this script.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear
close all

%% Question 1
% some constants
diameter = 20; %m
numoutnodes = 12;
numinnnodes = 6;

%define the origin of our reflector at the bottom of the paraboloid
origin = [0;0;0]; %m

% determine the position of the outer nodes. These are only dependent on the number of bars and the projected diameter
outnodes = zeros(3,numoutnodes);
outalph = 2*pi/numoutnodes;
for i = 0:1:numoutnodes-1
    outnodes(:,i+1) = getnode(diameter/2*cos(i*outalph),diameter/2*sin(i*outalph));
end





% Functions
function [ z ] = getz(x,y)
% This function gives the z coordinate of a node in with the given x,y position
x = getparax(x,y);
z = 1/(4*20)*x^2;
end

function [ d ] = getparax(x,y)
% This function gets the 'x' coordinate of the parabola based on the nodes x,y position
d = sqrt(x^2+y^2);
end

function [ node ] = getnode(x,y)
% This function returns the node on the parabola from an x,y position in the plane
z = getz(x,y);
node = [x;y;z];
end