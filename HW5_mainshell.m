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
outalph = 2*pi/numoutnodes;
inalph = 2*pi/numinnnodes;

%define the origin of our reflector at the bottom of the paraboloid
origin = [0;0;0]; %m

% determine the position of the outer nodes. These are only dependent on the number of bars and the projected diameter
outnodes = zeros(3,numoutnodes);
for i = 0:1:numoutnodes-1
    outnodes(:,i+1) = getnode(diameter/2*cos(i*outalph),diameter/2*sin(i*outalph));
end

% determine b5
b5 = nodedist(outnodes(:,1),outnodes(:,2));
disp(['b5 = ',num2str(b5),' m'])

% Find the value of b1 that makes the angular defect the same 
b1projlengths = linspace(0.1,9.9,10); % projected length of b1
b1lengths = zeros(1,length(b1projlengths));
b2lengths = b1lengths;
b3lengths = b1lengths;
b4lengths = b1lengths;
defori = b1lengths;
innnodes = zeros(3,numinnnodes);

for i = 1:1:length(b1projlengths)
    b1proj = b1projlengths(i); %assign a value to b1
    
    % calculate b2 based on b1
    for n = 0:1:numinnnodes-1
        innnodes(:,n+1) = getnode(b1proj*cos(n*inalph),b1proj*sin(n*inalph));
    end
    b1 = nodedist(origin,innnodes(:,1));
    b1lengths(i) = b1;
    b2 = nodedist(innnodes(:,1),innnodes(:,2));
    b2lengths(i) = b2; %store for reference
    
    % calculate the angular defect at the origin
    dorigin = 2*pi - 12*asin(b2/2/b1);
    defori(i) = dorigin;
    
    
end

figure
plot(b1projlengths,defori)

% Functions
function [ dist ] = nodedist(node1,node2)
% This function finds the euclidean distance between two nodes (ie a bar length)
dist = sqrt(sum((node1 - node2).^2));
end

function [ node ] = getnode(x,y)
% This function returns the node on the parabola from an x,y position in the plane
z = getz(x,y);
node = [x;y;z];
end

function [ z ] = getz(x,y)
% This function gives the z coordinate of a node in with the given x,y position
x = getparax(x,y);
z = 1/(4*20)*x^2;
end

function [ d ] = getparax(x,y)
% This function gets the 'x' coordinate of the parabola based on the nodes x,y position
d = sqrt(x^2+y^2);
end