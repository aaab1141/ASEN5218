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
rayunitvecs = raydirection(0,appgrid,paxis);

% seemyrays(rayunitvecs,appgrid)

%% Calculate N0, M, and P for each gridpoint
% N0 is the same for each gridpoint
N0 = -f*(1+ecc)*paxis;

% M is the same for each gridpoint
M = eye(3) - ecc^2*(paxis*paxis');

% P is different for each gridpoint so we will define P for each gridpoint
% because of the way we defined our coordinate system P is the gridpoint coordinate
% reshape the three layer matrix into things that are plottable
x = reshape(appgrid(:,:,1),1,size(appgrid,1)*size(appgrid,2));
y = reshape(appgrid(:,:,2),1,size(appgrid,1)*size(appgrid,2));
z = reshape(appgrid(:,:,3),1,size(appgrid,1)*size(appgrid,2));
P = [x;y;z];

%% Calculate L for each incident ray
rts = zeros(1,size(rayunitvecs,2));
flags = false(1,size(rayunitvecs,2));
% for each node
for ind = 1:1:size(rayunitvecs,2)
    % get P, and i for the gridpoint
    PP = P(:,ind);
    i = rayunitvecs(:,ind);
    
    % set up the equation and solve (basically solve the polynomial)
    one = i'*M*i;
    two = 2*i'*(M*PP + N0);
    three = PP'*(M*PP + 2*N0);
    temp = roots([one,two,three]);
    if numel(temp) == 1
        rts(1,ind) = temp;
        flags(1,ind) = false; %false flag indicates there was only one solution and there was no problem
    else
        if sign(temp(1)) ~= sign(temp(2)) %means it hit the parabola going forward and way far away going backwards
            flags(1,ind) = false;
            if temp(1) > 0
                rts(1,ind) = temp(1);
            else
                rts(1,ind) = temp(2);
            end
        else
            flags(1,ind) = true;
        end
    end
end

seemyroots(rts,n)

%% Determine which Ls are outside the circular aperture


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
figure
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

function [] = seemyroots(rts,n)
mygrid = reshape(rts,n,n);

figure
contour(mygrid)
colorbar
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