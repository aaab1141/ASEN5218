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

%% Calculate path lengths for each incident ray
point2mirror = zeros(1,size(rayunitvecs,2));
reflectpoint = zeros(size(rayunitvecs,1),size(rayunitvecs,2));
reflectvec = zeros(size(rayunitvecs,1),size(rayunitvecs,2));
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
        point2mirror(1,ind) = temp;
        flags(1,ind) = false; %false flag indicates there was only one solution and there was no problem
    else
        if sign(temp(1)) ~= sign(temp(2)) %means it hit the parabola going forward and way far away going backwards
            flags(1,ind) = false;
            if temp(1) > 0
                point2mirror(1,ind) = temp(1);
            else
                point2mirror(1,ind) = temp(2);
            end
        else
            flags(1,ind) = true;
        end
    end
    
    % find the point where the ray hit the mirror
    reflectpoint(:,ind) = [x(ind);y(ind);z(ind)] + point2mirror(ind)*rayunitvecs(:,ind);
    
    % calculate the normal at the reflection point
    N = N0 + M*reflectpoint(:,ind);
    Nhat = -sign(rayunitvecs(:,ind)'*N)*N/norm(N);
    
    % find the reflected vector direction
    R = eye(3) - 2*(Nhat*Nhat');
    reflectvec(:,ind) = R*rayunitvecs(:,ind);
    
    % find the path length of the reflected ray to the reference sphere
    % use the same thing as ray hitting the parabola except this time f = radius and e = 0
    % or maybe there is a simpler geometric way since the sphere is easily defined and we dont care about reflectance
    
    
    % find the path length of the reflected ray to the plane of the focus that is tangent to the vertex of the mirror
    % I have a feeling there is an easier geometric way here knowing where the plane is defined and that we dont care about reflectance
    
    
end

% seemyroots(point2mirror,n)
seemymirror(reflectpoint)
seemyreflectvecs(reflectvec,reflectpoint,appgrid)

%% Add the contributions to the path length of each ray


%% Determine which Ls are outside the circular aperture and set L = 0


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

function [] = seemymirror(reflectpoint)
% make the plot
figure
plot3(reflectpoint(1,:),reflectpoint(2,:),reflectpoint(3,:),'*b')
xlabel('x');ylabel('y');zlabel('z')
axis equal
end

function [] = seemyreflectvecs(reflectvec,reflectpoints,gridpoint3layer)
% reshape the three layer matrix into things that are plottable
x = reshape(gridpoint3layer(:,:,1),1,size(gridpoint3layer,1)*size(gridpoint3layer,2));
y = reshape(gridpoint3layer(:,:,2),1,size(gridpoint3layer,1)*size(gridpoint3layer,2));
z = reshape(gridpoint3layer(:,:,3),1,size(gridpoint3layer,1)*size(gridpoint3layer,2));

m = 25; %arbitrary multiplier for visual ease

figure; hold on
for i = 1:1:size(reflectvec,2)
    X = [x(i), reflectpoints(1,i), reflectpoints(1,i) + m*reflectvec(1,i)];
    Y = [y(i), reflectpoints(2,i), reflectpoints(2,i) + m*reflectvec(2,i)];
    Z = [z(i), reflectpoints(3,i), reflectpoints(3,i) + m*reflectvec(3,i)];
    
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