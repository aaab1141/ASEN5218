% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% This is the main shell script for ASEN 5218 Homework #4. It facilitates
% running the problems in the homework. To check this homework, the only
% thing that is required is to run this script.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear
close all
format
dashline = '--------------------------------------------------------------------------------';

% Build ray tracing code
%% Set up variables to make this easily modular
f = 20; %m focal length
ecc = 1; %eccentricity
% p0 = [0,0,0]; %vertex of the mirror
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

% seemyrays(rayunitvecs,appgrid) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
end

% seemyroots(point2mirror,n) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemymirror(reflectpoint) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemyreflectvecs(reflectvec,reflectpoint,appgrid) %%%%%%%%%%%%%%%%%%%%%%%

%% Find the path lengths to the reference sphere
% for this to work we have to pretend that the sphere vertex is at [0,0,0]
% Define the reference sphere
sphere.r = 1; %m sphere radius
sphere.x = -20;
sphere.y = 0;
sphere.z = 0;
sphere.ecc = 0; %define sphere eccentricity
sphere.f = sphere.r; %define the focus of the sphere which we know is at the center 1 radius away
% sphere.vertex = sphere.center + [sphere.f;0;0];
sphere.paxis = paxis; %should be the same as the mirror axis
sphere.p = sphere.f*(1 + sphere.ecc); 

% define N0 and M for the reference sphere
sphere.N0 = -sphere.p*sphere.paxis;
sphere.M = eye(3); %since e is zero then the principle axis doesnt matter

% move the reflect points such that they match our new CS where the vertex of the sphere is [0,0,0]
newreflectpoints = movemypoints_sphere(reflectpoint,sphere,0,0,0);

% find the path length of the reflected ray to the reference sphere
mirror2sphere = zeros(1,size(rayunitvecs,2));
spherepoints = zeros(size(rayunitvecs,1),size(rayunitvecs,2));
for ind = 1:1:size(rayunitvecs,2)
    % get P and i for the reflected vector
    PP = newreflectpoints(:,ind); %the origin point is the reflected point on the mirror
    i = reflectvec(:,ind); %the unit vector direction of the reflected ray is the unit vector of the incident ray
    
    % set up the equation to solve the polynomial
    one = i'*sphere.M*i;
    two = 2*i'*(sphere.M*PP + sphere.N0);
    three = PP'*(sphere.M*PP + 2*sphere.N0);
    temp = roots([one,two,three]);
    mirror2sphere(1,ind) = min(temp);
    
    % find the point where the ray hit the reference sphere
    spherepoints(:,ind) = newreflectpoints(:,ind) + mirror2sphere(ind)*reflectvec(:,ind);
   
    
end

% seemysphere(spherepoints) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Undo the coordinate transformation to figure out when it hits the reference sphere
newspherepoints = undomovemypoints_sphere(spherepoints,sphere,0,0,0);

% seemysphere(newspherepoints) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemyraypathstosphere(newspherepoints,reflectpoint,appgrid) %%%%%%%%%%%%%

%% Ray Tracing Code Questions
% calculate the ray path lengths to the focus. For this we assume the reference sphere
% can have a radius of 0 and its center is located at the focus.
[seg1,seg2] = raycodecheckup_focus();

% sum the two segments together
disp(dashline)
Ltotal = seg1+seg2;
disp('Test that all path lengths are the same if there is no angle offset and we calculate path lengths to the focus')
disp(Ltotal)
disp('We expect them all to be 60m and they all are so the ray tracing code works.')
disp(dashline)

% calculate the path lengths to the reference sphere of 1m radius
[seg1,seg2] = raycodecheckup_refsphere();
Ltotal = seg1+seg2;
disp('Test that all path lengths are the same at the reference sphere')
disp(Ltotal)
disp('We expect them to all be 59 at the reference sphere and they are so the tracing code works')
disp(dashline)

%% Question 2: Countour plot of optical path difference
% no tilt on optical path. contour plot using the difference from the path average.
[seg1,seg2,P] = changealph(0);

Ltotal = seg1 + seg2;

% remove the path lengths for the paths that are outside the cirular aperture
[X,Y,Z] = aperturetrim(Ltotal,P);

figure
contourf(X,Y,Z,[-2e-12,2e-12])
xlabel('Y, m');ylabel('Z, m')
c = colorbar;
c.Label.String = 'Path Length Deviation, m';
title({'Path Length of Aligned Rays at the Aligned Reference Sphere';'Tolerance: 2nm'})
save_fig_png('4.2a')

% tilt set to 0.01 degrees. same style contour plot
[seg1,seg2,P] = changealph(0.01);

Ltotal = seg1 + seg2;

% remove the path lengths for the paths that are outside the cirular aperture
[X,Y,Z] = aperturetrim(Ltotal,P);

figure; hold on
contourf(X,Y,Z)
plot(0,0,'k*')
xlabel('Y, m');ylabel('Z, m')
c = colorbar;
c.Label.String = 'Path Length Deviation, m';
title({'Average Path Length of 0.01^o Tilted Rays','at the Aligned Reference Sphere'})
save_fig_png('4.2b')

%% Question 3: Coma - the plot of optical path difference as a function of a
ns = -1:-1:-4;
ns = 10.^ns;
pathRMS = zeros(1,size(ns,2));

for ind = 1:1:size(ns,2)
    [seg1,seg2,P] = changealph(ns(1,ind));
    
    Ltotal = seg1 + seg2;
    
    % remove the path lengths for the paths that are outside the cirular aperture
    [~,~,~,pathRMS(1,ind)] = aperturetrim(Ltotal,P);
end

figure;
plot(ns,pathRMS,'-o')
xlabel('Tilt Angle, ^o');ylabel('Optical Path Difference RMS, m')
set(gca,'Xscale','log')
set(gca,'Yscale','log')
grid on
title('Optical Path Difference RMS over Tilt Angle')
save_fig_png('4.3a')

% plot the spot diagram when alpha is 0.01
spotdiagram(0.01);
save_fig_png('4.3b')

%% Question 4: Effect of moving the reference sphere
% we will cosider error to be rms deviation from the average optical path length
% moving the reference sphere along the X axis
delta = -1:-1:-9;
delta = 10.^delta;
pathRMS = zeros(1,size(delta,2));

for ind = 1:1:size(delta,2)
    [seg1,seg2,P] = changerefsphere(delta(1,ind),0,0);
    
    Ltotal = seg1 + seg2;
    
    % remove the path lengths for the paths that are outside the cirular aperture
    [~,~,~,pathRMS(1,ind)] = aperturetrim(Ltotal,P);
end

figure
plot(delta,pathRMS,'-o')
xlabel('\delta, m');ylabel('Optical Path Difference RMS, m')
title('Optical Path Difference RMS over X Axis Reference Sphere Delta')
set(gca,'Xscale','log')
set(gca,'Yscale','log')
axis square
grid on
save_fig_png('4.4a')

delta = -fliplr(delta);
for ind = 1:1:size(delta,2)
    [seg1,seg2,P] = changerefsphere(delta(1,ind),0,0);
    
    Ltotal = seg1 + seg2;
    
    % remove the path lengths for the paths that are outside the cirular aperture
    [~,~,~,pathRMS(1,ind)] = aperturetrim(Ltotal,P);
end

figure
plot(delta,pathRMS,'-o')
xlabel('\delta, m');ylabel('Optical Path Difference RMS, m')
title('Optical Path Difference RMS over X Axis Reference Sphere Delta')
set(gca,'Xscale','log')
set(gca,'Yscale','log')
axis square
grid on
save_fig_png('4.4b')

% moving the reference sphere along the Y axis
delta = -1:-1:-9;
delta = 10.^delta;
pathRMS = zeros(1,size(delta,2));

for ind = 1:1:size(delta,2)
    [seg1,seg2,P] = changerefsphere(0,delta(1,ind),0);
    
    Ltotal = seg1 + seg2;
    
    % remove the path lengths for the paths that are outside the cirular aperture
    [~,~,~,pathRMS(1,ind)] = aperturetrim(Ltotal,P);
end

figure
plot(delta,pathRMS,'-o')
xlabel('\delta, m');ylabel('Optical Path Difference RMS, m')
title('Optical Path Difference RMS over Y Axis Reference Sphere Delta')
set(gca,'Xscale','log')
set(gca,'Yscale','log')
axis square
grid on
save_fig_png('4.4c')

delta = -fliplr(delta);
for ind = 1:1:size(delta,2)
    [seg1,seg2,P] = changerefsphere(0,delta(1,ind),0);
    
    Ltotal = seg1 + seg2;
    
    % remove the path lengths for the paths that are outside the cirular aperture
    [~,~,~,pathRMS(1,ind)] = aperturetrim(Ltotal,P);
end

figure
plot(delta,pathRMS,'-o')
xlabel('\delta, m');ylabel('Optical Path Difference RMS, m')
title('Optical Path Difference RMS over Y Axis Reference Sphere Delta')
set(gca,'Xscale','log')
set(gca,'Yscale','log')
axis square
grid on
save_fig_png('4.4d')

% moving the reference sphere along the Z axis
delta = -1:-1:-9;
delta = 10.^delta;
pathRMS = zeros(1,size(delta,2));

for ind = 1:1:size(delta,2)
    [seg1,seg2,P] = changerefsphere(0,0,delta(1,ind));
    
    Ltotal = seg1 + seg2;
    
    % remove the path lengths for the paths that are outside the cirular aperture
    [~,~,~,pathRMS(1,ind)] = aperturetrim(Ltotal,P);
end

figure
plot(delta,pathRMS,'-o')
xlabel('\delta, m');ylabel('Optical Path Difference RMS, m')
title('Optical Path Difference RMS over Z Axis Reference Sphere Delta')
set(gca,'Xscale','log')
set(gca,'Yscale','log')
axis square
grid on
save_fig_png('4.4e')

delta = -fliplr(delta);
for ind = 1:1:size(delta,2)
    [seg1,seg2,P] = changerefsphere(0,0,delta(1,ind));
    
    Ltotal = seg1 + seg2;
    
    % remove the path lengths for the paths that are outside the cirular aperture
    [~,~,~,pathRMS(1,ind)] = aperturetrim(Ltotal,P);
end

figure
plot(delta,pathRMS,'-o')
xlabel('\delta, m');ylabel('Optical Path Difference RMS, m')
title('Optical Path Difference RMS over Z Axis Reference Sphere Delta')
set(gca,'Xscale','log')
set(gca,'Yscale','log')
axis square
grid on
save_fig_png('4.4f')

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
    Z = [z(i), z(i) + rayvecs(3,i)];
    
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

function [] = seemysphere(spherepoints)
% make the plot
figure
plot3(spherepoints(1,:),spherepoints(2,:),spherepoints(3,:),'*b')
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

function [] = seemyraypathstosphere(spherepoints,reflectpoints,gridpoint3layer)
% reshape the three layer matrix into things that are plottable
x = reshape(gridpoint3layer(:,:,1),1,size(gridpoint3layer,1)*size(gridpoint3layer,2));
y = reshape(gridpoint3layer(:,:,2),1,size(gridpoint3layer,1)*size(gridpoint3layer,2));
z = reshape(gridpoint3layer(:,:,3),1,size(gridpoint3layer,1)*size(gridpoint3layer,2));

figure; hold on
for i = 1:1:size(spherepoints,2)
    X = [x(i), reflectpoints(1,i), spherepoints(1,i)];
    Y = [y(i), reflectpoints(2,i), spherepoints(2,i)];
    Z = [z(i), reflectpoints(3,i), spherepoints(3,i)];
    
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

function [newpoints] = movemypoints_sphere(reflectpoint,sphere,deltax,deltay,deltaz)
newpoints = reflectpoint - [sphere.x + sphere.r + deltax;sphere.y + deltay;sphere.z + deltaz];
end

function [newpoints] = undomovemypoints_sphere(spherepoints,sphere,deltax,deltay,deltaz)
newpoints = spherepoints + [sphere.x + sphere.r + deltax;sphere.y + deltay;sphere.z + deltaz];
end

function [X,Y,Z,pathRMS] = aperturetrim(pathlengths,P)
appradius = 5; %m
flags = true(1,size(pathlengths,2));
n = sqrt(numel(pathlengths));
Z = zeros(1,size(pathlengths,2));

% if each ray with within the radius or not. track it if it is and exclude it if is not
for i = 1:1:size(pathlengths,2)
    if norm([P(1,i)+40;P(2,i);P(3,i)],2) <= appradius
        % track this point
        flags(1,i) = true;
    else
        % discard this point and set delta length to zero
        flags(1,i) = false;
    end
end
    
% calculate the average path length of all the points inside the aperture
pathavg = sum(pathlengths(1,flags==true))/sum(flags);

% Find the difference in the average path lengths for the relevant vectors
for i = 1:1:size(pathlengths,2)
    if flags(1,i)
        %find the difference
        Z(1,i) = pathavg - pathlengths(1,i);
    else
        %set to zero since its outside
        Z(1,i) = 0;
    end
end

% calculate path RMS
pathRMS = rms(Z(1,flags==true));

%reshape vector back into the matrix for the contour plot
Z = reshape(Z,n,n);

%meshgrid here
[X,Y] = meshgrid(linspace(-10,10,n),linspace(-10,10,n));
end

function [point2mirror,mirror2sphere] = raycodecheckup_focus()
% Set up variables to make this easily modular
f = 20; %m focal length
ecc = 1; %eccentricity
% p0 = [0,0,0]; %vertex of the mirror
paxis = [-1;0;0]; %principle axis of the mirror
d = 10; %m aperture diameter
al = [-2*f,0,0]; %aperture location (center of aperture)

% define the grid of rays to simulate
n = 5; %number of gridpoints per side
appgrid = zeros(n,n,3);
appgrid(:,:,1) = al(1,1); %all the x coordinates are at the aperture (2*f)
linvec = linspace(-d/2,d/2,n);
for i = 1:1:n
    appgrid(:,i,2) = linvec; %change x values
    appgrid(i,:,3) = linvec; %change y values
end

% seemygrid(appgrid)

% Define the direction of the incoming ray at each gridpoint
rayunitvecs = raydirection(0,appgrid,paxis);

% seemyrays(rayunitvecs,appgrid) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate N0, M, and P for each gridpoint
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

% Calculate path lengths for each incident ray
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
end

% seemyroots(point2mirror,n) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemymirror(reflectpoint) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemyreflectvecs(reflectvec,reflectpoint,appgrid) %%%%%%%%%%%%%%%%%%%%%%%

% Find the path lengths to the reference sphere
% for this to work we have to pretend that the sphere vertex is at [0,0,0]
% Define the reference sphere
sphere.r = 1e-6; %m sphere radius
sphere.x = -20;
sphere.y = 0;
sphere.z = 0;
sphere.ecc = 0; %define sphere eccentricity
sphere.f = sphere.r; %define the focus of the sphere which we know is at the center 1 radius away
% sphere.vertex = sphere.center + [sphere.f;0;0];
sphere.paxis = paxis; %should be the same as the mirror axis
sphere.p = sphere.f*(1 + sphere.ecc);

% define N0 and M for the reference sphere
sphere.N0 = -sphere.p*sphere.paxis;
sphere.M = eye(3); %since e is zero then the principle axis doesnt matter

% move the reflect points such that they match our new CS where the vertex of the sphere is [0,0,0]
newreflectpoints = movemypoints_sphere(reflectpoint,sphere,0,0,0);

% find the path length of the reflected ray to the reference sphere
mirror2sphere = zeros(1,size(rayunitvecs,2));
spherepoints = zeros(size(rayunitvecs,1),size(rayunitvecs,2));
for ind = 1:1:size(rayunitvecs,2)
    % get P and i for the reflected vector
    PP = newreflectpoints(:,ind); %the origin point is the reflected point on the mirror
    i = reflectvec(:,ind); %the unit vector direction of the reflected ray is the unit vector of the incident ray
    
    % set up the equation to solve the polynomial
    one = i'*sphere.M*i;
    two = 2*i'*(sphere.M*PP + sphere.N0);
    three = PP'*(sphere.M*PP + 2*sphere.N0);
    temp = roots([one,two,three]);
    mirror2sphere(1,ind) = min(temp);
    
    % find the point where the ray hit the reference sphere
    spherepoints(:,ind) = newreflectpoints(:,ind) + mirror2sphere(ind)*reflectvec(:,ind);
    
    
end

% seemysphere(spherepoints) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Undo the coordinate transformation to figure out when it hits the reference sphere
newspherepoints = undomovemypoints_sphere(spherepoints,sphere,0,0,0);

% seemysphere(newspherepoints) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemyraypathstosphere(newspherepoints,reflectpoint,appgrid) %%%%%%%%%%%%%
end

function [point2mirror,mirror2sphere] = raycodecheckup_refsphere()
% Set up variables to make this easily modular
f = 20; %m focal length
ecc = 1; %eccentricity
% p0 = [0,0,0]; %vertex of the mirror
paxis = [-1;0;0]; %principle axis of the mirror
d = 10; %m aperture diameter
al = [-2*f,0,0]; %aperture location (center of aperture)

% define the grid of rays to simulate
n = 5; %number of gridpoints per side
appgrid = zeros(n,n,3);
appgrid(:,:,1) = al(1,1); %all the x coordinates are at the aperture (2*f)
linvec = linspace(-d/2,d/2,n);
for i = 1:1:n
    appgrid(:,i,2) = linvec; %change x values
    appgrid(i,:,3) = linvec; %change y values
end

% seemygrid(appgrid)

% Define the direction of the incoming ray at each gridpoint
rayunitvecs = raydirection(0,appgrid,paxis);

% seemyrays(rayunitvecs,appgrid) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate N0, M, and P for each gridpoint
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

% Calculate path lengths for each incident ray
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
end

% seemyroots(point2mirror,n) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemymirror(reflectpoint) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemyreflectvecs(reflectvec,reflectpoint,appgrid) %%%%%%%%%%%%%%%%%%%%%%%

% Find the path lengths to the reference sphere
% for this to work we have to pretend that the sphere vertex is at [0,0,0]
% Define the reference sphere
sphere.r = 1; %m sphere radius
sphere.x = -20;
sphere.y = 0;
sphere.z = 0;
sphere.ecc = 0; %define sphere eccentricity
sphere.f = sphere.r; %define the focus of the sphere which we know is at the center 1 radius away
% sphere.vertex = sphere.center + [sphere.f;0;0];
sphere.paxis = paxis; %should be the same as the mirror axis
sphere.p = sphere.f*(1 + sphere.ecc);

% define N0 and M for the reference sphere
sphere.N0 = -sphere.p*sphere.paxis;
sphere.M = eye(3); %since e is zero then the principle axis doesnt matter

% move the reflect points such that they match our new CS where the vertex of the sphere is [0,0,0]
newreflectpoints = movemypoints_sphere(reflectpoint,sphere,0,0,0);

% find the path length of the reflected ray to the reference sphere
mirror2sphere = zeros(1,size(rayunitvecs,2));
spherepoints = zeros(size(rayunitvecs,1),size(rayunitvecs,2));
for ind = 1:1:size(rayunitvecs,2)
    % get P and i for the reflected vector
    PP = newreflectpoints(:,ind); %the origin point is the reflected point on the mirror
    i = reflectvec(:,ind); %the unit vector direction of the reflected ray is the unit vector of the incident ray
    
    % set up the equation to solve the polynomial
    one = i'*sphere.M*i;
    two = 2*i'*(sphere.M*PP + sphere.N0);
    three = PP'*(sphere.M*PP + 2*sphere.N0);
    temp = roots([one,two,three]);
    mirror2sphere(1,ind) = min(temp);
    
    % find the point where the ray hit the reference sphere
    spherepoints(:,ind) = newreflectpoints(:,ind) + mirror2sphere(ind)*reflectvec(:,ind);
    
    
end

% seemysphere(spherepoints) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Undo the coordinate transformation to figure out when it hits the reference sphere
newspherepoints = undomovemypoints_sphere(spherepoints,sphere,0,0,0);

% seemysphere(newspherepoints) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemyraypathstosphere(newspherepoints,reflectpoint,appgrid) %%%%%%%%%%%%%
end

function [point2mirror,mirror2sphere,P] = changealph(alph)
% Set up variables to make this easily modular
f = 20; %m focal length
ecc = 1; %eccentricity
% p0 = [0,0,0]; %vertex of the mirror
paxis = [-1;0;0]; %principle axis of the mirror
d = 10; %m aperture diameter
al = [-2*f,0,0]; %aperture location (center of aperture)

% define the grid of rays to simulate
n = 50; %number of gridpoints per side
appgrid = zeros(n,n,3);
appgrid(:,:,1) = al(1,1); %all the x coordinates are at the aperture (2*f)
linvec = linspace(-d/2,d/2,n);
for i = 1:1:n
    appgrid(:,i,2) = linvec; %change x values
    appgrid(i,:,3) = linvec; %change y values
end

% seemygrid(appgrid)

% Define the direction of the incoming ray at each gridpoint
rayunitvecs = raydirection(alph,appgrid,paxis);

% seemyrays(rayunitvecs,appgrid) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate N0, M, and P for each gridpoint
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

% Calculate path lengths for each incident ray
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
end

% seemyroots(point2mirror,n) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemymirror(reflectpoint) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemyreflectvecs(reflectvec,reflectpoint,appgrid) %%%%%%%%%%%%%%%%%%%%%%%

% Find the path lengths to the reference sphere
% for this to work we have to pretend that the sphere vertex is at [0,0,0]
% Define the reference sphere
sphere.r = 1; %m sphere radius
sphere.x = -20;
sphere.y = 0;
sphere.z = 0;
sphere.ecc = 0; %define sphere eccentricity
sphere.f = sphere.r; %define the focus of the sphere which we know is at the center 1 radius away
% sphere.vertex = sphere.center + [sphere.f;0;0];
sphere.paxis = paxis; %should be the same as the mirror axis
sphere.p = sphere.f*(1 + sphere.ecc);

% define N0 and M for the reference sphere
sphere.N0 = -sphere.p*sphere.paxis;
sphere.M = eye(3); %since e is zero then the principle axis doesnt matter

% move the reflect points such that they match our new CS where the vertex of the sphere is [0,0,0]
newreflectpoints = movemypoints_sphere(reflectpoint,sphere,0,0,0);

% find the path length of the reflected ray to the reference sphere
mirror2sphere = zeros(1,size(rayunitvecs,2));
spherepoints = zeros(size(rayunitvecs,1),size(rayunitvecs,2));
for ind = 1:1:size(rayunitvecs,2)
    % get P and i for the reflected vector
    PP = newreflectpoints(:,ind); %the origin point is the reflected point on the mirror
    i = reflectvec(:,ind); %the unit vector direction of the reflected ray is the unit vector of the incident ray
    
    % set up the equation to solve the polynomial
    one = i'*sphere.M*i;
    two = 2*i'*(sphere.M*PP + sphere.N0);
    three = PP'*(sphere.M*PP + 2*sphere.N0);
    temp = roots([one,two,three]);
    mirror2sphere(1,ind) = min(temp);
    
    % find the point where the ray hit the reference sphere
    spherepoints(:,ind) = newreflectpoints(:,ind) + mirror2sphere(ind)*reflectvec(:,ind);
    
    
end

% seemysphere(spherepoints) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Undo the coordinate transformation to figure out when it hits the reference sphere
newspherepoints = undomovemypoints_sphere(spherepoints,sphere,0,0,0);

% seemysphere(newspherepoints) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemyraypathstosphere(newspherepoints,reflectpoint,appgrid) %%%%%%%%%%%%%
end

function [] = spotdiagram(alph)
% Set up variables to make this easily modular
f = 20; %m focal length
ecc = 1; %eccentricity
% p0 = [0,0,0]; %vertex of the mirror
paxis = [-1;0;0]; %principle axis of the mirror
d = 10; %m aperture diameter
al = [-2*f,0,0]; %aperture location (center of aperture)

% define the grid of rays to simulate
n = 50; %number of gridpoints per side
appgrid = zeros(n,n,3);
appgrid(:,:,1) = al(1,1); %all the x coordinates are at the aperture (2*f)
linvec = linspace(-d/2,d/2,n);
for i = 1:1:n
    appgrid(:,i,2) = linvec; %change x values
    appgrid(i,:,3) = linvec; %change y values
end

% seemygrid(appgrid)

% Define the direction of the incoming ray at each gridpoint
rayunitvecs = raydirection(alph,appgrid,paxis);

% seemyrays(rayunitvecs,appgrid) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate N0, M, and P for each gridpoint
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

% Calculate path lengths for each incident ray
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
end

% seemyroots(point2mirror,n) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemymirror(reflectpoint) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemyreflectvecs(reflectvec,reflectpoint,appgrid) %%%%%%%%%%%%%%%%%%%%%%%

% Find the path lengths to the plane of the focus for a spot diagram
% For this we can just use a multiplier defined by the x coordinate of the reflection point
% and the location of the focus plane to find the point where that vector would hit the plane
spot = zeros(size(reflectpoint,1),size(reflectpoint,2));
for ind = 1:1:size(reflectvec,2)
    % identify the point we are starting at
    p = reflectpoint(:,ind);
    
    % calculate the proper multiplier
    dtogo = -f - p(1,1);
    m = dtogo/reflectvec(1,ind);
    
    % calculate the point in the focus plane where that vector crosses
    spot(:,ind) = p + m*reflectvec(:,ind);
end

% plot the spot diagram
figure; hold on
plot(spot(2,:),spot(3,:),'bx')
title(['Spot Diagram for \alpha=',num2str(alph),'^o'])
xlabel('Y, m');ylabel('Z, m')
grid on
axis equal
axis square

[dist,indexes] = pdist2([spot(2,:)',spot(3,:)'],[spot(2,:)',spot(3,:)'],'euclidean','largest',1); %https://www.mathworks.com/matlabcentral/answers/349856-calculating-the-distance-between-every-two-points-of-a-list
[maxdist,ind2] = max(dist);
ind1 = indexes(ind2);

% plot the line of max distance
plot([spot(2,ind1),spot(2,ind2)],[spot(3,ind1),spot(3,ind2)],'mo-','LineWidth',2,'MarkerSize',9);
text(spot(2,ind2),spot(3,ind2),['Max Distance is ',num2str(maxdist),' m  '],'horizontalalignment','right')
end

function [point2mirror,mirror2sphere,P] = changerefsphere(deltax,deltay,deltaz)
% Set up variables to make this easily modular
f = 20; %m focal length
ecc = 1; %eccentricity
% p0 = [0,0,0]; %vertex of the mirror
paxis = [-1;0;0]; %principle axis of the mirror
d = 10; %m aperture diameter
al = [-2*f,0,0]; %aperture location (center of aperture)

% define the grid of rays to simulate
n = 50; %number of gridpoints per side
appgrid = zeros(n,n,3);
appgrid(:,:,1) = al(1,1); %all the x coordinates are at the aperture (2*f)
linvec = linspace(-d/2,d/2,n);
for i = 1:1:n
    appgrid(:,i,2) = linvec; %change x values
    appgrid(i,:,3) = linvec; %change y values
end

% seemygrid(appgrid)

% Define the direction of the incoming ray at each gridpoint
rayunitvecs = raydirection(0,appgrid,paxis);

% seemyrays(rayunitvecs,appgrid) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate N0, M, and P for each gridpoint
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

% Calculate path lengths for each incident ray
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
end

% seemyroots(point2mirror,n) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemymirror(reflectpoint) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemyreflectvecs(reflectvec,reflectpoint,appgrid) %%%%%%%%%%%%%%%%%%%%%%%

% Find the path lengths to the reference sphere
% for this to work we have to pretend that the sphere vertex is at [0,0,0]
% Define the reference sphere
sphere.r = 1; %m sphere radius
sphere.x = -20;
sphere.y = 0;
sphere.z = 0;
sphere.ecc = 0; %define sphere eccentricity
sphere.f = sphere.r; %define the focus of the sphere which we know is at the center 1 radius away
% sphere.vertex = sphere.center + [sphere.f;0;0];
sphere.paxis = paxis; %should be the same as the mirror axis
sphere.p = sphere.f*(1 + sphere.ecc);

% define N0 and M for the reference sphere
sphere.N0 = -sphere.p*sphere.paxis;
sphere.M = eye(3); %since e is zero then the principle axis doesnt matter

% move the reflect points such that they match our new CS where the vertex of the sphere is [0,0,0]
newreflectpoints = movemypoints_sphere(reflectpoint,sphere,deltax,deltay,deltaz);

% find the path length of the reflected ray to the reference sphere
mirror2sphere = zeros(1,size(rayunitvecs,2));
spherepoints = zeros(size(rayunitvecs,1),size(rayunitvecs,2));
for ind = 1:1:size(rayunitvecs,2)
    % get P and i for the reflected vector
    PP = newreflectpoints(:,ind); %the origin point is the reflected point on the mirror
    i = reflectvec(:,ind); %the unit vector direction of the reflected ray is the unit vector of the incident ray
    
    % set up the equation to solve the polynomial
    one = i'*sphere.M*i;
    two = 2*i'*(sphere.M*PP + sphere.N0);
    three = PP'*(sphere.M*PP + 2*sphere.N0);
    temp = roots([one,two,three]);
    mirror2sphere(1,ind) = min(temp);
    
    % find the point where the ray hit the reference sphere
    spherepoints(:,ind) = newreflectpoints(:,ind) + mirror2sphere(ind)*reflectvec(:,ind);
    
    
end

% seemysphere(spherepoints) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Undo the coordinate transformation to figure out when it hits the reference sphere
newspherepoints = undomovemypoints_sphere(spherepoints,sphere,deltax,deltay,deltaz);

% seemysphere(newspherepoints) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seemyraypathstosphere(newspherepoints,reflectpoint,appgrid) %%%%%%%%%%%%%
end