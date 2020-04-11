% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% This is the main shell script for ASEN 5218 Homework #4. It facilitates
% running the problems in the homework. To check this homework, the only
% thing that is required is to run this script.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear
close all

%% Question 1
disp('Question 1')
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
b1projlengths = linspace(0.01,9.99,1000); % projected length of b1
b1lengths = zeros(1,length(b1projlengths));
b2lengths = b1lengths;
b3lengths = b1lengths;
b4lengths = b1lengths;
defori = b1lengths;
definn = b1lengths;
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
    a = asin(b2/2/b1);
    dorigin = 2*pi - 12*a;
    defori(i) = dorigin;
    
    % calculate phi
    phi = acos(b2/2/b1);
    
    % calculate b3
    b3 = nodedist(innnodes(:,2),outnodes(:,4));
    b3lengths(i) = b3;
    
    % calculate psi
    psi = acos(b2/2/b3);
    
    % calculate b4
    b4 = nodedist(innnodes(:,1),outnodes(:,1));
    b4lengths(i) = b4;
    
    % calculate theta with law of cosines
    theta = acos((b5^2-b4^2-b3^2)/(-2*b4*b3));
    
    % calculate angular defect of the inner ring joints
    dinner = 2*pi - 2*(phi+theta+psi);
    definn(i) = dinner;
end

% Find the intersection point
[~,minin] = min(abs(defori-definn));
angledef = defori(minin);
disp('The bar lengths should be:')
disp(['b1 = ',num2str(b1lengths(minin)),' m | b2 = ',num2str(b2lengths(minin)),' m | b3 = ',num2str(b3lengths(minin)),' m | b4 = ',num2str(b4lengths(minin)),' m | b5 = ',num2str(b5)]);

figure
plot(b1lengths,defori,b1lengths,definn)
title('Origin and Inner Nodes Angular Defect')
xlabel('Length of b1, m');ylabel('Angular Defect, rad')
grid on
legend('Origin node','Inner Circle Node','location','northwest')
save_fig_png('Q1 Angular Defect')

% plot the reflector
[nodes,bars] = makereflect(b1projlengths(minin),outnodes,numinnnodes,inalph);
plotmytruss(nodes,bars,'Reflector Truss','m')
save_fig_png('Q1 Reflector')

%% Question 3
disp('Question 3')
qr = 1350; %w/m^2 radiation from the sun
a = .9; %absortivity
e = .8; %emissivity
c = 700; %J/K/kg specifig heat
alph = 7e-6; %m/m/K coefficient of thermal expansion
rho = 1800; %kg/m^3 density of bars
E = 50e9; %Pa modulus of elasticity
L0 = .5; %m length of battens
r = .01; %m radius of solid cylinder bars
sig = 5.670374419e-8; %Stephan - Boltzman constant

% Find T0
T0 = (a*2*qr/sig/e/pi)^(1/4);
disp(['T0 = ',num2str(round(T0,2))])

% Transient response
tspan = [0,1000]; %s
sol = ode45(@(t,T) -2*sig*e*T^4/rho/c/r,tspan,T0);
t = 1:1:1000;
T = deval(sol,t);

figure
plot(t,T,t,T0*ones(size(T)))
title('Temperature of Truss Bars After Maneuver')
xlabel('Time, s');ylabel('Temperature, K')
grid on
legend('L_{left}, L_{batt} & L_{diag}','L_{right}','location','east')
save_fig_png('Transient_Response')

% calculate the change in lengths of the left bars
L = L0*(1+alph*(T-T0));
figure
plot(t,L)
title('Length of L_{left} After Maneuver')
xlabel('Time, s');ylabel('L_{left}, m')
grid on
save_fig_png('Length_of_L')

l = L0*(1+alph*(deval(sol,[10,100,1000])-T0));

disp(['Length at 10s = ',num2str(l(1)),' m'])
disp(['Length at 100s = ',num2str(l(2)),' m'])
disp(['Length at 1000s = ',num2str(l(3)),' m'])

% Estimate the curvature based on the length of the bars
roc = L*L0./(L0-L)-L/2;
figure
plot(t,roc)
title('Radius of Curvature After Maneuver')
xlabel('Time, s');ylabel('Radius of Curvature, m')
set(gca,'Yscale','log')
set(gca,'Xscale','log')
grid on
save_fig_png('Radius_of_Curvature')

ll = l*L0./(L0-l)-l/2;

disp(['Radius of Curvature at 10s = ',num2str(ll(1)),' m'])
disp(['Radius of Curvature at 100s = ',num2str(ll(2)),' m'])
disp(['Radius of Curvature at 1000s = ',num2str(ll(3)),' m'])

% N = 10; %number of bays

% Functions
function [nodes,bars] = makereflect(b1proj,outnodes,numinnnodes,inalph)
% calcualte inner circle nodes
for n = 0:1:numinnnodes-1
    innnodes(:,n+1) = getnode(b1proj*cos(n*inalph),b1proj*sin(n*inalph));
end
% make the nodes matrix
nodes = [outnodes,innnodes,[0;0;0]]';
% make the bars matrix
bars = [1,2;2,3;3,4;4,5;5,6;6,7;7,8;8,9;9,10;10,11;11,12;12,1;...
    13,1;13,2;14,2;14,3;14,4;15,4;15,5;15,6;16,6;16,7;16,8;...
    17,8;17,9;17,10;18,10;18,11;18,12;13,12;13,14;14,15;15,16;16,17;17,18;18,13;...
    13,19;14,19;15,19;16,19;17,19;18,19];   
end

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