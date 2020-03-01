% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This is the main shell script for ASEN 5218 Homework #3. It facilitates
% running the problems in the homework. To check this homework, the only 
% thing that is required is to run this script.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear
close all

% Question 1
theta = 0:0.01:2*pi;
%   a) triangle-on-triangle lattice
R = 2/3*(cos(theta).^4 + sin(theta).^4);
figure
polarplot(theta,R)
title('Bending Efficiency: Square-on-Square 2 Layer Lattice')

%   b) square-on-square lattice
R = 2/3 ; R = repelem(R,length(theta));
figure
polarplot(theta,R)
title('Bending Efficiency: Triangle-on-Triangle 2 Layer Lattice')

% Question 2
EE = 105e9; %Pa
rho = .0015/100^3; %kg/m^3
d = 20; %m
A = 0.005^2*pi; %m Seems like a reasonable value since it's not given
Nr = linspace(3,20,50);
eta = linspace(0.05,0.5,50);
nf = zeros(length(Nr),length(eta));
for i = 1:length(Nr)
    for j = 1:length(eta)
        nf(i,j) = natfreq(EE,rho,A,Nr(i),d,eta(j));
    end
end
[Nr,eta] = meshgrid(Nr,eta);

figure
contour(Nr,eta,nf)
xlabel('Number of Rings, N_r');ylabel('Parasitc Mass Ratio, \eta')
colorbar; grid on
title('Contour Plot of Natural Frequency')

figure
surf(Nr,eta,nf)
xlabel('Number of Rings, N_r');ylabel('Parasitc Mass Ratio, \eta');zlabel('Natural Frequency, Hz')
title('Surface Plot of Natural Frequency')

% Question 3
% assume same values as Q2 for E and A
theta_0 = pi/4;
% choose reasonable values for d (thickness), l (length of bars), and r_0
d = .01; %10cm
l = .01; %10cm
r_0 = 1;
oor = linspace(1/100,10,50);
theta = linspace(0,pi,50);
for i = 1:length(oor)
    for j = 1:length(theta)
        U(i,j) = strainenergy(d,A,EE,l,theta(i),theta_0,oor(j),r_0);
    end
end
[oor,theta] = meshgrid(oor,theta);

figure
contour(oor,theta,U)
xlabel('1/r');ylabel('\theta, radians')
title('Countour Plot of Strain Energy')
colorbar; grid on

figure
surf(oor,theta,U)
xlabel('1/r');ylabel('\theta, radians');zlabel('Strain Energy, J')
title('Surface Plot of Strain Energy')

% Functions:
function [f] = natfreq(EE,rho,A,Nr,d,eta)
f = (3^0.25*(1)*sqrt(EE)*sqrt(d)) / (3*pi*sqrt(192)*(2*Nr+1)^(5/2)*sqrt(A)*sqrt(rho)*sqrt(1+eta));
end

function [u] = strainenergy(d,A,EE,l,theta,theta_0,oor,r_0)
u = d^2*A*EE/(4*l)*(oor^(2)*(cos(theta)^4+cos(theta+pi/2)^4) - r_0^(-2)*(cos(theta_0)^4+cos(theta_0+pi/2)^4));
end