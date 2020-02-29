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
Nr = linspace(3,20,100);
eta = linspace(0.05,0.5,100);
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

% Functions:
function [f] = natfreq(EE,rho,A,Nr,d,eta)
f = (3^0.25*(1)*sqrt(EE)*sqrt(d)) / (3*pi*sqrt(192)*(2*Nr+1)^(5/2)*sqrt(A)*sqrt(rho)*sqrt(1+eta));
end