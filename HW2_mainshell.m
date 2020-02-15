% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This is the main shell script for ASEN 5218 Homework #2. It facilitates
% running the problems in the homework. To check this homework, the only 
% thing that is required is to run this script.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear
close all

% Question 1 Plot
figure
hold on
theta = pi/40:.01:pi/2;
for n = 2:1:16
    temp = (sin(pi/n)/n)^(1/2) * (5*sqrt(2*cot(theta)) + sqrt(5./sin(theta)).*1./cos(theta));
    plot(theta,temp)
end
axis([0,pi/2,0,50]); grid on
title('$\mu_n \mu_{\theta}$ for Increasing Number of Longerons','interpreter','latex')
ylabel('\mu_n \mu_{\theta}');xlabel('\theta, radians')
legend('n = 2','n = 3','n = 4','n = 5','n = 6','n = 7','location','northwest')

% Question 2

% Question 3
% Looking to plot the linear density of both truss designs
E = 110e9; %Pa
rho = 1500; %kg/m^3
thickness_min = 0.5e-3; %m
num_longerons = [3 6 12];
radius = [0.1 0.3 1]; %m
loadP_min = 10; %N
loadP_max = 10000; %N


% Question 4 is a derivation

% Question 5
