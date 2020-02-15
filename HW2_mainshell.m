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
num_longerons = 2:1:7;
num_longerons = num_longerons';
for n = 2:1:7
    temp = (sin(pi/n)/sqrt(n)) * ((10*cot(theta)) + sqrt(10)./sin(theta).*1./cos(theta));
    [value_of_mu(n-1),i] = min(temp);
    optimal_theta(n-1,1) = theta(i);
    plot(theta,temp)
end
plot(optimal_theta',value_of_mu,'r','linewidth',2);
axis([0,pi/2,0,50]); grid on
title('$\mu_n \mu_{\theta}$ for Increasing Number of Longerons','interpreter','latex')
ylabel('\mu_n \mu_{\theta}');xlabel('\theta, radians')
legend('n = 2','n = 3','n = 4','n = 5','n = 6','n = 7','Minimum \theta','location','north')

% table of min thetas
table(num_longerons,optimal_theta)
% 
% % Question 2
% figure
% hold on
% theta = pi/40:.01:pi/2;
% num_longerons = 2:1:7;
% num_longerons = num_longerons';
% for n = 2:1:7
%     temp = sin(pi/n)^(2/3)/n^(1/3) * (5*2^(2/3)*cot(theta) + 
%     [value_of_mu(n-1),i] = min(temp);
%     optimal_theta(n-1,1) = theta(i);
%     plot(theta,temp)
% end
% plot(optimal_theta',value_of_mu,'r','linewidth',2);
% axis([0,pi/2,0,50]); grid on
% title('$\mu_n \mu_{\theta}$ for Increasing Number of Longerons','interpreter','latex')
% ylabel('\mu_n \mu_{\theta}');xlabel('\theta, radians')
% legend('n = 2','n = 3','n = 4','n = 5','n = 6','n = 7','Minimum \theta','location','north')
% 
% % table of min thetas
% table(num_longerons,optimal_theta)

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
