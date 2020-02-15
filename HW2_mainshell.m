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
    [value_of_mu.cyl(n-1),i] = min(temp);
    optimal_theta.cyl(n-1,1) = theta(i);
    plot(theta,temp)
end
% table of min thetas
table(num_longerons,optimal_theta.cyl)

%finish plot
optimal_theta.cyl = [optimal_theta.cyl;optimal_theta.cyl(1)];
value_of_mu.cyl = [value_of_mu.cyl,50];
plot(optimal_theta.cyl',value_of_mu.cyl,'r','linewidth',2);
axis([0,pi/2,0,50]); grid on
title('Solid Cylinder $\mu_n \mu_{\theta}$ for Increasing Number of Longerons','interpreter','latex')
ylabel('\mu_n \mu_{\theta}');xlabel('\theta, radians')
legend('n = 2','n = 3','n = 4','n = 5','n = 6','n = 7','Minimum \theta','location','north')
h = text(optimal_theta.cyl(1)+.03,30,['Minimum \theta = ',num2str(optimal_theta.cyl(1)*180/pi),'^o'],'horizontalalignment','center');
set(h,'Rotation',90)
 
% Question 2
figure
hold on
theta = pi/40:.01:pi/2;
num_longerons = 2:1:7;
num_longerons = num_longerons';
for n = 2:1:7
    temp = sin(pi/n)^(2/3)/n^(1/3) * (5*2^(2/3)*cot(theta) + (2^(1/3)*5^(2/3))./(sin(theta).^(2/3).*cos(theta)));
    [value_of_mu.tube(n-1),i] = min(temp);
    optimal_theta.tube(n-1,1) = theta(i);
    plot(theta,temp)
end
% table of min thetas
table(num_longerons,optimal_theta.tube)

%finish plot
optimal_theta.tube = [optimal_theta.tube;optimal_theta.tube(1)];
value_of_mu.tube = [value_of_mu.tube,50];
plot(optimal_theta.tube',value_of_mu.tube,'r','linewidth',2);
axis([0,pi/2,0,50]); grid on
title('Solid Cylinder $\mu_n \mu_{\theta}$ for Increasing Number of Longerons','interpreter','latex')
ylabel('\mu_n \mu_{\theta}');xlabel('\theta, radians')
legend('n = 2','n = 3','n = 4','n = 5','n = 6','n = 7','Minimum \theta','location','north')
h = text(optimal_theta.tube(1)+.03,30,['Minimum \theta = ',num2str(optimal_theta.tube(1)*180/pi),'^o'],'horizontalalignment','center');
set(h,'Rotation',90)
 

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
