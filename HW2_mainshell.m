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
h = text(optimal_theta.cyl(1)+.03,30,['Minimum \theta = ',num2str(optimal_theta.cyl(1)*180/pi,3),'^o'],'horizontalalignment','center');
set(h,'Rotation',90)
 
% Question 2 Plot
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
title('Hollow Tube $\mu_n \mu_{\theta}$ for Increasing Number of Longerons','interpreter','latex')
ylabel('\mu_n \mu_{\theta}');xlabel('\theta, radians')
legend('n = 2','n = 3','n = 4','n = 5','n = 6','n = 7','Minimum \theta','location','north')
h = text(optimal_theta.tube(1)+.03,30,['Minimum \theta = ',num2str(optimal_theta.tube(1)*180/pi,3),'^o'],'horizontalalignment','center');
set(h,'Rotation',90)

% Question 3
E = 110e9; %Pa
rho = 1500; %kg/m^3
t_min = 0.5e-3; %m minimum thickness of the tube
loadP_min = .1; %N
loadP_max = 100; %N
theta_solid = optimal_theta.cyl(1); %radians
theta_tube = optimal_theta.tube(1); %radians

% Plot M/L as a function of the load P for fied parameters
P = loadP_min:1:loadP_max;
figure; hold on; grid on
xlabel('Load P, N');ylabel('M/L, kg/m')
title('Le Title')
for n = [3 6 12]
    for R = [0.1 0.3 1] %m global truss radius
        % For solid cylinder
        K.cyl = 2/5/sqrt(pi);
        muR.cyl = R;
        muP.cyl = sqrt(P);
        muM.cyl = rho/E^(1/2);
        muN.cyl = sin(pi/n)/sqrt(n);
        muTheta.cyl = 10*cot(theta_solid) + sqrt(10)/(sin(theta_solid)*cos(theta_solid));
        store.cyl.(['longeron',num2str(n)]).(['R',num2str(10*R)]) = K.cyl*muR.cyl*muP.cyl*muM.cyl*muN.cyl*muTheta.cyl;
        plot(P,store.cyl.(['longeron',num2str(n)]).(['R',num2str(10*R)]),'linewidth',2,'color','g')
        text(.7*loadP_max,store.cyl.(['longeron',num2str(n)]).(['R',num2str(10*R)])(round(.7*loadP_max)),['Cyl R=',num2str(R),' n=',num2str(n)],'horizontalalignment','center')
        
        % For hollow tube
        K.tube = 1/5;
        muR.tube = R^(2/3);
        muP.tube = P.^(1/3);
        muM.tube = rho*t_min^(2/3)/E^(1/3);
        muN.tube = sin(pi/n)^(2/3)/n^(1/3);
        mutheta_tube.tube = (5*2^(2/3)*cot(theta_tube) + (2^(1/3)*5^(2/3))/(sin(theta_tube)^(2/3)*cos(theta_tube)));
        store.tube.(['longeron',num2str(n)]).(['R',num2str(10*R)]) = K.tube*muR.tube*muP.tube*muM.tube*muN.tube*mutheta_tube.tube;
        plot(P,store.tube.(['longeron',num2str(n)]).(['R',num2str(10*R)]),'linewidth',2,'color','b')
        text(.5*loadP_max,store.tube.(['longeron',num2str(n)]).(['R',num2str(10*R)])(round(.5*loadP_max)),['tube R=',num2str(R),' n=',num2str(n)],'horizontalalignment','center')
        
        % Find the value of M/L where these two lines cross
        [store.intersects.(['longeron',num2str(n)]).(['R',num2str(10*R)]).x,store.intersects.(['longeron',num2str(n)]).(['R',num2str(10*R)]).y] = polyxpoly(P,store.cyl.(['longeron',num2str(n)]).(['R',num2str(10*R)]),P,store.tube.(['longeron',num2str(n)]).(['R',num2str(10*R)]));
        
        
    end
end

% Question 4 is a derivation

% Question 5
% For a poisson ratio of 1/3 we derived the formulation in class for bars
% at an arbitrary angle alpha, so we use that here.
% Square Lattice: alpha=0 and alpha=90
disp('Square Lattice')
a = 0;
A1 = [cos(a)^2; sin(a)*cos(a); sin(a)^2] * [cos(a)^2 sin(a)*cos(a) sin(a)^2]
a = pi/2;
A2 = [cos(a)^2; sin(a)*cos(a); sin(a)^2] * [cos(a)^2 sin(a)*cos(a) sin(a)^2]
A = A1 + A2

% Polar Plot
theta = 0:0.01:2*pi;
R = 1/2*(cos(theta).^4 + sin(theta).^4);
figure
polarplot(theta,R)
title('Equivalent Modulus for a Square Lattice')

% Triangle Lattice: alpha=0 and alpha=60 and alpha=120
disp('Triangle Lattice')
a = 0;
A1 = [cos(a)^2; sin(a)*cos(a); sin(a)^2] * [cos(a)^2 sin(a)*cos(a) sin(a)^2]
a = pi/3;
A2 = [cos(a)^2; sin(a)*cos(a); sin(a)^2] * [cos(a)^2 sin(a)*cos(a) sin(a)^2]
a = 2*pi/3;
A3 = [cos(a)^2; sin(a)*cos(a); sin(a)^2] * [cos(a)^2 sin(a)*cos(a) sin(a)^2]
A = A1 + A2 + A3

% Polar Plot
theta = 0:0.01:2*pi;
R = 2.25*sqrt(3)/6; R = repelem(R,length(theta));
figure
polarplot(theta,R)
title('Equivalent Modulus for a Equilateral Triangle Lattice')