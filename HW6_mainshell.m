% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% This is the main shell script for ASEN 5218 Homework #6. It facilitates
% running the problems in the homework. To check this homework, the only
% thing that is required is to run this script.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear
close all
clc

%% Question 1
disp("Question 1")
% -pi/4 structure
a = -pi/4;
nodes = [.5  .5  0 1 1 1;
    -.5 .5  0 1 1 1;
    -.5 -.5 0 1 1 1;
    .5  -.5 0 1 1 1;
    (cos(a)*.5-sin(a)*.5)   (sin(a)*.5+cos(a)*.5)   1 0 0 0;
    (cos(a)*-.5-sin(a)*.5)  (sin(a)*-.5+cos(a)*.5)  1 0 0 0;
    (cos(a)*-.5-sin(a)*-.5) (sin(a)*-.5+cos(a)*-.5) 1 0 0 0;
    (cos(a)*.5-sin(a)*-.5)  (sin(a)*.5+cos(a)*-.5)  1 0 0 0];
bars = [3 6;3 7;2 5;2 6;
    1 8;1 5;4 7;4 8;
    8 7;7 6;6 5;5 8];

plotmytruss(nodes,bars,'Structure 1','m')
save_fig_png('Q1.Structure 1')

% Get the equilibrium matrix and do the svd
A = equilibrium_matrix(nodes,bars);
[u,~,w] = svd(A);
r = rank(A);
m = size(u,2)-r;
s = size(w',1)-r;
disp('Structure: alpha = -pi/4')
disp(['Rank = ',num2str(r),' | m = ',num2str(m),' | s = ',num2str(s)])

stresses = w(:,end);
highlightmytruss(nodes,bars,'Structure 1 Bar Stresses','m',stresses)
save_fig_png('Q1.Structure 1 Tensegrity')
displs = u(:,end);

% Set the displacements of the fixed nodes to zero and add to displacmeents vector
displs = [zeros(size(displs));displs];

% Build the B matrix
B = makeB(displs,bars);

% get little g
g = B*stresses;

GtD = g'*displs;
disp(['GtD = ',num2str(GtD),' which is positive thus we have an infinitesimal mechanism that can be stabilized.'])


%% +pi/4 structure
a = pi/4;
nodes = [.5  .5  0 1 1 1;
    -.5 .5  0 1 1 1;
    -.5 -.5 0 1 1 1;
    .5  -.5 0 1 1 1;
    (cos(a)*.5-sin(a)*.5)   (sin(a)*.5+cos(a)*.5)   1 0 0 0;
    (cos(a)*-.5-sin(a)*.5)  (sin(a)*-.5+cos(a)*.5)  1 0 0 0;
    (cos(a)*-.5-sin(a)*-.5) (sin(a)*-.5+cos(a)*-.5) 1 0 0 0;
    (cos(a)*.5-sin(a)*-.5)  (sin(a)*.5+cos(a)*-.5)  1 0 0 0];
bars = [3 6;3 7;2 5;2 6;
    1 8;1 5;4 7;4 8;
    8 7;7 6;6 5;5 8];

plotmytruss(nodes,bars,'Structure 2','m')
save_fig_png('Q1.Structure 2')

% Get the equilibrium matrix and do the svd
A = equilibrium_matrix(nodes,bars);
[u,~,w] = svd(A);
r = rank(A);
m = size(u,2)-r;
s = size(w',1)-r;
disp('Structure: alpha = +pi/4')
disp(['Rank = ',num2str(r),' | m = ',num2str(m),' | s = ',num2str(s)])

stresses = w(:,end);
highlightmytruss(nodes,bars,'Structure 2 Bar Stresses','m',stresses)
save_fig_png('Q1.Structure 2 Tensegrity')
displs = u(:,end);

% Set the displacements of the fixed nodes to zero and add to displacmeents vector
displs = [zeros(size(displs));displs];

% Build the B matrix
B = makeB(displs,bars);

% get little g
g = B*stresses;

GtD = g'*displs;
disp(['GtD = ',num2str(round(GtD,4),'%f'),' which means we have an inextensional mechanism that cannot be stabilized.'])

%% Question 2
disp("Question 2")
warning('off','MATLAB:ode45:IntegrationTolNotMet')
% A = 1.5;
Amin = .5; Amax = 4;
% H = 2;
Hmin = .5; Hmax = 3;
% rc = .4;
rcmin = .2; rcmax = .5;

% opts = odeset('RelTol',1e-3,'AbsTol',1e-6); %default
% [z,drdz] = ode45(@(z,drdz) balloon(z,drdz,A,H,rc),[0,H],[.02;10000],opts);

disp('If we don''t set a threshold on the throat radius:')
divs = 20; rtthresh = 0;
allA = linspace(Amin,Amax,divs);
allH = linspace(Hmin,Hmax,divs);
allrc = linspace(rcmin,rcmax,divs);
V = zeros(divs,divs,divs);
tic
for i = 1:1:divs
    for j = 1:1:divs
        for k = 1:1:divs
            V(i,j,k) = balloonvolumeonly(allrc(i),allA(j),allH(k),rtthresh);
        end
    end
end
toc
[~,in] = min(V(:));
[i,j,k] = ind2sub(size(V),in);

[z,r,rt,l,vol] = balloonvolume(allrc(i),allA(j),allH(k));

figure
plot(z,r)
xlabel('z coordinate, m');ylabel('Radius, m')
title(['Optimal Balloon Profile if Rt > ',num2str(rtthresh)])
axis equal; grid on
disp(['Balloon Throat = ',num2str(rt),' m | Balloon Arclength = ',num2str(l),' m | Crown Radius = ',num2str(r(1)),' m'])
disp(['Balloon Volume = ',num2str(vol),' m^3'])
save_fig_png('Q2.Balloon small rt')

n = 10;
[tphi,ttheta] = balloonstress(z,r,allA(j),allH(k),n);

% find the maximum stress
[maxstress,inms] = max(tphi);

figure
plot(z,tphi,'k',z,ttheta,'--r')
title('Pressure Differential Normalized Stresses Along Balloon Surface')
xlabel('Z Position Along Balloon Axis, m');ylabel('Stress')
legend('T_\phi','T_\theta') 
text(z(inms),tphi(inms),['\leftarrow Maximum Stress at z = ',num2str(round(z(inms),3)),' r = ',num2str(round(r(inms),4))])
grid on
save_fig_png('Q2.Balloon small rt stress')

rtthresh = .010;
disp(['If we set a threshold that the throat radius >',num2str(rtthresh)])
divs = 20; 
allA = linspace(Amin,Amax,divs);
allH = linspace(Hmin,Hmax,divs);
allrc = linspace(rcmin,rcmax,divs);
V = zeros(divs,divs,divs);
tic
for i = 1:1:divs
    for j = 1:1:divs
        for k = 1:1:divs
            V(i,j,k) = balloonvolumeonly(allrc(i),allA(j),allH(k),rtthresh);
        end
    end
end
toc
[~,in] = min(V(:));
[i,j,k] = ind2sub(size(V),in);

[z,r,rt,l,vol] = balloonvolume(allrc(i),allA(j),allH(k));

figure
plot(z,r)
xlabel('z coordinate, m');ylabel('Radius, m')
title(['Optimal Balloon Profile if Rt > ',num2str(rtthresh)])
axis equal; grid on
disp(['Balloon Throat = ',num2str(rt),' m | Balloon Arclength = ',num2str(l),' m | Crown Radius = ',num2str(r(1)),' m'])
disp(['Balloon Volume = ',num2str(vol),' m^3'])
save_fig_png('Q2.Balloon big rt')

n = 10;
[tphi,ttheta] = balloonstress(z,r,allA(j),allH(k),n);

% find the maximum stress
[maxstress,inms] = max(tphi);

figure
plot(z,tphi,'k',z,ttheta,'--r')
title('Pressure Differential Normalized Stresses Along Balloon Surface')
xlabel('Z Position Along Balloon Axis, m');ylabel('Stress')
legend('T_\phi','T_\theta')
text(z(inms),tphi(inms),['\leftarrow Maximum Stress at z = ',num2str(round(z(inms),3)),' r = ',num2str(round(r(inms),4))])
grid on
save_fig_png('Q2.balloon big rt stress')

disp('There is no optimal balloon profile above a Rt = ~.0105 m')

%% Question 4
theta = linspace(0,pi/2,100);
a = 1; b = a;

for gam = [pi/3,pi/4,2*pi/3]
    H = a*sin(gam).*sin(theta);
    twoL = 2*sqrt((a*cos(gam))^2 + (a*sin(gam).*cos(theta)).^2);
    squiggly = tan((a*sin(gam).*cos(theta))./(a*cos(gam)));
    twoS = 2*cos(pi/2 - squiggly);
    
    figure
    plot(theta,H,theta,twoL,theta,twoS)
    xlabel('\Theta, radians');ylabel('Length, units')
    title(['Lengths assuming a = b & \gamma = ',num2str(round(gam,2))])
    legend('H','2L','2S')
    grid on
    save_fig_png(['Q4.gam is ',num2str(round(gam,2))])
end

%% Functions
function [] = highlightmytruss(nodes,bars,plottitle,units,stresses)
% Create a figure to work on
figure

% turn the hold on so that we can plot each bar individually
hold on

% loop through the bars and plot them all
for i = 1:size(bars,1)
    % barnum = i;
    firstnodenum = bars(i,1);
    secondnodenum = bars(i,2);
    x = [nodes(firstnodenum,1);nodes(secondnodenum,1)];
    y = [nodes(firstnodenum,2);nodes(secondnodenum,2)];
    z = [nodes(firstnodenum,3);nodes(secondnodenum,3)];
    if stresses(i) > 0
        plot3(x,y,z,'r') %plot cable elements in tension red
    else
        plot3(x,y,z,'k') %plot solid elements in compression black
    end
    
    % Add the bar number to the bar
    % find the middle of the bar
    midx = mean(x);
    midy = mean(y);
    midz = mean(z);
    
    % Add the text on the plot
    text(midx,midy,midz,num2str(i),'fontsize',10,'fontangle','italic','color','b')
end

% Add the node numbers
for i = 1:size(nodes,1)
    text(nodes(i,1),nodes(i,2),nodes(i,3),num2str(i),'fontsize',12,'fontweight','bold','color','k')
end

% Label the plot and change the viewing angle
title(plottitle)
xlabel(['X, ',units])
ylabel(['Y, ',units])
zlabel(['Z, ',units])
view(45,35.264)
axis equal
grid on
end

function [B] = makeB(displs,bars)
% reshape the mechanisms vector to work with the algo
displs = reshape(displs,3,length(displs)/3)';

% Get the size of the B matrix
nr = size(displs,1)*3; %there are 3 DOF for each node (3 equations)
nc = size(bars,1); %there is a column for each bar because there is a bar force

% Allocate the A matrix
B = zeros(nr,nc);

% Go through each bar and get the information into the A matrix
for b = 1:1:size(bars,1) %b is the bar number
    % Node geometery of bar
    fromnode = bars(b,1);
    tonode = bars(b,2);
    
    % Coordinates of nodes
    fromnodepos = displs(fromnode,1:3);
    tonodepos = displs(tonode,1:3);
    
    % Unit vector of the bar
    barunitvec = (tonodepos' - fromnodepos');
    barunitvec = barunitvec/norm(barunitvec,2);
    
    % Add the contribution of each bar to the appropriate node locations in A
    % The appropriate row location is a set of rows corresponding to the
    % two nodes attached to the particular bar force, and the appropriate
    % column location is the column associated with the bar force in question
    B(fromnode*3-2:fromnode*3,b) = B(fromnode*3-2:fromnode*3,b) - barunitvec;
    B(tonode*3-2:tonode*3,b) = B(tonode*3-2:tonode*3,b) + barunitvec;
end
end

function [drdz] = balloon(z,drdz,A,H,rc)
r = drdz(1);
rp = drdz(2); %rprime

r1 = rp;
r2 =((1+rp^2)^(3/2))/H*((H*(1-A*sqrt(r/H)))/(r*sqrt(1+rp^2)) - 2/rc*(H-z)*exp(2*A*sqrt(r/H)));

drdz = [r1;r2];
end

function [l] = getarclength(x,y)
l = 0;
for i = 2:1:length(x)
    l = l+norm([x(i-1);y(i-1);0]-[x(i);y(i);0]);
end
end

function [z,r,rt,l,vol] = balloonvolume(rc,A,H)
% integrate balloon shape
[z,drdz] = ode45(@(z,drdz) balloon(z,drdz,A,H,rc),[0,H],[.02;10000]);

% check if its valid by searching for that local minimum
lm = islocalmin(drdz(:,1)); %find the local mins of the balloon profile to find the throat
copylm = lm;
lm(lm == 0) = [];
nummin = length(lm);

if nummin > 1
    % find the index of the local minima
    minindex = find(copylm == 1);
    
    % crop to the first minimum
    z = z(1:minindex(1),1);
    r = drdz(1:minindex(1),1);
    
    % find the throat radius
    rt = r(end);
    
    % find the arclength
    l = getarclength(z,r);
    
    % set the valid flag so that volume is calculated
    validflag = 1;
else
    % we have an invalid balloon shape
    validflag = 0;
end

% determine volume
if validflag == 1
    vol = 0;
    % calculate the volume numerically by method of rings
    for i = 2:1:length(z)
        rmid = (r(i-1) + r(i))/2; %midpoint radius of ring
        thickness = z(i) - z(i-1); %thickness of ring
        vol = vol + pi*rmid^2*thickness; %add volume of ring
    end
else
    % say everything is zero
    vol = 0;
    l = 0;
    rt = 0;
end
end

function [vol] = balloonvolumeonly(rc,A,H,rtthresh)
% integrate balloon shape
[z,drdz] = ode45(@(z,drdz) balloon(z,drdz,A,H,rc),[0,H],[.02;10000]);

% check if its valid by searching for that local minimum
checkreal = isreal(drdz(:,1));
if checkreal == 1
    lm = islocalmin(drdz(:,1)); %find the local mins of the balloon profile to find the throat
    copylm = lm;
    lm(lm == 0) = [];
    nummin = length(lm);
    
    if nummin > 1
        % find the index of the local minima
        minindex = find(copylm == 1);
        
        % crop to the first minimum
        z = z(1:minindex(1),1);
        r = drdz(1:minindex(1),1);
        
        % find the throat radius
        rt = r(end);
        
        % find the arclength
        l = getarclength(z,r);
        
        % set the valid flag so that volume is calculated
        validflag = 1;
    else
        % we have an invalid balloon shape
        validflag = 0;
    end
    
    % determine volume
    if validflag == 1
        if l > .27
            vol = 0;
        elseif rt < rtthresh
            vol = 0;
        else
            vol = 0;
            % calculate the volume numerically by method of rings
            for i = 2:1:length(z)
                rmid = (r(i-1) + r(i))/2; %midpoint radius of ring
                thickness = z(i) - z(i-1); %thickness of ring
                vol = vol + pi*rmid^2*thickness; %add volume of ring
            end
        end
    else
        % say everything is zero
        vol = 0;
        %     l = 0;
        %     rt = 0;
    end
else
    vol = 0;
end
vol = -vol; %make it negative so the biggest volume is the minumum
end

function [Tphi,Ttheta] = balloonstress(z,r,A,H,n)
dp = 1;
for i = 1:1:length(z)
    R = r(i);
    ra = R*H;
    C = dp/(H+z(i))*H*ra/2;
    Tphi(i) = C*exp(-A*R^n/n);
    Ttheta(i) =  C*(1-A*R^n)*exp(-A*R^n/n);
end
end