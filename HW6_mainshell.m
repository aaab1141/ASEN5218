% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% This is the main shell script for ASEN 5218 Homework #6. It facilitates
% running the problems in the homework. To check this homework, the only
% thing that is required is to run this script.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear
close all
clc

%% Question 1
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
displs = u(:,end);

% Set the displacements of the fixed nodes to zero and add to displacmeents vector
displs = [zeros(size(displs));displs];

% Build the B matrix
B = makeB(displs,bars);

% get little g
g = B*stresses;

GtD = g'*displs;
disp(['GtD = ',num2str(round(GtD,4),'%f'),' which means we have an inextensional mechanism that cannot be stabilized.'])


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