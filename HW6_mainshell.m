% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% This is the main shell script for ASEN 5218 Homework #6. It facilitates
% running the problems in the homework. To check this homework, the only
% thing that is required is to run this script.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear
close all

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

% Build G and D matrices to test the kind of mechanism we have


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