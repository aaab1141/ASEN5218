function [] = plotmytruss(nodes,bars,plottitle,units)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% This function generates a 3D plot of the truss structure based on the
% nodes and bars information.
% 
% INPUTS:
% nodes = j x 6 matrix where j is the number of nodes in the truss. Each
%         row represents one node and columns 1,2,3 are the x,y,z
%         coordinate and columns 4,5,6 are either 0 or 1 where 0 indicates
%         that DOF is free and 1 indicates that DOF is fixed. The DOF is not
%         used in the plotting.
% bars  = b x 2 matrix where b is the number of bars. The first column is
%         the first node number of bar b and the second column is the
%         second node number of bar b.
% title = string that you want the title of the plot to be.
% units = string that you want the units of the plot to be.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

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
    plot3(x,y,z)
    
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

