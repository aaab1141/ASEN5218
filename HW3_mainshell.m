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
R = 3/4 ; R = repelem(R,length(theta));
figure
polarplot(theta,R)
title('Bending Efficiency: Triangle-on-Triangle 2 Layer Lattice')

% Question 2