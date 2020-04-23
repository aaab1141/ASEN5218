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

% +pi/4 structure
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


