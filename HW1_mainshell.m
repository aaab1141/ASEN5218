% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This is the main shell script for ASEN 5218 Homework #1. It facilitates
% running the problems in the homework. For each problem a single script
% file is used. To check this homework, the only thing that is required is
% to run this script.
% 
% Written 2020-01-23 | Aaron Aboaf
% Modified 2020-01-23 | Aaron Aboaf
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear
close all

% Run Problem 1
% Build a function that creates the A matrix from nodes and bars matrices.
testnodes = [0  0  1  0 0 0;
             -1 0  0  1 1 1;
             0  -1 0  1 1 1;
             1  0  0  1 1 1;
             0  1  0  1 1 1];
testbars = [1 2;1 3;1 4;1 5];

A = equilibrium_matrix(testnodes,testbars);

testnodes2 = [0 0 0 1 1 1;
              1 0 0 1 1 1;
              0 1 0 0 0 1;
              1 1 0 0 0 1];
testbars2 = [1 3;3 4;4 1;4 2];

A = equilibrium_matrix(testnodes2,testbars2);

testnodes3 = [0 0 0 1 1 1;
              1 0 0 1 1 1;
              0 1 0 0 0 1;
              1 1 0 0 0 1];
testbars3 = [1 3;3 4;4 1;4 2;2 3];

A = equilibrium_matrix(testnodes3,testbars3);

testnodes4 = [0 0 0 1 1 1;
              1 0 0 1 1 1;
              0 1 0 0 0 1;
              1 1 0 0 0 1];
testbars4 = [1 3;3 4;4 2];

A = equilibrium_matrix(testnodes4,testbars4);

% Run Problem 2


% Run Problem 3


% Run Problem 4


% Run Problem 5


% End of Script