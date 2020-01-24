% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This is the main shell script for ASEN 5218 Homework #1. It facilitates
% running the problems in the homework. For each problem a single script
% file is used. To check this homework, the only thing that is required is
% to run this script.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear
close all

% Problem 1
% Build a function that creates the A matrix from nodes and bars matrices.
%{
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
%}

% Problem 2
% Equilateral Triangle case 1
nodes1 = [0   0         0 0 0 0;
          1   0         0 0 0 0;
          1/2 sqrt(3)/2 0 0 0 0];
bars1 = [1 2;2 3;3 1];

A1 = equilibrium_matrix(nodes1,bars1);

[U1,V1,W1] = svd(A1);

% Regular Tetrahedron case 2
nodes2 = [0   0         0   0 0 0;
          1   1         0   0 0 0;
          1/2 sqrt(3)/2 0   0 0 0;
          1/2 sqrt(3)/4 3/4 0 0 0];
bars2 = [1 2;2 4;4 1;4 3;3 1;3 2];

A2 = equilibrium_matrix(nodes2,bars2);

[U2,V2,W2] = svd(A2);

% Problem 3


% Problem 4


% Problem 5


% End of Script