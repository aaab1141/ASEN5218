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
testbars = [1 2;1 4;1 3;1 5];

test.case1.A = equilibrium_matrix(testnodes,testbars);

testnodes2 = [0 0 0 1 1 1;
              1 0 0 1 1 1;
              0 1 0 0 0 1;
              1 1 0 0 0 1];
testbars2 = [1 3;3 4;4 1;4 2];

test.case2.A = equilibrium_matrix(testnodes2,testbars2);

testnodes3 = [0 0 0 1 1 1;
              1 0 0 1 1 1;
              0 1 0 0 0 1;
              1 1 0 0 0 1];
testbars3 = [1 3;3 4;4 1;4 2;2 3];

test.case3.A = equilibrium_matrix(testnodes3,testbars3);

testnodes4 = [0 0 0 1 1 1;
              1 0 0 1 1 1;
              0 1 0 0 0 1;
              1 1 0 0 0 1];
testbars4 = [1 3;3 4;4 2];

test.case4.A = equilibrium_matrix(testnodes4,testbars4);
%}

% Problem 2
% Equilateral Triangle case 1
triangle.nodes = [0   0         0 0 0 0;
          1   0         0 0 0 0;
          1/2 sqrt(3)/2 0 0 0 0];
triangle.bars = [1 2;2 3;3 1];

triangle.A = equilibrium_matrix(triangle.nodes,triangle.bars);
[triangle.U,triangle.V,triangle.W] = svd(triangle.A);

% Regular Tetrahedron case 2
tet.nodes = [0   0         0   0 0 0;
          1   1         0   0 0 0;
          1/2 sqrt(3)/2 0   0 0 0;
          1/2 sqrt(3)/4 3/4 0 0 0];
tet.bars = [1 2;2 4;4 1;4 3;3 1;3 2];

tet.A = equilibrium_matrix(tet.nodes,tet.bars);

[tet.U,tet.V,tet.W] = svd(tet.A);

% Problem 3
a = 0; %case where alpha=0 degrees
twobay.a0.nodes = [.5  .5  0 1 1 1;
         -.5 .5  0 1 1 1;
         -.5 -.5 0 1 1 1;
         .5  -.5 0 1 1 1;
         (cos(a)*.5-sin(a)*.5)   (sin(a)*.5+cos(a)*.5)   1 0 0 0;
         (cos(a)*-.5-sin(a)*.5)  (sin(a)*-.5+cos(a)*.5)  1 0 0 0;
         (cos(a)*-.5-sin(a)*-.5) (sin(a)*-.5+cos(a)*-.5) 1 0 0 0;
         (cos(a)*.5-sin(a)*-.5)  (sin(a)*.5+cos(a)*-.5)  1 0 0 0;
         .5  .5  2 0 0 0;
         -.5 .5  2 0 0 0;
         -.5 -.5 2 0 0 0;
         .5  -.5 2 0 0 0];
twobay.a0.bars = [3 6;3 7;2 5;2 6;
        1 8;1 5;4 7;4 8;
        7 12;7 11;6 11;6 10;
        5 10;5 9;8 9;8 12;
        8 7;7 6;6 5;5 8;
        12 11;11 10;10 9;9 12];

twobay.a0.A = equilibrium_matrix(twobay.a0.nodes,twobay.a0.bars);
[twobay.a0.U,twobay.a0.V,twobay.a0.W] = svd(twobay.a0.A);
plotmytruss(twobay.a0.nodes,twobay.a0.bars,'Two-bay Truss','m')

a = pi/4; %case where alpha-45 or -45 degrees
twobay.a45.nodes = [.5  .5  0 1 1 1;
         -.5 .5  0 1 1 1;
         -.5 -.5 0 1 1 1;
         .5  -.5 0 1 1 1;
         (cos(a)*.5-sin(a)*.5)   (sin(a)*.5+cos(a)*.5)   1 0 0 0;
         (cos(a)*-.5-sin(a)*.5)  (sin(a)*-.5+cos(a)*.5)  1 0 0 0;
         (cos(a)*-.5-sin(a)*-.5) (sin(a)*-.5+cos(a)*-.5) 1 0 0 0;
         (cos(a)*.5-sin(a)*-.5)  (sin(a)*.5+cos(a)*-.5)  1 0 0 0;
         .5  .5  2 0 0 0;
         -.5 .5  2 0 0 0;
         -.5 -.5 2 0 0 0;
         .5  -.5 2 0 0 0];
twobay.a45.bars = twobay.a0.bars;
     
twobay.a45.A = equilibrium_matrix(twobay.a45.nodes,twobay.a45.bars);
[twobay.a45.U,twobay.a45.V,twobay.a45.W] = svd(twobay.a45.A);
plotmytruss(twobay.a45.nodes,twobay.a0.bars,'Two-bay Truss','m')

% Tracking th singluar values over variation in alpha
alph = -45:.1:45;
bars = [3 6;3 7;2 5;2 6;
        1 8;1 5;4 7;4 8;
        7 12;7 11;6 11;6 10;
        5 10;5 9;8 9;8 12;
        8 7;7 6;6 5;5 8;
        12 11;11 10;10 9;9 12];
% evol = zeros(size(bars,1),length(alph));
for i = 1:length(alph)
    nodes = [.5  .5  0 1 1 1;
         -.5 .5  0 1 1 1;
         -.5 -.5 0 1 1 1;
         .5  -.5 0 1 1 1;
         (cos(alph(i))*.5-sin(alph(i))*.5)   (sin(alph(i))*.5+cos(alph(i))*.5)   1 0 0 0;
         (cos(alph(i))*-.5-sin(alph(i))*.5)  (sin(alph(i))*-.5+cos(alph(i))*.5)  1 0 0 0;
         (cos(alph(i))*-.5-sin(alph(i))*-.5) (sin(alph(i))*-.5+cos(alph(i))*-.5) 1 0 0 0;
         (cos(alph(i))*.5-sin(alph(i))*-.5)  (sin(alph(i))*.5+cos(alph(i))*-.5)  1 0 0 0;
         .5  .5  2 0 0 0;
         -.5 .5  2 0 0 0;
         -.5 -.5 2 0 0 0;
         .5  -.5 2 0 0 0];
     A = equilibrium_matrix(nodes,bars);
     evol(:,i) = svd(A);
     ms(:,i) = size(A)'-rank(A);
end

figure
plot(alph,evol)
title('Variation of Singular Values over \alpha')
xlabel('\alpha, ^o');ylabel('Singular Value');grid on

figure
plot(alph,ms)
title('Number of Mechanisms and Surface Stresses over \alpha');grid on
xlabel('\alpha, ^o');ylabel('Count'),legend('Mechanisms','Surface Stresses')
    
% Problem 4

% Problem 5


% End of Script