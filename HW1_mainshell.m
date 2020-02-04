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
test.case1.nodes = [0  0  1  0 0 0;
             -1 0  0  1 1 1;
             0  -1 0  1 1 1;
             1  0  0  1 1 1;
             0  1  0  1 1 1];
test.case1.bars = [1 2;1 4;1 3;1 5];

test.case1.A = equilibrium_matrix(test.case1.nodes,test.case1.bars);
plotmytruss(test.case1.nodes,test.case1.bars,'Test Case 1 Truss','m')

test.case2.nodes = [0 0 0 1 1 1;
              1 0 0 1 1 1;
              0 1 0 0 0 1;
              1 1 0 0 0 1];
test.case2.bars = [1 3;3 4;4 1;4 2];

test.case2.A = equilibrium_matrix(test.case2.nodes,test.case2.bars);
plotmytruss(test.case2.nodes,test.case2.bars,'Test Case 2 Truss','m')

test.case3.nodes = [0 0 0 1 1 1;
              1 0 0 1 1 1;
              0 1 0 0 0 1;
              1 1 0 0 0 1];
test.case3.bars = [1 3;3 4;4 1;4 2;2 3];

test.case3.A = equilibrium_matrix(test.case3.nodes,test.case3.bars);

test.case4.nodes = [0 0 0 1 1 1;
              1 0 0 1 1 1;
              0 1 0 0 0 1;
              1 1 0 0 0 1];
test.case4.bars = [1 3;3 4;4 2];

test.case4.A = equilibrium_matrix(test.case4.nodes,test.case4.bars);

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

node_num = 1:1:12;
node_num = node_num';
node_num = repelem(node_num,3);
a = pi/4; %case where alpha = 45 degrees
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
plotmytruss(twobay.a45.nodes,twobay.a0.bars,'Twisted Two-bay Truss','m')

% make a table for the node displacments
node_num = 5:1:12;
node_num = node_num';
node_num = repelem(node_num,3);
direction = ['x';'y';'z';'x';'y';'z';'x';'y';'z';'x';'y';'z';'x';'y';'z';'x';'y';'z';'x';'y';'z';'x';'y';'z'];
disp_1 = twobay.a45.U(:,end);
p45distab = table(node_num,direction,disp_1);
disp('45 Degree Twist Node Displacements');disp(p45distab)

%make a table for the bar stresses
bar_num = 1:1:size(twobay.a0.bars,1);
bar_num = bar_num';
stress_1 = twobay.a45.W(end,:)';
p45strtab = table(bar_num,stress_1);
disp('45 Degree Twist Bar Stresses');disp(p45strtab)

a = -pi/4; %case where alpha = -45 degrees
twobay.na45.nodes = [.5  .5  0 1 1 1;
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
twobay.na45.bars = twobay.a0.bars;
     
twobay.na45.A = equilibrium_matrix(twobay.na45.nodes,twobay.na45.bars);
[twobay.na45.U,twobay.na45.V,twobay.na45.W] = svd(twobay.na45.A);

% make a table for the node displacments
node_num = 5:1:12;
node_num = node_num';
node_num = repelem(node_num,3);
direction = ['x';'y';'z';'x';'y';'z';'x';'y';'z';'x';'y';'z';'x';'y';'z';'x';'y';'z';'x';'y';'z';'x';'y';'z'];
disp_1 = twobay.na45.U(:,end-1);
disp_2 = twobay.na45.U(:,end);
n45distab = table(node_num,direction,disp_1,disp_2);
disp('-45 Degree Twist Node Displacements');disp(n45distab)

%make a table for the bar stresses
bar_num = 1:1:size(twobay.a0.bars,1);
bar_num = bar_num';
stress_1 = twobay.na45.W(end-1,:)';
stress_2 = twobay.na45.W(end,:)';
n45strtab = table(bar_num,stress_1,stress_2);
disp('-45 Degree Twist Bar Stresses');disp(n45strtab)

% Tracking the singluar values over variation in alpha
alph = -pi/4:pi/180:pi/4;
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
plot(alph*180/pi,evol,'-o')
title('Variation of Singular Values over \alpha')
xlabel('\alpha, ^o');ylabel('Singular Value');grid on

figure; hold on
plot(alph*180/pi,ms(1,:),'b-o','linewidth',2)
plot(alph*180/pi,ms(2,:),'r-o')
title('Number of Mechanisms and Surface Stresses over \alpha');grid on
xlabel('\alpha, ^o');ylabel('Count'),legend('Mechanisms','Surface Stresses')
    
% Problem 4
gr = (1+sqrt(5))/2;
iso.coordinates = [-gr/2 -.5 0;gr/2 -.5 0;gr/2  .5 0;-gr/2  .5 0;
         0 -gr/2  -.5;0  gr/2  -.5;0  gr/2   .5;0 -gr/2   .5;
         -.5 0 -gr/2;-.5 0  gr/2;.5 0  gr/2;.5 0 -gr/2];
iso.nodes = zeros(length(iso.coordinates),6);
a = pi/2 - atan(1/gr);
Rx = [1 0 0;0 cos(a) -sin(a); 0 sin(a) cos(a)];
for i = 1:length(iso.coordinates)
    iso.nodes(i,1:3) = (Rx*iso.coordinates(i,:)')';
end
iso.bars = [9 5;1 9;5 1;8 1;8 2;
            2 12;2 11;11 10;10 8;10 1;
            10 4;4 9;4 7;4 1;8 5;9 12;
            9 6;4 6;12 3;12 5;12 6;5 2;
            2 3;11 3;3 6;3 7;11 7;10 7;7 6;11 8];
         
plotmytruss(iso.nodes,iso.bars,'Icosahedron','m')
                   
iso.A = equilibrium_matrix(iso.nodes,iso.bars);
[iso.U,iso.V,iso.W] = svd(iso.A);
sum(iso.A,1)
size(iso.A)' - rank(iso.V)'

% Problem 5
line.nodes = [0 0 1 0 0 0; 0 0 2 0 0 0; 0 0 3 0 0 0; 0 0 4 0 0 0];
line.bars = [1 2;2 3;3 4;4 3];
line.A = equilibrium_matrix(line.nodes,line.bars);
plotmytruss(line.nodes,line.bars,'Line','m')

% End of Script