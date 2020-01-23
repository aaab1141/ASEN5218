function [A] = equilibrium_matrix(nodes,bars)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This function creates the equilibrium matrix A from information given
% in two matrices, [nodes] and [bars].
% 
% INPUTS:
% nodes = j x 6 matrix where j is the number of nodes in the truss. Each
%         row represents one node and columns 1,2,3 are the x,y,z
%         coordinate and columns 4,5,6 are either 0 or 1 where 0 indicates
%         that DOF is free and 1 indicates that DOF is fixed.
% bars  = b x 2 matrix where b is the number of bars. The first column is
%         the first node number of bar b and the second column is the
%         second node number of bar b.
% 
% OUTPUTS: 
% A     = Equilibrium Matrix
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Get the size of the A matrix
nr = size(nodes,1)*3; %number of rows in A is going to be three per node given in the nodes matrix
nc = size(bars,1); %number of bars is the number of bars in the bars matrix

% Allocate the A matrix
A = zeros(nr,nc);

% Go through each bar and get the information into the A matrix
for b = 1:1:nc %b is the bar number
    % Node geometery of bar
    fromnode = bars(b,1); 
    tonode = bars(b,2);
    
    % Coordinates of nodes
    fromnodepos = nodes(fromnode,1:3);
    tonodepos = nodes(tonode,1:3);
    
    % Unit vector of the bar
    barunitvec = (tonodepos' - fromnodepos');
    barunitvec = barunitvec/norm(barunitvec,2);
    
    % Add the contribution of each bar to the appropriate node locations in A
    A(b*3-2:b*3,fromnode) = A(b*3-2:b*3,fromnode) + barunitvec;
    A(b*3-2:b*3,tonode) = A(b*3-2:b*3,tonode) + barunitvec;
end
end

