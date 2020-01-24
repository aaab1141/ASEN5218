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
nr = size(nodes,1)*3; %there are 3 DOF for each node (3 equations)
nc = size(bars,1); %there is a column for each bar because there is a bar force

% Allocate the A matrix
A = zeros(nr,nc);

% Go through each bar and get the information into the A matrix
for b = 1:1:size(bars,1) %b is the bar number
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
    % The appropriate row location is a set of rows corresponding to the
    % two nodes attached to the particular bar force, and the appropriate
    % column location is the column associated with the bar force in question
    A(fromnode*3-2:fromnode*3,b) = A(fromnode*3-2:fromnode*3,b) - barunitvec;
    A(tonode*3-2:tonode*3,b) = A(tonode*3-2:tonode*3,b) + barunitvec;    
end

% Find rows that correspond with kinematic DOF constraints (info from nodes)
kinc = zeros(size(A,1),1);
for i = 1:1:size(nodes,1)
    kinc(i*3-2:i*3,1) = nodes(i,4:6)';
end

% Remove redundant rows that are locked by the constraints on the nodes
A(kinc == 1,:) = [];

end

