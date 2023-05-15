function [ DD ] = FD1d4o( D,d )
%% License

% Written by K.J. Groot
%% Description
% This function determines the first derivatives of a function/field that
% is given on a uniform grid. Note: the derivative is taken in the along
% the rows of the data, assuming the corresponding coordinate ascends when
% the column index decreases. A check is performed to see if enough data-
% points are present to perform the calculation.
%
%       <-- x;  x_2 > x_1
%      _ _ _ _ _ _ _ _ _ _
%     | f(x_2) ... f(x_1) |
% D = |                   |
%     :                   :

%% Error
if length(D(1,:)) < 6
    % .....< backward/forward difference on most inner node
    % >..... forward/backward difference on most inner node
    % ...
    % ^ ^ outer node done with central differences
    % | boundary node
    error('--- Function length too small; less than 6')
end

%% Perform differentiation
% set size of derivative matrix DD
DD = 0*D; 

% using cental differences for the "center nodes":
DD(:,  (1+2):(end-2)) =  (1/12 * D(:,  (1+4):(end-0))...
                      -    2/3 * D(:,  (1+3):(end-1))...
                      +      0 * D(:,  (1+2):(end-2))...
                      +    2/3 * D(:,  (1+1):(end-3))...
                      -   1/12 * D(:,  (1+0):(end-4)))/ d;

% using forward differences for the first (smallest coordinate values) two 
% columns:
DD(:,(end-1):(end-0)) =(-25/12 * D(:,(end-1):(end-0))...
                      +      4 * D(:,(end-2):(end-1))...
                      -      3 * D(:,(end-3):(end-2))...
                      +    4/3 * D(:,(end-4):(end-3))...
                      -    1/4 * D(:,(end-5):(end-4)))/ d;
        
% using backward differences for the last (largest coordinate values) two 
% columns:
DD(:,  (1+0):  (1+1)) =(-25/12 * D(:,  (1+0):  (1+1))...
                      +      4 * D(:,  (1+1):  (1+2))...
                      -      3 * D(:,  (1+2):  (1+3))...
                      +    4/3 * D(:,  (1+3):  (1+4))...
                      -    1/4 * D(:,  (1+4):  (1+5)))/(-d);                  

end

