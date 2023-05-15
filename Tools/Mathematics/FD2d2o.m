function [ DD ] = FD2d2o( D,d )
% This function determines the second order derivatives of a function/field 
% that is given on a uniform grid. Note: the derivative is taken in the 
% along the rows of the data, assuming the corresponding coordinate ascends 
% when the column index decreases. 
%
%       <-- x;  x_2 > x_1
%      _ _ _ _ _ _ _ _ _ _
%     | f(x_2) ... f(x_1) |
% D = |                   |
%     :                   :

% set size of derivative matrix DD
DD = 0*D; 

% using cental differences for the "center nodes":
DD(:,  (1+1):(end-1)) =(    D(:,(1+2):(end-0))...
                      - 2 * D(:,(1+1):(end-1))...
                      + 1 * D(:,(1+0):(end-2)))./d^2;

% using forward differences for the first (smallest coordinate values) two 
% columns:
DD(:,        (end-0)) =(2 * D(:,      (end-0))...
                      - 5 * D(:,      (end-1))...
                      + 4 * D(:,      (end-2))...
                      - 1 * D(:,      (end-3)))./d^2;
        
% using backward differences for the last (largest coordinate values) two 
% columns:
DD(:,  (1+0)        ) =(2 * D(:,(1+0)        )...
                      - 5 * D(:,(1+1)        )...
                      + 4 * D(:,(1+2)        )...
                      - 1 * D(:,(1+3)        ))./d^2;                  

end

