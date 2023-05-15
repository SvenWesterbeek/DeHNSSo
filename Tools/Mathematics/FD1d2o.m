function [ DD ] = FD1d2o( D,d )
% This function determines the first derivatives of a function/field that
% is given on a uniform grid. Note: the derivative is taken in the along
% the rows of the data, assuming the corresponding coordinate ascends when
% the column index decreases. 
%
%       <-- x;  x_2 > x_1
%      _ _ _ _ _ _ _ _ _ _
%     | f(x_2) ... f(x_1) |
% D = |                   |
%     :                   :

% set size of derivative matrix DD
DD = 0*D; 

% using cental differences for the "center nodes":
DD(:,  (1+1):(end-1)) =  (-1/2 * D(:,  (1+2):(end-0))...
                      +      0 * D(:,  (1+1):(end-1))...
                      +    1/2 * D(:,  (1+0):(end-2)))/  d;

% using forward differences for the first (smallest coordinate value)  
% column:
DD(:,   end)          =  (-3/2 * D(:,        (end-0))...
                      +      2 * D(:,        (end-1))...
                      -    1/2 * D(:,        (end-2)))/  d;
        
% using backward differences for the last (largest coordinate value)  
% column:
DD(:,   1)            =  (-3/2 * D(:,          (1+0))...
                      +      2 * D(:,          (1+1))...
                      -    1/2 * D(:,          (1+2)))/(-d);                  

end

