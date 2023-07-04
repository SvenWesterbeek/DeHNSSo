function [ HB ] = Hbalancing(M, N, Mmat, Nmat, Modevec )
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2023 S.H.J. Westerbeek, S. Hulshoff, H. Schuttelaars
%                          & M. Kotsonis
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%     Please cite this code and the paper if you have used this for your
%     publication:

% DOI:
% DOI:

%% Description
% Hbalancing creates a 3D matrix that will be used to calculate
% the source terms for the ILPSE solver. In each matrix HB(:,:,j), a 1 
% means that the modes corresponding to the row and column index of that 
% entry interact to force that mode j. For example, for N = 2, M = 0, 
% HB(:,:,3) corresponding to the MFD will look like:

%    n,m -2,0   -1,0    0,0     1,0     2,0        
%  n,m    __________________________________
% -2,0   | 0      0      0       0       1
% -1,0   | 0      0      0       1       0
%  0,0   | 0      0      1       0       0
% -1,0   | 0      1      0       0       0
% -2,0   | 1      0      0       0       0

% As the interactions (-2,0)x(2,0), (-1,0)x(1,0), (0,0)x(0,0), (1,0)x(-1,0
% and (2,0)x(-2,0) all force the MFD (0,0). In practice, the interaction
% (0,0)x(0x0) is excluded from forcing term calculations however.
%%  Inputs
% Name    size                      unit explanation
% M       (1,1)                     [-]  Number of modes in t
% N       (1,1)                     [-]  Number of modes in z
% Mmat    ((M+1)*(N+1),(M+1)*(N+1)) [-]  Contains the m part of the modes affected by the mode interactions
% Nmat    ((M+1)*(N+1),(M+1)*(N+1)) [-]  Contains the n part of the modes affected by the mode interactions
% Modevec ((M+1)*(N+1),1)           [-]  Contains the m and n of the modecounter l

%%  Outputs
%   HB (Harmonic Balancing) is a collection of [2M+1]x[2N+1] matrices of size [2M+1]x[2N+1]
%   stored as a 3D matrix where every matrix is a binary matrix that serves
%   as a registry for the source amplitudes that affect that respective mode.

% Name    size                                          unit explanation
% HB      ((2M+1)x(2N+1),(2M+1)x(2N+1),(2M+1)x(2N+1))   [-]  Harmonic balancing matrix

%% 
L = (2*N+1)*(2*M+1);
for i = 1:L
% Read desired combination of m,n for mode i (mode l)
m=Modevec(1,i); 
n=Modevec(2,i);
A = Mmat == m;      %A tells us which entry has the correct mode for m
B = Nmat == n;      %B tells us which entry has the correct mode for n
HB(:,:,i) = A.*B;   %Simulates LOGIC AND, given that both are required.
end 

end

