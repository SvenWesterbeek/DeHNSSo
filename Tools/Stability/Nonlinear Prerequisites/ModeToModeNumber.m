function [mode] = ModeToModeNumber(n,m,N,M)
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
% Calculates the mode number for a given (n,m) for a truncation given by (N,M)

%% Inputs
% Name  size    unit explanation
% M     (1,1)   [-]  Modal truncation in t
% N     (1,1)   [-]  Modal truncation in z
% m     (1,1)   [-]  Mode number in t
% n     (1,1)   [-]  Mode number in z

%% Error check
if n>N || m>M
    error('Input exceeds truncation (m,n) > (M,N)')
end

%% mint
L = (2*N+1)*(2*M+1); %Counter

%initialize basic building blocks (m,n)
M1 = -M:1:M; 
N1 = -N.*ones(1,2*N+1); 

% Set up omega vector
for i = 1:L 
 Mvec(i) = M1(mod(i-1,length(M1))+1);
end

% Set up beta vector
for k = 1:L/(2*M+1)
 Nvec([(k-1)*(2*M+1)+1:k*(2*M+1)]) = -N+(k-1)*ones(1,2*M+1);
end

% Find mode number
[~,mindex] = find(Mvec==m);
[~,nindex] = find(Nvec==n);
mode = intersect(mindex,nindex);


end

