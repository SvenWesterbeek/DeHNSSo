function [Nmat, Mmat, Modevec,Mvec, Nvec] = Mint(M,N)
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
% mint (Mode INTeractions) is a function that sets up a matrix containing the
% modes affected by the interactions of all currently represented modes

% Outputs Mmat and Nmat are mode numbers for beta and omega respectively.
% Output Modevec is the mode vector [m(1:2M+1);n(1:2M+1)] for every mode
% counter l_{m,n}.

%% Inputs
% Name  size    unit explanation
% M     (1,1)   [-]  Modal truncation in t
% N     (1,1)   [-]  Modal truncation in z

%% Outputs
% Name    size                      unit explanation
% Mmat    (M+1)*(N+1) x (M+1)*(N+1) [-]  Contains the m part of the modes affected by the mode interactions
% Nmat    (M+1)*(N+1) x (M+1)*(N+1) [-]  Contains the n part of the modes affected by the mode interactions
% Modevec (M+1)*(N+1) x 2           [-]  Contains the m and n of the modecounter l
% Mvec    (M+1)*(N+1) x 1           [-]  Contains the m of the modecounter l
% Nvec    (M+1)*(N+1) x 1           [-]  Contains the n of the modecounter l
        

%% mint
L = (2*N+1)*(2*M+1); %Counter

%initialize basic building block
M1 = -M:1:M; 

% Set up omega vector
for i = 1:L 
 Mvec(i) = M1(mod(i-1,length(M1))+1);
end

% Set up beta vector
for k = 1:L/(2*M+1)
 Nvec((k-1)*(2*M+1)+1:k*(2*M+1)) = -N+(k-1)*ones(1,2*M+1);
end
Modevec = [Mvec;Nvec];

% Find recipient mode for interaction of mode i with mode j
for i = 1:L %Multiplying powers is addition in power
    for j = 1:L
        Nmat(i,j) = Nvec(i)+Nvec(j); 
    end
end
for i = 1:L %Multiplying powers is addition in power
    for j = 1:L
        Mmat(i,j) = Mvec(i)+Mvec(j);
    end
end
end

