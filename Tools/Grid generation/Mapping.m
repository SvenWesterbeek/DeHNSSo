function [ ynd,D1y,D2y ] = Mapping(Grid,eta,D1eta,D2eta)
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
%     publication using:

% DOI:
% DOI:

if Grid.ytype ~= "step"
%% Mapping Malik (see Malik 1990)
% Transform Chebyshev to physical grid
ynd = (Grid.y_i*Grid.H*(1+eta)./(Grid.H - eta*(Grid.H - 2*Grid.y_i)))';

if nargout > 1
    % Compute scale factors w.r.t. differentiation
    detady   = (Grid.H - eta*(Grid.H - 2*Grid.y_i)).^2/(2*Grid.y_i*Grid.H*(Grid.H - Grid.y_i));
    d2etady2 = -2*detady.^2*(Grid.H - 2*Grid.y_i)./(Grid.H - eta*(Grid.H - 2*Grid.y_i));

    D1y = diag(detady)*D1eta;    
    D2y = diag(d2etady2)*D1eta + diag(detady.^2)*D2eta;
end

else
%% Stretched mapping (
% This function maps the Chebyshev grid on the domain [-1,0,1] to
% [0,Grid.y_i,Grid.H]. Furthermore, it delivers the pseudo-spectral 
% differentiation matrices that correspond to this transformation

%Correction to Grid.y_i for stretching
Grid.y_i = 0.5*0.43655*Grid.y_i; 
n = length(eta);

% Set Distribution of the first Constant part
dy1 = ones(1,floor(n/2)).*Grid.y_i/floor(n/2);

% Calculate the distance that is covered by the quadratic function
A3 = Grid.H-dy1(1)*n;

% Calculate alpha numerically
alpha = 0;
stepSize = 0.00001; % Adjust this value for finer discretization

for x = floor(n/2):stepSize:n
    alpha = alpha + (x - floor(n/2))^2;
end

alpha = A3 / alpha * (1/stepSize);


for j = 1:ceil(n/2)
    if j == 1
    dy2(j) = dy1(1);
    else
    dy2(j) = alpha*(j)^2+dy1(1);
    end
end

%% Stretch the distribution in various ways 
k = 2;
yndStretch = k*(Grid.y_i*Grid.H*(1+eta)./(Grid.H - eta*(Grid.H - 2*k*Grid.y_i)))';

k = 1;
yndOG = k*(Grid.y_i*Grid.H*(1+eta)./(Grid.H - eta*(Grid.H - 2*k*Grid.y_i)))';

% Define f as the stretch parameter
f= Grid.ystretch;

ynd = (1+eta')./(f+eta')*(f+1)/2*Grid.H;
ynd = ynd-ynd(end);
ynd = ynd/max(ynd)*Grid.H;

% Sum the various ynd to create a near constant delta y near the wall
ynd = (ynd + 20*yndOG + yndStretch)/22;

figure
title('Wall-normal collocation point distribution')
hold on
plot(ynd,'.')
grid on
xlabel('Collocation point','FontName','Times New Roman','FontSize',10)
ylabel('\eta','FontName','Times New Roman','FontSize',10)
set(gcf, 'Position', [1500 50 400 400]);

%% Calculate elements of the derivative matrix E (malik, 1990)
%Set up c (eq 3.22)
c = ones(1,n);
c(1) = 2;
c(end) = 2;

%Set up E (eq 3.23)
for j = 1:n
    for k = 1:n
       if k~=j
           E(j,k) = c(j)/c(k) * (-1)^(k+j) / (eta(j)-eta(k));
       else
           E(j,j) = - (eta(j) / (2*(1-eta(j)^2)));
       end
    end
end

E(1,1) = (2*n^2+1)/6; %
E(n,n) = -E(1,1);

% Calculate second derivative matrix
E2 = zeros(size(E));
for j = 1:n
    for k = 1:n
        for m = 1:n
E2(j,k) = E2(j,k) + E(j,m)*E(m,k);
        end
    end
end

%% Set up D1y and D2y
% Calculate the scaling factor for transformation physical -> comp
S = FD1d2o_uneven(ynd,eta')';

S2 = FD2d2o_uneven(ynd,eta')';

% Calculate first derivative matrix 
D1y = diag(S)*D1eta;

D2y = diag(S2)*D1eta + diag(S.^2)*D2eta;

end % grid mode
end % of function


