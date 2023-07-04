function [StabRes] = IC_HNS(j,i,Re,Ur,Wr,dyUr,dyWr,D1,D2,StabRes,StabGrid)
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
%   The IC_HNS code calls ILST (that finds a solution to the local 
%   Orr-Sommerfeld problem). The result is normalized and stored in
%   StabRes.

%   Note: In HNS, ILST should only be used to introduce modes at linear 
%   amplitudes at the inflow therefore no nonlinear introduction methods 
%   are included here. To initialize nonlinear shape functions, please use 
%   the Stab.IC = "LOAD" option.

%% Overview of inputs
% Name              size               unit explanation
% i                 (1,1)               [-] Streamwise station
% j                 (1,1)               [-] Mode counter
% Re                (1,1)               [-] Reynolds number (global)
% Ur                (ny,nx)             [-] Streamwise velocity profile of the base flow
% Wr                (ny,nx)             [-] Spanwise velocity profile of the base flow
% dyUr              (ny,nx)             [-] Wall-normal derivative of the streamwise velocity profile of the base flow
% dyWr              (ny,nx)             [-] Wall-normal derivative of the spanwise velocity profile of the base flow
% D1                (4*ny, 4*ny)        [-] First-order spectral derivative operator
% D2                (4*ny, 4*ny)        [-] Second-order spectral derivative operator
% StabRes.omegavec  ((2*N+1)*(2*M+1),1) [-] Omega per mode
% StabRes.betavec   ((2*N+1)*(2*M+1),1) [-] Beta per mode

%% Overview of outputs
% Name          size      unit    explanation
% StabRes.u     (ny,nx)   [-]     Streamwise perturbation velocity profile
% StabRes.v     (ny,nx)   [-]     Wall-normal perturbation velocity profile
% StabRes.w     (ny,nx)   [-]     Spanwise perturbation velocity profile
% StabRes.p     (ny,nx)   [-]     Perturbation pressure profile
% StabRes.phi   (4*ny,nx) [-]     Perturbation vector of state variables (u,v,w,p)^T

%% pre define some values from inputs

ny = length(StabGrid.etaun);
omega = StabRes.omegavec(j);
beta = StabRes.betavec(j);

%% Find initial condition as solution to local stability eigenvalue problem

% Solve local eigenvalue problem ILST
fprintf(' Calculating IC for nonzero modes with ILST: ')
[eig,u,v,w,p] = solver_ILST(Re,Ur(:,i),Wr(:,i),dyUr(:,i),dyWr(:,i),omega,beta,ny,D1,D2);

% Find d99 for filtering purposes
k = find(Ur(:,i)<=0.99*max(Ur(:,i)),1,'first');   
d99(i)=StabGrid.etaun(k);

% Filter eigenfunctions for correct eigenvalue
[alpha,index]   = EVfilter(eig,v,StabGrid.etaun,omega,beta,d99);



%% Normalize output by max(u)

% Calculate new maximum amplitude after spline interpolation
% Select correct u
u                 = u(:,index);

% Calculate amplitude
y_mdInt           = linspace(StabGrid.etaun(1),StabGrid.etaun(end),4000);
u1int             = interp1(StabGrid.etaun,abs(u),y_mdInt,'spline');
[ampltd, ~]       = max(abs(u1int));

% Normalize result by u_max
StabRes.u(j,:,i)  = u/ampltd;
StabRes.v(j,:,i)  = v(:,index)/ampltd;
StabRes.w(j,:,i)  = w(:,index)/ampltd;
StabRes.p(j,:,i)  = p(:,index)/ampltd;

% Create phi from eigenvectors
StabRes.phi(j,:,i) = ([ StabRes.u(j,:,i).';
                        StabRes.v(j,:,i).';
                        StabRes.w(j,:,i).';
                        StabRes.p(j,:,i).']);

%% Visual Check for Mode Selection

% Open figure
figure
tiledlayout(1,2);

% Eigenvaluespectrum plot
nexttile
plot(eig,'k.')
axis([-2 2 -2 2]) % make a function of the eig or alpha
hold on
plot(alpha,'ro')
grid on
xlabel('Real','FontName','Times New Roman','FontSize',10)
ylabel('Imaginary','FontName','Times New Roman','FontSize',10)
legend('Eigenvalues','Selected eigenvalue','Location','NorthOutside','FontName','Times New Roman','FontSize',10)
sgtitle('Inflow LST solution: eigenvalue and shape function','FontName','Times New Roman','FontSize',10)

% Shape function
nexttile
title('Shape function')
plot(abs(squeeze(StabRes.u(j,:,1))),StabGrid.etaun)
xlabel('$\hat{u}$','Interpreter','Latex','FontName','Times New Roman','FontSize',10)
ylabel('\eta','FontName','Times New Roman','FontSize',10)
grid on
set(gcf, 'Position', [1000 50 700 400]);
pause(0.05)
end



