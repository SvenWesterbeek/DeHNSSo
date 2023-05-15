%% Description
% This example caller can be used to perform linear and nonlinear 
% simulations of Tollmien-Schlichting waves in a Blasius boundary layer
% using 

% ██████████            █████   █████ ██████   █████  █████████   █████████          
%░░███░░░░███          ░░███   ░░███ ░░██████ ░░███  ███░░░░░███ ███░░░░░███         
% ░███   ░░███  ██████  ░███    ░███  ░███░███ ░███ ░███    ░░░ ░███    ░░░   ██████ 
% ░███    ░███ ███░░███ ░███████████  ░███░░███░███ ░░█████████ ░░█████████  ███░░███
% ░███    ░███░███████  ░███░░░░░███  ░███ ░░██████  ░░░░░░░░███ ░░░░░░░░███░███ ░███
% ░███    ███ ░███░░░   ░███    ░███  ░███  ░░█████  ███    ░███ ███    ░███░███ ░███
% ██████████  ░░██████  █████   █████ █████  ░░█████░░█████████ ░░█████████ ░░██████ 
%░░░░░░░░░░    ░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░  ░░░░░░░░░   ░░░░░░░░░   ░░░░░░  
                                                                                          
% Banner created using: https://manytools.org/hacker-tools/ascii-banner/     

%% License (GNU GENERAL PUBLIC LICENSE v3)
%                    Delft Harmonic Navier-Stokes Solver
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

% Paper DOI:
% Code  DOI:

%% Clean up and load paths
clear all
close all

%Load paths
addpath(genpath('../../..'))

%% BF: Base Flow data, base flow grid and reference values

% Name         size  unit    explanation
% BF.X   (nxbl,1)    [-]     Base Flow grid streamwise locations 
% BF.Y   (1,nybl)    [-]     Base Flow grid streamwise locations
% BF.U   (nxbl,nybl) [-]     Base Flow Streamwise Velocity
% BF.V   (nxbl,nybl) [-]     Base Flow Wall-normal Velocity
% BF.W   (nxbl,nybl) [-]     Base Flow Spanwise Velocity
% BF.dxU (nxbl,nybl) [-]     Streamwise Gradient of Base Flow Streamwise Velocity
% BF.dxV (nxbl,nybl) [-]     Streamwise Gradient of Base Flow Wall-normal Velocity
% BF.dxW (nxbl,nybl) [-]     Streamwise Gradient of Base Flow Spanwise Velocity
% BF.dyU (nxbl,nybl) [-]     Wall-normal Gradient of Base Flow Streamwise Velocity
% BF.dyV (nxbl,nybl) [-]     Wall-normal Gradient of Base Flow Wall-normal Velocity
% BF.dyW (nxbl,nybl) [-]     Wall-normal Gradient of Base Flow Spanwise Velocity

% BF.lref   (1)      [m]     reference length (Blasius lengthscale)
% BF.Uref   (1)      [m/s]   reference velocity 
% BF.nu     (1)      [m^2/s] kinematic viscocity 
% BF.Re     (1)      [-]     Reynolds number

% Load Blasius boundary layer base flow and reference values
load('BF_Blasius.mat','BF')

%% Grid: Numerical domain specifications OR Numerical domain grid points  
% Name         size         units explanation
% Grid.nx       (1)         [-]   Number of streamwise stations of the numerical grid
% Grid.ny       (1)         [-]   Number of wall-normal collocation points of the numerical grid
% Grid.wall     (nxwall,2)  [-]   Wall definition x and y locations on arbitrary grid
% Grid.H        (1)         [-]   Domain height
% Grid.y_i      (1)         [-]   Median collocation point height
% Grid.S        (1)         [-]   Stability grid domain start
% Grid.L        (1)         [-]   Stability grid domain length (wall length)
% Grid.mode     (string)    [-]   Grid generation mode (see grid_gen)
%%%% equidistant    - an equidistant streamwise distribution. Wall-refined. 
%%%% xrefined       - streamwise gaussian distribution. Wall-refined. eta parallel to y
%%%% fanned         - equidistant streamwise distribution. Wall-normal eta axes to 
%%%%                  account for wall curvature. wall-refined. 
%%%% wallorthogonal - Locally wall-orthogonal grid. eta is curved.
% Grid.mug      (1)         [-]   Streamwise grid refinement location, domain [S S+L]
% Grid.sig      (1)         [-]   Streamwise grid refinement variance (Gaussian), domain [0 1] 
% Grid.ag       (1)         [-]   Streamwise grid refinement strength, domain [0 1]
% Grid.StepX    (1)         [-]   Step location
% Grid.StepH    (1)         [-]   Step Height
% Grid.ystretch (1)         [-]   Wall-normal distribution stretching factor
% Grid.StepType (string)    [-]   Sharp geometry type FFS, BFS (not implemented), GAP(not implemented), SHUMP(not implemented) (default = flat)

% Grid information
Grid.nx = 1500; % [-] # of streamwise (xi) stations
Grid.ny = 50;   % [-] # of wall-normal (eta) collocation points

% Bottom wall coordinates
xw = linspace(BF.X(1),BF.X(end),5000);  % [-] wall x-locations, defined here at 5000 locations, interpolated in DeHNSSo                
yw = 0*xw;                              % [-] Wall is flat
Grid.wall = [xw;yw];                    % [-] Wall description;

% Domain description
Grid.H = max(BF.Y);         % [-] Domain height
Grid.y_i = Grid.H/20;       % [-] Median collocation point height
Grid.S = BF.X(1);           % [-] Start of the domain in wall-coordinate
Grid.L = BF.X(end)-BF.X(1); % [-] Length of the domain in wall-coordinate
Grid.type = "equidistant";  % Select grid type (see above)

%% Stab: Perturbation specifications and boundary conditions 
% Name          size                    units  explanation
% Stab.N        (1)                     [-]    Spectral truncation of beta modes
% Stab.M        (1)                     [-]    spectral truncation of omega modes
% Stab.A0       ((2N+1)x(2M+1),1)       [-]    Initial amplitudes of all modes
% Stab.omega_0  (1)                     [-]    Fundamental angular frequency
% Stab.beta_0   (1)                     [-]    Fundamental spanwise wavenumber
% Stab.IC       (string)                [-]    Initialization method ILST, WALL, LOAD
% Stab.bcw      (any, 3x(2N+1)x(2M+1))  [-]    Inhomogeneous boundary conditions wall per mode (x,u,v,w)^T
% Stab.bcf      (any, 3x(2N+1)x(2M+1))  [-]    Inhomogeneous boundary conditions top per mode (x,u,v,w)^T
% Stab.y0       (1,ny0)                 [-]    Wall-normal distribution of inflow perturbation data
% Stab.u0       (3x(2N+1)x(2M+1)),ny0)  [-]    Normalized streamwise perturbation velocity at inflow 
% Stab.v0       (3x(2N+1)x(2M+1)),ny0)  [-]    Normalized wall-normal perturbation velocity at inflow 
% Stab.w0       (3x(2N+1)x(2M+1)),ny0)  [-]    Normalized spanwise perturbation velocity at inflow 
% Stab.p0       (3x(2N+1)x(2M+1)),ny0)  [-]    Normalized perturbation pressures at inflow 

% Spectral Truncation
Stab.N = 0; % [-]   # of beta modes
Stab.M = 5; % [-]   # of omega modes

% Fundamental frequency
f = 90.58454735;                         % [Hz]     Fundamental frequency
F = (2*pi*f*BF.nu)/(BF.Uref^2)*1e6;      % [-]      Reduced frequency (not used)
Stab.omega_0 = 2*pi*f*(BF.lref/BF.Uref); % [-]      Omega 

% Fundametal spanwise wavenumber
lambda=inf;                        % [m] Spanwise wavelength of primary mode            
Stab.beta_0 = 2*pi*BF.lref/lambda; % [-] Spanwise wavenumber

% Mode initialization
Stab.IC = "ILST"; % Method for Primary mode introduction @ inflow  
Stab.A0=zeros(1,(2*Stab.N+1)*(2*Stab.M+1)); % [-] initialization amplitude vector
Stab.A0(ModeToModeNumber( 1,0,Stab.M,Stab.N))=0.00125*sqrt(2); % [-] Half the desired inflow amplitude
Stab.A0(ModeToModeNumber(-1,0,Stab.M,Stab.N))=0.00125*sqrt(2); % [-] Add a complex conjugate for NonLinear simulations

% Explanation for the indexing above:
% Mode counter vector has size 1 x (2*Stab.N+1)(2*Stab.M+1) 
% which includes complex conjugate modes.
% For example, if Stab.M = 1 and Stab.N = 1, 9 modes are accounted for.
% These modes are counted as
%    M  -1  0  1
% N      _______
% -1|    1  2  3
%  0|    4  5  6  
%  1|    7  8  9
% where the contents of the table are the mode numbers.
% For your convenience, use ModeToModeNumber(m,n,Stab.M,Stab.n) to 
% translate mode notation to a mode counter (e.g. mode (1,0) = mode 6).

%% Opt

% Opt.xb       (1) Buffer starting location
% Opt.kappa    (1) Buffer strength (default = 6)
% Opt.nltbufxb (1) Nonlinear term buffer starting location (=xb by default)
% Opt.Th       (1) Nonlinear introduction threshold
% Opt.Conv     (1) Convergence criterion
% Opt.ConvF    (1) Convergence criterion relaxation factor during ramping
% Opt.Sweep    (1) Output intermediate results flag (true=1, false=0)
% Opt.AFg      (1) Amplitude factor growth rate (default = 1.1)

Opt.xb = 85; % [-] Buffer start as a % of numerical domain

%% Run CHNS
[StabRes,StabGrid,BF] = DeHNSSo(BF,Grid,Stab,Opt);

%% Plotting
% Reset figures
close all 

% Define markers and marker locations
marker = 'os*d>';
xmark = linspace(StabGrid.xun(1),StabGrid.xun(end),20);

% Amplitudes
% Open figure
figure(1)

% Plot HNS results
semilogy(StabGrid.xun,StabRes.A(1,:)/sqrt(2),'k--','linewidth',1.5)
hold on
for j = 2:length(StabRes.omegavec)
    % plot amplitudes
semilogy(StabGrid.xun,StabRes.A(j,:)/sqrt(2),'k-','linewidth',1.5)
    % Interpolate marker locations
ymark = interp1(StabGrid.xun,StabRes.A(j,:)/sqrt(2),xmark);
    % plot markers
semilogy(xmark,ymark,marker(j-1),'color','k','markersize',5)
end

% Plot (N)PSE results
if Stab.M>1 %NPSE
    load('StabRes_Blasius_NPSE.mat')
    semilogy(StabGridNPSE.xun,StabResNPSE.A(1,:)/sqrt(2),'r--','linewidth',1.5)
    hold on
    for j = 2:length(StabResNPSE.omegavec)
        % plot amplitudes
    semilogy(StabGridNPSE.xun,StabResNPSE.A(j,:)/sqrt(2),'r-','linewidth',1.5)
        % Interpolate marker locations
    ymark = interp1(StabGridNPSE.xun,StabResNPSE.A(j,:)/sqrt(2),xmark);
        % plot markers
    semilogy(xmark,ymark,marker(j-1),'color','r','markersize',5)
    end
    ylim([1e-11 1e-1])

else %LPSE
    load('StabRes_Blasius_LPSE.mat')
    semilogy(StabGridLPSE.xun,StabResLPSE.A(1,:)/sqrt(2),'r--','linewidth',1.5)
    hold on
    for j = 2:length(StabResLPSE.omegavec)
        % plot amplitudes
    semilogy(StabGridLPSE.xun,StabResLPSE.A(j,:)/sqrt(2),'r-','linewidth',1.5)
        % Interpolate marker locations
    ymark = interp1(StabGridLPSE.xun,StabResLPSE.A(j,:)/sqrt(2),xmark);
        % plot markers
    semilogy(xmark,ymark,marker(j-1),'color','r','markersize',5)
    end
    ylim([1e-3 5e-2])
end
xlim([400 2200])

% Set label and figure title text
hYLabel = ylabel('$u^{\prime}_{rms}$', 'interpreter', 'latex','rotation',0);
hXLabel = xlabel('$x$', 'interpreter', 'latex');
hTitle = title (['\color[rgb]{0.3 0.3 0.3} Amplitude Development F = ', num2str(StabRes.omegavec(2)/BF.Re*1e6)]);

% Fonts and font sizes
set( gca,'FontName','Times' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName', 'Times');
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 10,  'Rotation',0, 'VerticalAlignment', 'cap');
set( hTitle                    , ...
    'FontSize'   , 10          , ...
    'FontWeight' , 'bold'      );

%Turn on grid
grid on
hold off
clear marker xmark ymark

