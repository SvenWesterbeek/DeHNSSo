%% Description
% This file is the caller for the example case that considers the evolution
% of crossflow instabilities in a swept-wing boundary layer using

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

%% Clean up and load paths
clear all
close all

%Load paths
addpath(genpath('../../..'))

%% BF: Base Flow data, base flow grid and reference values

% Name         size  unit    explanation
% BF.X   (nxbl,1)    [-]     Base flow grid streamwise locations 
% BF.Y   (1,nybl)    [-]     Base flow grid wall-normal locations
% BF.U   (nxbl,nybl) [-]     Base flow streamwise velocity
% BF.V   (nxbl,nybl) [-]     Base flow wall-normal velocity
% BF.W   (nxbl,nybl) [-]     Base flow spanwise velocity
% BF.dxU (nxbl,nybl) [-]     Streamwise gradient of base flow streamwise velocity
% BF.dxV (nxbl,nybl) [-]     Streamwise gradient of base flow wall-normal velocity
% BF.dxW (nxbl,nybl) [-]     Streamwise gradient of base flow spanwise velocity
% BF.dyU (nxbl,nybl) [-]     Wall-normal gradient of base flow streamwise velocity
% BF.dyV (nxbl,nybl) [-]     Wall-normal gradient of base flow wall-normal velocity
% BF.dyW (nxbl,nybl) [-]     Wall-normal gradient of base flow spanwise velocity

% BF.lref   (1)      [m]     Reference length (Blasius lengthscale)
% BF.Uref   (1)      [m/s]   Reference velocity 
% BF.nu     (1)      [m^2/s] Kinematic viscocity 
% BF.Re     (1)      [-]     Reynolds number

% Load Base Flow data
load('BF_SweptWing_flat.mat')

%% Grid: Numerical domain specifications OR Numerical domain grid points  

% Name         size         units explanation
% Grid.nx       (1)         [-]   Number of streamwise stations of the numerical grid
% Grid.ny       (1)         [-]   Number of wall-normal collocation points of the numerical grid
% Grid.wall     (nxwall,2)  [-]   Wall definition x and y locations 
% Grid.H        (1)         [-]   Domain height
% Grid.y_i      (1)         [-]   Median collocation point height
% Grid.S        (1)         [-]   Numerical grid domain start
% Grid.L        (1)         [-]   Numerical grid domain length (wall length)
% Grid.mode     (string)    [-]   Grid generation mode (see grid_gen)
%   * equidistant    - an equidistant streamwise distribution. Wall-refined. 
%   * xrefined       - streamwise gaussian distribution. Wall-refined. eta parallel to y
%   * fanned         - equidistant streamwise distribution. Wall-normal eta axes to 
%                      account for wall curvature. Wall-refined. 
%   * wallorthogonal - Locally wall-orthogonal grid. Eta is curved.
% Grid.mug      (1)         [-]   Streamwise grid refinement center, domain [S S+L]
% Grid.sig      (1)         [-]   Streamwise grid refinement variance (Gaussian), domain [0 1] 
% Grid.ag       (1)         [-]   Streamwise grid refinement strength, domain [0 1]
% Grid.StepX    (1)         [-]   Step location
% Grid.StepH    (1)         [-]   Step Height
% Grid.ystretch (1)         [-]   Wall-normal distribution stretching factor
% Grid.StepType (string)    [-]   Sharp geometry type "flat", "FFS" (default = flat)

% Define grid size
Grid.nx = 1200; % [-] # of streamwise (xi) stations
Grid.ny = 40;   % [-] # of wall-normal (eta) collocation points
 
% Define bottom wall
xw = linspace(BF.X(1),BF.X(end),5000); % [-] Wall x-locations
yw = 0*xw;                             % [-] Wall is flat
Grid.wall = [xw;yw];                   % [-] Wall description [-]

% Set domain
Grid.H = max(BF.Y)*0.99;        % [-] Domain height
Grid.y_i = Grid.H/10;           % [-] Median collocation point height
Grid.S = BF.X(1);               % [-] Start of the domain in wall-coordinate
Grid.L = BF.X(end-1)-BF.X(1);   % [-] Length of the domain in wall-coordinate

Grid.type = "equidistant";  % Select the type of grid

%% Stab: Perturbation specifications and boundary conditions 

% Name          size                    units  explanation
% Stab.N        (1)                     [-]    Spectral truncation of beta modes
% Stab.M        (1)                     [-]    Spectral truncation of omega modes
% Stab.A0       ((2N+1)x(2M+1),1)       [-]    Initial amplitudes of all modes
% Stab.omega_0  (1)                     [-]    Fundamental angular frequency
% Stab.beta_0   (1)                     [-]    Fundamental spanwise wavenumber
% Stab.IC       (string)                [-]    Initialization method ILST, WALL, LOAD
% Stab.bcw      (any, 3x(2N+1)x(2M+1))  [-]    Inhomogeneous boundary conditions wall per mode (x,u,v,w)^T
% Stab.bcf      (any, 3x(2N+1)x(2M+1))  [-]    Inhomogeneous boundary conditions top per mode (x,u,v,w)^T
% Stab.y0       (1,ny0)                 [-]    Wall-normal distribution of inflow perturbation data
% Stab.u0       (3x(2N+1)x(2M+1)),ny0)  [-]    Streamwise perturbation velocity shape function at inflow 
% Stab.v0       (3x(2N+1)x(2M+1)),ny0)  [-]    Wall-normal perturbation velocity shape function at inflow 
% Stab.w0       (3x(2N+1)x(2M+1)),ny0)  [-]    Spanwise perturbation velocity shape function at inflow 
% Stab.p0       (3x(2N+1)x(2M+1)),ny0)  [-]    Perturbation pressures shape function at inflow 

% Spectral Truncation
Stab.N = 5; % [-] # of beta modes
Stab.M = 0; % [-] # of omega modes

% Fundamental frequency
f = 0;                                   % [Hz]     Fundamental frequency
F = (2*pi*f*BF.nu)/(BF.Uref^2)*1e6;      % [-]      Reduced frequency (not used)
Stab.omega_0 = 2*pi*f*(BF.lref/BF.Uref); % [-]      Omega 

% Fundametal spanwise wavenumber
lambda=7.5e-3;                     % [m] Spanwise wavelength of primary mode            
Stab.beta_0 = 2*pi*BF.lref/lambda; % [-] Spanwise wavenumber

% Mode initialization
Stab.IC = "ILST"; % Method for Primary mode introduction @ inflow  
Stab.A0=zeros(1,(2*Stab.N+1)*(2*Stab.M+1));               % [-] initialization amplitude vector
Stab.A0(ModeToModeNumber(0, 1,Stab.M,Stab.N))=3.5e-3/2;   % [-] Half the desired inflow amplitude
Stab.A0(ModeToModeNumber(0,-1,Stab.M,Stab.N))=3.5e-3/2; % [-] Add a complex conjugate for NonLinear simulations

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

% Name         size units explanation
% Opt.xb       (1)  [-]   Buffer starting location as a % of the domain (default = 85)
% Opt.kappa    (1)  [-]   Buffer strength (default = 6)
% Opt.nltbufxb (1)  [-]   Nonlinear term buffer starting location (=xb by default)
% Opt.Th       (1)  [-]   Nonlinear introduction threshold (default = 1e-11)
% Opt.Conv     (1)  [-]   Convergence criterion (default = 1e-4)
% Opt.ConvF    (1)  [-]   Convergence criterion relaxation factor during ramping (default = 100)
% Opt.Sweep    (1)  [-]   Output intermediate results flag (true=1,false=0) (default = 0)
% Opt.AFg      (1)  [-]   Amplitude factor growth rate (default = 1.1)

Opt.xb = 85;       % [-] Buffer start as a % of numerical domain
Opt.nltbufxb = 80; % [-] NLT buffer start as a % of numerical domain


%% Run CHNS
tstart = tic;
[StabRes,StabGrid,BF] = DeHNSSo(BF,Grid,Stab,Opt);
time = toc(tstart);

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
if any(StabRes.A(1,:)) %Plot MFD
semilogy(StabGrid.xun,StabRes.A(1,:)/sqrt(2),'k--','linewidth',1.5)
hold on
end
for j = 2:length(StabRes.betavec)
    % plot amplitudes
semilogy(StabGrid.xun,StabRes.A(j,:)/sqrt(2),'k-','linewidth',1.5, 'HandleVisibility', 'off')
hold on
    % Interpolate marker locations
ymark = interp1(StabGrid.xun,StabRes.A(j,:)/sqrt(2),xmark);
    % plot markers
semilogy(xmark,ymark,marker(j-1),'color','k','markersize',5)
end

% Plot Reference results
if Stab.N > 1
    load('StabRes_SweptWing_flat_NPSE.mat')
if any(StabRes2.A(1,:)) %Plot MFD
semilogy(StabGrid2.x,StabRes2.A(1,:)/sqrt(2),'r--','linewidth',1.5)
hold on
end
for j = 2:size(StabRes2.A,1)
    % plot amplitudes
semilogy(StabGrid2.x,StabRes2.A(j,:)/sqrt(2),'r-','linewidth',1.5, 'HandleVisibility', 'off')
hold on
    % Interpolate marker locations
ymark = interp1(StabGrid2.x,StabRes2.A(j,:)/sqrt(2),xmark);
    % plot markers
semilogy(xmark,ymark,marker(j-1),'color','r','markersize',5)
end

load('StabRes_SweptWing_flat_DNS.mat')
semilogy(StabGridDNS.xun, StabResDNS.A(:,1)/sqrt(2),'m--')
for j = 2:size(StabRes.A,1)
    % plot amplitudes
semilogy(StabGridDNS.xun, StabResDNS.A(:,j)/sqrt(2),'m-', 'HandleVisibility', 'off')
    % Interpolate marker locations
ymark = interp1(StabGridDNS.xun,StabResDNS.A(:,j)/sqrt(2),xmark);
    % plot markers
semilogy(xmark,ymark,marker(j-1),'color','m','markersize',5)
end
ylim([1e-12 10])

else %Stab.N = 1
load('StabRes_SweptWing_flat_LPSE.mat')
semilogy(StabGridLPSE.xun,StabResLPSE.A(1,:)/sqrt(2),'r--','linewidth',1.5)
ylim([1e-3 10])
end

xlim([218.9 1400])

% Set label and figure title text
hYLabel = ylabel('$u^{\prime}_{rms}$', 'interpreter', 'latex','rotation',0);
hXLabel = xlabel('$x$', 'interpreter', 'latex');
hTitle = title ('Amplitude development CFI in a swept-wing BL','FontName','Times New Roman','FontSize',10);

% Fonts and font sizes
set( gca,'FontName','Times' );
set([ hXLabel, hYLabel], ...
    'FontName', 'Times');
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 10,  'Rotation',0, 'VerticalAlignment', 'cap');              

%Turn on grid
grid on
hold off
clear marker xmark ymark

% HNS modes
hnsmodeslegend = [];
for jj = 1:length(StabRes.betavec)
    if any(StabRes.A(jj,:))
% hnsmodeslegend = [hnsmodeslegend, convertCharsToStrings(['HNS mode (m,n) = (' num2str(StabRes.betavec(jj)/StabRes.betavec(2)) ',' num2str(0) ')']),convertCharsToStrings(['LPSE mode (m,n) = (' num2str(StabRes.betavec(jj)/StabRes.betavec(2)) ',' num2str(0) ') marker'])  ];
hnsmodeslegend = [hnsmodeslegend, convertCharsToStrings(['LPSE mode (m,n) = (' num2str(StabRes.betavec(jj)/StabRes.betavec(2)) ',' num2str(0) ')'])  ];
    end
end


if Stab.N>1 %NPSE
% PSE modes
psemodeslegend = [];
for jj = 1:length(StabRes.betavec)
    if any(StabRes2.A(jj,:))
psemodeslegend = [psemodeslegend, convertCharsToStrings(['NPSE mode (m,n) = (' num2str(StabRes.betavec(jj)/StabRes.betavec(2)) ',' num2str(0) ')'])  ];
    end
end

% DNS modes
dnsmodeslegend = [];
for jj = 1:length(StabRes.betavec)
    if any(StabRes.A(jj,:))
dnsmodeslegend = [dnsmodeslegend, convertCharsToStrings(['HNS mode (m,n) = (' num2str(StabRes.betavec(jj)/StabRes.betavec(2)) ',' num2str(0) ')'])];
    end
end
else
% PSE modes
psemodeslegend = [convertCharsToStrings(['LPSE mode (m,n) = (' num2str(StabRes.betavec(2)/StabRes.betavec(2)) ',' num2str(0) ')'])];
dnsmodeslegend = [];
end


legend([hnsmodeslegend psemodeslegend dnsmodeslegend],'Location','EastOutside')


