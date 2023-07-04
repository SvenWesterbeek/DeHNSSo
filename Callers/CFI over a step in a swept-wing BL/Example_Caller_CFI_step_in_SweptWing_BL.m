%% Description
% This file is the caller for the example case that considers the evolution
% of crossflow instabilities over a step in a swept-wing boundary layer using

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
 load('BF_SweptWing_Step.mat')



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

% Note: Only FFS StepType currently supported
% Note: Current implementation features EBM. Supply flat wall data and
% seperate step specifications

Grid.nx = 1500; % # of streamwise, xi, stability grid stations
Grid.ny = 150;   % # of wall-normal, eta, stability grid collocation points

% Note: This case requires a high refinement to converge. 
% Run with nx = 4000, ny = 350 for a converged result with 1.5Tb of free RAM
    % Should take ~8 hours to solve on HPC

% On PC:
% Run with nx = 1500, ny = 150 for a decent result with ~30Gb of free RAM
    % Should take ~8 minutes to set up LHS and ~30 minutes to solve.

% Run with nx = 1200, ny = 120 for a test and ~10Gb of free RAM
    % should take ~3 minutes to set up LHS and ~7 minutes to solve.

% Run with nx = 600, ny = 100 for a quick test
    % should take ~1 minute to set up LHS and ~1 minute to solve.

xw = linspace(0.1464,0.2125,5000)/BF.lref; % wall x-locations
yw = 0*xw;           % Wall is flat
Grid.wall = [xw;yw]; % Wall description [-]

Grid.H   = 0.02/BF.lref; % Domain height [-] 

Grid.S = 0.1464/BF.lref;     % Start of the domain in wall-coordinate [-]
Grid.L = 0.0661/BF.lref; % Length of the domain in wall-coordinate [-]
Grid.xtype = "xrefined";  % Select streamwise distribution
Grid.ytype = "step"; % Select wall-normal distribution

Grid.StepX = 0.1837/BF.lref; %[m]
Grid.StepH = 7.4710e-04/BF.lref; %[m]

Grid.mug = Grid.StepX; %[m]
Grid.sig = 0.1; %[-]
Grid.ag = 0.33; %[-]

Grid.y_i = 4.823;         % Median collocation point height [-]
Grid.ystretch = 1000; 
Grid.StepType = 'FFS';

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
Stab.N = 1; % [-] # of beta modes
Stab.M = 0; % [-] # of omega modes

% Fundamental frequency
f = 0;                                   % [Hz]     Fundamental frequency
F = (2*pi*f*BF.nu)/(BF.Uref^2)*1e6;      % [-]      Reduced frequency (not used)
Stab.omega_0 = 2*pi*f*(BF.lref/BF.Uref); % [-]      Omega 

% Fundametal spanwise wavenumber
lambda=7.5e-3;                     % [m] Spanwise wavelength of primary mode            
Stab.beta_0 = 2*pi*BF.lref/lambda; % [-] Spanwise wavenumber

% Mode initialization
Stab.IC = "ZERO"; % Method for Primary mode introduction @ inflow

    Stab.bcw=zeros(3,5000,ModeToModeNumber(Stab.M,Stab.N,Stab.M,Stab.N));
    Stab.bct=zeros(3,5000,ModeToModeNumber(Stab.M,Stab.N,Stab.M,Stab.N));
    Stab.bcx(1,:) = xw;
    xbcwstart = 709;
    xbcwend = 738;
    for j = 1:(2*Stab.N+1)*(2*Stab.M+1)
        if j == ModeToModeNumber_v2(1,0,Stab.N,Stab.M) % only blowing and suction at (1,0) mode
            for i = 1:length(xw)
                if xw(i)>= xbcwstart && xw(i)<=xbcwend
                Stab.bcw(2,i,j) =  1.414*0.5*1e-5*(4*((xw(i)-xbcwstart)*(xbcwend-xw(i)))/(xbcwend-xbcwstart)^2)^3;
                else
                Stab.bcw(2,i,j)  = 0;
                end
            end
        end
    end

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

Opt.xb = 85; % Buffer start in % of numerical domain
Opt.nltbufxb = 80; % NLT buffer start in % of numerical domain



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

% Plot DNS results
load('StabRes_SweptWing_step_DNS.mat')
semilogy(StabResDNS.X,StabResDNS.A,'m-')
xlim([440 916])
ylim([1e-12 10])

% Set label and figure title text
hYLabel = ylabel('$u^{\prime}_{rms}$', 'interpreter', 'latex','rotation',0);
hXLabel = xlabel('$x$', 'interpreter', 'latex');
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

% Shape Function
