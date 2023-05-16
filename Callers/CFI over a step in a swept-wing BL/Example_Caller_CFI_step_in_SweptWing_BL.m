clear all
close all

%Load paths
addpath(genpath('../../..'))

%% General inputs

% Domain definitions
S = 0.1464; % inflow start
E = 0.2125;  % domain end
H = 0.02;   % domain height

% Chebyshev node coordinate median 
y_i = H/20; 

% Kinematic viscocity
nu = 1.4711e-5;

% Freestream velocity [m/s] at inlet
V  = 15.1;

% External velocity
X  = linspace(S,E,5000);
Ue = V*(0.0023*log(X).^4+0.0377*log(X).^3+0.1752*log(X).^2+0.5303*log(X)+1.8574);

% Spanwise velocity
We = -18.7379;

%% NonDim: NonDimensionalization reference values and Re
% Name         size explanation
% NonDim.Re     (1) Reynolds number (=U_ref * l_ref / nu)
% Nondim.Uref   (1) Reference velocity
% NonDim.lref   (1) Reference length
% Nondim.nu     (1) Kinematic viscosity

NonDim.nu   = nu;
NonDim.Uref = 15.1; % reference velocity
NonDim.lref = sqrt(0.0468*nu/NonDim.Uref);% reference length scale (blasius length)
NonDim.Re   = NonDim.Uref*NonDim.lref/nu;

%% BF: Base Flow data and grid

% Name         size explanation
% BF.X   (nxbl,nybl) Base Flow grid streamwise locations 
% BF.Y   (nxbl,nybl) Base Flow grid streamwise locations
% BF.U   (nxbl,nybl) Base Flow Streamwise Velocity
% BF.V   (nxbl,nybl) Base Flow Wall-normal Velocity
% BF.W   (nxbl,nybl) Base Flow Spanwise Velocity
% BF.dxU (nxbl,nybl) Streamwise Gradient of Base Flow Streamwise Velocity
% BF.dxV (nxbl,nybl) Streamwise Gradient of Base Flow Wall-normal Velocity
% BF.dxW (nxbl,nybl) Streamwise Gradient of Base Flow Spanwise Velocity
% BF.dyU (nxbl,nybl) Wall-normal Gradient of Base Flow Streamwise Velocity
% BF.dyV (nxbl,nybl) Wall-normal Gradient of Base Flow Wall-normal Velocity
% BF.dyW (nxbl,nybl) Wall-normal Gradient of Base Flow Spanwise Velocity

% Note: All values in BF should be nondimensionalized by the respective
% reference value matching NonDim.

% Load Base Flow data
 load('BF.mat')



%% Grid: Numerical domain specifications OR Numerical domain grid points
% Name         size explanation
% Grid.nx       (1) Number of streamwise stations of the stability grid
% Grid.ny       (1) Number of wall-normal collocation points of the stability grid
% Grid.wall     (nxwall,2) Wall definition x and y locations on arbitrary grid
% Grid.H        (1) Domain height
% Grid.y_i      (1) Median collocation point height
% Grid.S        (1) Stability grid domain start
% Grid.L        (1) Stability grid domain length (wall length)
% Grid.mode     (string) Grid generation mode (see grid_gen)
    % equidistant    - an equidistant streamwise distribution. Wall-refined. 
    % xrefined       - streamwise gaussian distribution. Wall-refined. eta parallel to y
    % fanned         - equidistant streamwise distribution. Wall-normal eta axes to 
    %                  account for wall curvature. wall-refined. 
    % wallorthogonal - Locally wall-orthogonal grid. eta is curved.
% Grid.mug      (1) Streamwise grid refinement location [S S+L]
% Grid.sig      (1) Streamwise grid refinement variance [0 1] (Gaussian)
% Grid.ag       (1) Streamwise grid refinement strength [0 1]
% Grid.StepX    (1) Step location
% Grid.StepH    (1) Step Height
% Grid.ystretch (1) wall-normal distribution stretching factor
% Grid.StepType (string) Sharp geometry type FFS, BFS (not implemented), GAP(not implemented), SHUMP(not implemented) (default = flat)

% Note: Only BFS StepType currently supported
% Note: Current implementation features EBM. Supply flat wall data and
% seperate step specifications

Grid.nx = 600; % # of streamwise, xi, stability grid stations
Grid.ny = 100;   % # of wall-normal, eta, stability grid collocation points

% Note: This case requires a high refinement to converge. 
% Run with nx = 4000, ny = 350 for a converged result with 1.5Tb of free RAM
    % Should take ~8 hours to solve on HPC

% On PC:
% Run with nx = 1500, ny = 150 for a decent result with ~19Gb of free RAM
    % Should take ~8 minutes to set up LHS and ~30 minutes to solve.

% Run with nx = 1200, ny = 120 for a test and ~10Gb of free RAM
    % should take ~3 minutes to set up LHS and ~7 minutes to solve.

% Run with nx = 600, ny = 100 for a quick test
    % should take ~1 minute to set up LHS and ~1 minute to solve.

xw = linspace(S,E,5000)/NonDim.lref; % wall x-locations
yw = 0*xw;           % Wall is flat
Grid.wall = [xw;yw]; % Wall description [-]

Grid.H   = H/NonDim.lref; % Domain height [-] 

Grid.S = S/NonDim.lref;     % Start of the domain in wall-coordinate [-]
Grid.L = (E-S)/NonDim.lref; % Length of the domain in wall-coordinate [-]
Grid.xtype = "xrefined";  % Select streamwise distribution
Grid.ytype = "step"; % Select wall-normal distribution

Grid.StepX = 0.1837/NonDim.lref; %[m]
Grid.StepH = 7.4710e-04/NonDim.lref; %[m]

Grid.mug = Grid.StepX; %[m]
Grid.sig = 0.1; %[-]
Grid.ag = 0.33; %[-]

Grid.y_i = 10*0.4823;         % Median collocation point height [-]
Grid.ystretch = 1000; 
Grid.StepType = 'FFS';
%% Stab
% Add How to nondimensionalize and that it IS
% Check if bct works inhomogeneously

% Name         size explanation
% Stab.N        (1) Spectral truncation of beta modes
% Stab.M        (1) spectral truncation of omega modes
% Stab.A0       ((2N+1)x(2M+1),1) Initial amplitudes of all modes
% Stab.omega_0  (1) Fundamental angular frequency
% Stab.beta_0   (1) Fundamental spanwise wavenumber
% Stab.IC       (string) Inflow condition ILST, ZERO, LOAD
% Stab.bcw      (any, 3x(2N+1)x(2M+1)) Inhomogeneous boundary conditions wall per mode (x,u,v,w)^T
% Stab.bct      (any, 3x(2N+1)x(2M+1)) Inhomogeneous boundary conditions top per mode (x,u,v,w)^T
% Stab.bcx      (any, 1) Streamwise locations of the boundary condition definitions
%   if Stab.IC == "LOAD", supply:
% Stab.y0       (1,ny0) Wall-normal distribution of inflow perturbation data
% Stab.u0       (3x(2N+1)x(2M+1)),ny0) Normalized streamwise perturbation velocity at inflow 
% Stab.v0       (3x(2N+1)x(2M+1)),ny0) Normalized wall-normal perturbation velocity at inflow 
% Stab.w0       (3x(2N+1)x(2M+1)),ny0) Normalized spanwise perturbation velocity at inflow 
% Stab.p0       (3x(2N+1)x(2M+1)),ny0) Normalized perturbation pressures at inflow 

% Modes
Stab.M = 0; %# of omega modes
Stab.N = 1; %# of beta modes

f = 0; %hz
Stab.omega_0=f*NonDim.lref/NonDim.Uref; % Primary mode omega

lambda = 7.5e-3; % spanwise wavelength (m) of primary mode
Stab.beta_0 = 2*pi*NonDim.lref/lambda; % rad dimensionless

Stab.A0 = zeros(1,(2*Stab.N+1)*(2*Stab.M+1)); % initialization amplitude vector
Stab.A0(ceil((2*Stab.N+1)*(2*Stab.M+1)/2)+1) = 0*2.5e-3/2; % Linear
Stab.A0(ceil((2*Stab.N+1)*(2*Stab.M+1)/2)-1) = 0*2.5e-3/2; % Add if NonLinear

Stab.IC = "ZERO"; % Method for Primary mode introduction @ inflow

    Stab.bcw=zeros(3,5000,ModeToModeNumber(m,n,Stab.M,Stab.N));
    Stab.bct=zeros(3,5000,ModeToModeNumber(m,n,Stab.M,Stab.N));
    Stab.bcx(1,:) = xw;
    xbcwstart = 709+1;
    xbcwend = 738+1;
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
%% Opt

% Opt.xb       (1) Outflow buffer starting location
% Opt.kappa    (1) Outflow buffer strength
% Opt.nltbufxb (1) Nonlinear term outflow buffer start (=xb by default)
% Opt.Th       (1) Nonlinear introduction threshold
% Opt.Sweep    (1) Output intermediate results flag (true=1, false=0)
% Opt.AFg      (1) Amplitude factor growth rate (default = 1.1)
% Opt.Conv     (1) Convergence criterion (default = 1e-4)
% Opt.AMAX     (1) Amplitude limit for linear development (default = 0.1)

Opt.xb = 85; % Buffer start in % of numerical domain
Opt.nltbufxb = 80; % NLT buffer start in % of numerical domain



%% Run CHNS
[StabRes,StabGrid,BF] = DeHNSSo(BF,Grid,Stab,NonDim,Opt);

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
load('StabResDNS.mat')
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
