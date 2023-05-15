function [StabRes,StabGrid,BF] = DeHNSSo(BF,Grid,Stab,Opt)

% ██████████            █████   █████ ██████   █████  █████████   █████████          
%░░███░░░░███          ░░███   ░░███ ░░██████ ░░███  ███░░░░░███ ███░░░░░███         
% ░███   ░░███  ██████  ░███    ░███  ░███░███ ░███ ░███    ░░░ ░███    ░░░   ██████ 
% ░███    ░███ ███░░███ ░███████████  ░███░░███░███ ░░█████████ ░░█████████  ███░░███
% ░███    ░███░███████  ░███░░░░░███  ░███ ░░██████  ░░░░░░░░███ ░░░░░░░░███░███ ░███
% ░███    ███ ░███░░░   ░███    ░███  ░███  ░░█████  ███    ░███ ███    ░███░███ ░███
% ██████████  ░░██████  █████   █████ █████  ░░█████░░█████████ ░░█████████ ░░██████ 
%░░░░░░░░░░    ░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░  ░░░░░░░░░   ░░░░░░░░░   ░░░░░░  
                                                                                          
% Made using: https://manytools.org/hacker-tools/ascii-banner/                                                   
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

% Note: Once this is complete and checked i will paste it in every subfunction.

%% Description
 
% Executes Nonlinear Harmonic Navier-Stokes simulation through:
% 1. Creating a numerical grid for for a domain defined in "Grid"
% 2. Set up mode interaction matrices for spectral content defined in "Stab"
% 3. Load a base flow presented in "BF"
% 4. Set up outflow buffer and embedded boundaries
% 5. Interpolate base flow onto the numerical grid
% 6. Load inflow perturbation or solve local eigenvalue problem at inflow
% 7. Solve HNS
% 8. Post-process results

% Everything should be presented in nondimensional quantities, with 
% reference values defined in BF. Solver options should be presented in 
% Opt.

% The output contains two structs: StabGrid and StabRes. StabGrid contains
% the global and numerical domains as well as the transformation
% coefficients. StabRes contains the stability results, providing shape
% functions, amplitude development and wavenumbers.

% Authors: Sven Westerbeek and Marios Kotsonis
% last update: March 2023

% Article: 
% DOI: 
% 

% For comments, questions, suggestions, ideas or collaborations, 
% please contact us at:
% S.H.J.Westerbeek@tudelft.nl; svenwesterbeek@gmail.nl
% M.Kotsonis@tudelft.nl



%% Updates
% 04-2023: V1 Published 

%% Overview of inputs

% BF: Base Flow data + grid and Reference values + Reynolds number
% Name    size        unit    explanation
% BF.X    (nxbl,nybl) [-]     Base Flow grid streamwise locations 
% BF.Y    (nxbl,nybl) [-]     Base Flow grid streamwise locations
% BF.U    (nxbl,nybl) [-]     Base Flow Streamwise Velocity
% BF.V    (nxbl,nybl) [-]     Base Flow Wall-normal Velocity
% BF.W    (nxbl,nybl) [-]     Base Flow Spanwise Velocity
% BF.dxU  (nxbl,nybl) [-]     Streamwise Gradient of Base Flow Streamwise Velocity
% BF.dxV  (nxbl,nybl) [-]     Streamwise Gradient of Base Flow Wall-normal Velocity
% BF.dxW  (nxbl,nybl) [-]     Streamwise Gradient of Base Flow Spanwise Velocity
% BF.dyU  (nxbl,nybl) [-]     Wall-normal Gradient of Base Flow Streamwise Velocity
% BF.dyV  (nxbl,nybl) [-]     Wall-normal Gradient of Base Flow Wall-normal Velocity
% BF.dyW  (nxbl,nybl) [-]     Wall-normal Gradient of Base Flow Spanwise Velocity
% BF.Re   (1)         [-]     Reynolds number (=U_ref * l_ref / nu)
% BF.Uref (1)         [m/s]   Reference velocity
% BF.lref (1)         [m]     Reference length
% BF.nu   (1)         [m^2/s] Kinematic viscosity

% Grid: Stability domain and grid
% Name          size     unit explanation
% Grid.nx       (1)      [-]  Number of streamwise stations of the numerical grid
% Grid.ny       (1)      [-]  Number of wall-normal collocation points of the numerical grid
% Grid.wall     (any,2)  [-]  Wall definition x and y locations on arbitrary grid
% Grid.H        (1)      [-]  Domain height
% yi            (1)      [-]  Median collocation point height
% Grid.S        (1)      [-]  Stability grid domain start
% Grid.L        (1)      [-]  Stability grid domain length (wall length)
% Grid.mode     (string) [-]  Grid generation mode (see grid_gen)
% Grid.ft       (1)      [-]  Flat top boundary flag 0=no, 1=yes (default = 1)
% Grid.mug      (1)      [-]  Streamwise grid refinement location [S S+L]
% Grid.sig      (1)      [-]  Streamwise grid refinement variance [0 1] (Gaussian)
% Grid.ag       (1)      [-]  Streamwise grid refinement strength [0 1]
% Grid.StepX    (1)      [-]  Step location
% Grid.StepH    (1)      [-]  Step Height
% Grid.ystretch (1)      [-]  wall-normal distribution stretching factor
% Grid.StepType (string) [-]  Sharp geometry type FFS, BFS, GAP, HUMP

% Stab: Stability specifications
% Name          size                    unit explanation
% Stab.N        (1)                     [-]  Spectral truncation of beta modes
% Stab.M        (1)                     [-]  Spectral truncation of omega modes
% Stab.A0       ((2N+1)x(2M+1),1)       [-]  Initial amplitudes of all modes
% Stab.omega_0  (1)                     [-] Fundamental angular frequency
% Stab.beta_0   (1)                     [-] Fundamental spanwise wavenumber
% Stab.IC       (string)                [-] Initialization method ILST, ZERO, LOAD
% Stab.bcw      (any, 3x(2N+1)x(2M+1))  [-] Inhomogeneous boundary conditions wall per mode (x,u,v,w)^T
% Stab.bct      (any, 3x(2N+1)x(2M+1))  [-] Inhomogeneous boundary conditions top per mode (x,u,v,w)^T
% Stab.y0       (1,ny0)                 [-] Wall-normal distribution of inflow perturbation data
% Stab.u0       (3x(2N+1)x(2M+1)),ny0)  [-] Normalized streamwise perturbation velocity at inflow 
% Stab.v0       (3x(2N+1)x(2M+1)),ny0)  [-] Normalized wall-normal perturbation velocity at inflow 
% Stab.w0       (3x(2N+1)x(2M+1)),ny0)  [-] Normalized spanwise perturbation velocity at inflow 
% Stab.p0       (3x(2N+1)x(2M+1)),ny0)  [-] Normalized perturbation pressures at inflow 

% Opt: Solver options
% Name         size  unit explanation
% Opt.xb       (1)   [-]  Buffer starting location
% Opt.kappa    (1)   [-]  Buffer strength (default = 6)
% Opt.nltbufxb (1)   [-]  Nonlinear term buffer starting location (=xb by default)
% Opt.Th       (1)   [-]  Nonlinear introduction threshold
% Opt.Conv     (1)   [-]  Convergence criterion
% Opt.ConvF    (1)   [-]  Convergence criterion relaxation factor during ramping
% Opt.Sweep    (1)   [-]  Output intermediate results flag (true=1, false=0)
% Opt.AFg      (1)   [-]  Amplitude factor growth rate (default = 1.1)
%

%% Overview of outputs

% StabGrid: Numerical domain and grid
% Name             size    unit explanation
% StabGrid.x       (nx,ny) [-]  x-locations
% StabGrid.y       (nx,ny) [-]  y-locations
% StabGrid.xi      (nx,ny) [-]  xi-locations
% StabGrid.eta     (nx,ny) [-]  eta-locations
% StabGrid.etaun   (nx,ny) [-]  eta-locations of first column
% StabGrid.xiun    (nx,ny) [-]  xi-locations of the last row
% StabGrid.yun     (nx,ny) [-]  y-locations of first column
% StabGrid.xun     (nx,ny) [-]  x-locations of the last row
% StabGrid.J       (nx,ny) [-]  Jacobian
% StabGrid.xix     (nx,ny) [-]  dx/dxi
% StabGrid.etax    (nx,ny) [-]  deta/dx
% StabGrid.xxi     (nx,ny) [-]  dx/dxi
% StabGrid.xeta    (nx,ny) [-]  dx/deta
% StabGrid.yxi     (nx,ny) [-]  dy/dxi
% StabGrid.yeta    (nx,ny) [-]  dy/deta
% StabGrid.xxixi   (nx,ny) [-]  d^2x/dxi^2
% StabGrid.yxixi   (nx,ny) [-]  d^2y/dxi^2
% StabGrid.xetaeta (nx,ny) [-]  d^2x/deta^2
% StabGrid.yetaeta (nx,ny) [-]  d^2y/deta^2
% StabGrid.xixx    (nx,ny) [-]  d^2xi/dx^2
% StabGrid.etaxx   (nx,ny) [-]  d^2eta/dx^2
% StabGrid.xiyy    (nx,ny) [-]  d^2xi/dy^2
% StabGrid.etayy   (nx,ny) [-]  d^2eta/dy^2

% StabRes: Stability results
% Name              size                    unit explanation
% StabRes.omegavec  (1,(2N+1)x(2M+1))       [-]  Radial frequency
% StabRes.betavec   (1,(2N+1)x(2M+1))       [-]  Spanwise wavenumber
% StabRes.phi       ((2N+1)x(2M+1),4ny,nx)  [-]  State vector of variables (u,v,w,p)^T
% StabRes.alpha     ((2N+1)x(2M+1),nx)      [-]  Streamwise wavenumber
% StabRes.A         ((2N+1)x(2M+1),nx)      [-]  Amplitude, max((abs(u'))
% StabRes.u         ((2N+1)x(2M+1),ny,nx)   [-]  Streamwise perturbation velocity
% StabRes.v         ((2N+1)x(2M+1),ny,nx)   [-]  Wall-normal perturbation velocity
% StabRes.w         ((2N+1)x(2M+1),ny,nx)   [-]  Spanwise perturbation velocity
% StabRes.p         ((2N+1)x(2M+1),ny,nx)   [-]  Perturbation Pressure

% BF: Base flow (Additions to input struct only, interpolated on num. grid)
% Name      size    unit explanation
% BF.Ur     (nx,ny) [-]  Streamwise velocity
% BF.Vr     (nx,ny) [-]  Wall-normal velocity
% BF.Wr     (nx,ny) [-]  Spanwise velocity
% BF.dxUr   (nx,ny) [-]  x-derivative of streamwise velocity
% BF.dxVr   (nx,ny) [-]  x-derivative of wall-normal velocity
% BF.dxWr   (nx,ny) [-]  x-derivative of spanwise velocity
% BF.dyUr   (nx,ny) [-]  y-derivative of streamwise velocity
% BF.dyVr   (nx,ny) [-]  y-derivative of wall-normal velocity
% BF.dyWr   (nx,ny) [-]  y-derivative of spanwise velocity


%% Load default settings

% Stab: Homogeneous BC's by default
if ~isfield(Stab,'bcx'); Stab.bcx = linspace(Grid.S,Grid.S+Grid.L,10); end
if ~isfield(Stab,'bcw'); Stab.bcw = zeros(3,10,3*(2*Stab.N+1)*(2*Stab.M+1)+1); end
if ~isfield(Stab,'bct'); Stab.bct = zeros(3,10,3*(2*Stab.N+1)*(2*Stab.M+1)+1); end
if ~isfield(Stab,'A0'); Stab.A0 = []; end

% Grid: Straight geometry without streamwise refinement by default
if ~isfield(Grid,'StepType'), Grid.StepType = "FLAT"; end
if ~isfield(Grid,'ytype'), Grid.ytype ="malik"; end
if ~isfield(Grid,'xtype'), Grid.xtype ="equidistant"; end
if ~isfield(Grid,'StepX'), Grid.StepX = 0; end
if ~isfield(Grid,'StepH'), Grid.StepH = 0; end
if ~isfield(Grid,'mug'), Grid.mug = 0; end
if ~isfield(Grid,'sig'), Grid.sig = 1; end
if ~isfield(Grid,'ag'), Grid.ag = 0; end
if ~isfield(Grid,'ft'), Grid.ft = 1; end

% Opt: Set default solver specifications
if ~isfield(Opt,'xb'), Opt.xb = 85; end
if ~isfield(Opt,'AFg'),  Opt.AFg = 1.1; end
if ~isfield(Opt,'Sweep'), Opt.Sweep = 0; end
if ~isfield(Opt,'nltbufxb'); Opt.nltbufxb = Opt.xb; end
if ~isfield(Opt,'TH'); Opt.TH = 1e-11; end
if ~isfield(Opt,'Conv'); Opt.Conv = 1e-4 ; end
if ~isfield(Opt,'kappa'); Opt.kappa = 6; end
if ~isfield(Opt,'AMAX'); Opt.AMAX = 0.1; end
if ~isfield(Opt,'ConvF'); Opt.ConvF = 100; end

% Define imaginary unit
iu=sqrt(-1);

%% Setup mode vector and harmonic balancing matrices.

% Find all mode interactions
[Nmat, Mmat, Modevec,Mvec,Nvec] = Mint(Stab.M,Stab.N);

% Set up Harmonic Balancing matrix 
[HB] = Hbalancing(Stab.M, Stab.N, Mmat, Nmat, Modevec ); 

% Create vectors of omegas and betas per mode
StabRes.omegavec = Stab.omega_0*Mvec;
StabRes.betavec = Stab.beta_0*Nvec;

% Define mode counter L and count total number of modes nf
L = 1:(2*Stab.N+1)*(2*Stab.M+1); 
nf = max(L); 


%% Grid Generation
fprintf('Generating grid. \n')
[StabGrid,D1,D2]=grid_gen(Grid);

%% Base Flow and BC Interpolation on stability grd

fprintf('Interpolating base flow on the numerical grid. \n')
% Check DNS base flow from step loading

% Interpolate base flow on numerical grid
% fillmissing corrects possible numerical differences in the wall description
BF.Ur = fillmissing(griddata(BF.X,BF.Y,BF.U,StabGrid.x,StabGrid.y,'cubic'),'pchip');
BF.Vr = fillmissing(griddata(BF.X,BF.Y,BF.V,StabGrid.x,StabGrid.y,'cubic'),'pchip');
BF.Wr = fillmissing(griddata(BF.X,BF.Y,BF.W,StabGrid.x,StabGrid.y,'cubic'),'pchip');

%Enforce no-slip condition
BF.Ur(end,:)=0; BF.Vr(end,:)=0; BF.Wr(end,:)=0;

% Interpolate base flow derivatives on numerical grid
% fillmissing corrects possible numerical differences in the wall description
BF.dxUr = fillmissing(griddata(BF.X,BF.Y,BF.dxU,StabGrid.x,StabGrid.y,'cubic'),'pchip');
BF.dxVr = fillmissing(griddata(BF.X,BF.Y,BF.dxV,StabGrid.x,StabGrid.y,'cubic'),'pchip');
BF.dxWr = fillmissing(griddata(BF.X,BF.Y,BF.dxW,StabGrid.x,StabGrid.y,'cubic'),'pchip');

BF.dyUr = fillmissing(griddata(BF.X,BF.Y,BF.dyU,StabGrid.x,StabGrid.y,'cubic'),'pchip');
BF.dyVr = fillmissing(griddata(BF.X,BF.Y,BF.dyV,StabGrid.x,StabGrid.y,'cubic'),'pchip');
BF.dyWr = fillmissing(griddata(BF.X,BF.Y,BF.dyW,StabGrid.x,StabGrid.y,'cubic'),'pchip');

% Define Inhomogeneous BC's on numerical grid
for k = 1:3 % loop over u v w
    for j = 1:nf % loop over all modes
    bcw(k,:,j) = interp1(Stab.bcx,Stab.bcw(k,:,j),StabGrid.xun,"cubic");
    bct(k,:,j) = interp1(Stab.bcx,Stab.bct(k,:,j),StabGrid.xun,"cubic");
    end
end


%% Inflow BC definition

nonzeromodeswall = [];  % initialize nonzeromodeswall
nonzeromodesic = [];    % initialize nonzeromodesic

% Initialize solution vector phi = [u v w p]'(x,y) and alpha = alpha(x)
StabRes.phi    = zeros(nf,4*Grid.ny,Grid.nx);
StabRes.alpha  = zeros(nf,Grid.nx);
StabRes.A      = zeros(nf,Grid.nx);
StabRes.u      = zeros(nf,Grid.ny,Grid.nx);
StabRes.v      = zeros(nf,Grid.ny,Grid.nx);
StabRes.w      = zeros(nf,Grid.ny,Grid.nx);
StabRes.p      = zeros(nf,Grid.ny,Grid.nx);
f              = zeros(nf,4*Grid.ny*Grid.nx);


if Stab.IC =="LOAD"
    [~,nonzeromodesic]= find(Stab.A0>0);  

    % Load normalized shape functions and impose inflow amplitude
    for j = 1:nf
    StabRes.u(j,:,1) = interp1(Stab.y0,Stab.A0(j)*Stab.u0(j,:),StabGrid.yun,'pchip');
    StabRes.v(j,:,1) = interp1(Stab.y0,Stab.A0(j)*Stab.v0(j,:),StabGrid.yun,'pchip');
    StabRes.w(j,:,1) = interp1(Stab.y0,Stab.A0(j)*Stab.w0(j,:),StabGrid.yun,'pchip');
    StabRes.p(j,:,1) = interp1(Stab.y0,Stab.A0(j)*Stab.p0(j,:),StabGrid.yun,'pchip');
    StabRes.phi(j,:,1) = [StabRes.u(j,:,1) StabRes.v(j,:,1) StabRes.w(j,:,1) StabRes.p(j,:,1)];
    end

elseif Stab.IC=="ILST"
    [~,nonzeromodesic]= find(Stab.A0>0);  

    for j = flip(nonzeromodesic) %Run initial condition for all nonzero modes 
    if j >= round(nf/2)% Mode is physical
    [StabRes]=IC_HNS(j,1,BF.Re,BF.Ur,BF.Wr,BF.dyUr,BF.dyWr,D1,D2,StabRes,StabGrid); 
           
    %output from CIC is normalized by max(u'), superimpose Amplitude
    StabRes.u(j,:,1) = Stab.A0(j).*StabRes.u(j,:,1);
    StabRes.v(j,:,1) = Stab.A0(j).*StabRes.v(j,:,1);
    StabRes.w(j,:,1) = Stab.A0(j).*StabRes.w(j,:,1);
    StabRes.p(j,:,1) = Stab.A0(j).*StabRes.p(j,:,1);
    StabRes.phi(j,:,1) = [StabRes.u(j,:,1) StabRes.v(j,:,1) StabRes.w(j,:,1) StabRes.p(j,:,1)];
    
    
    else % Apply symmetry for conjugate mode
    j2 = round(nf/2)-(j-round(nf/2));% Defines conjugate mode
    
    StabRes.u(j,:,1) = conj(StabRes.u(j2,:,1));
    StabRes.v(j,:,1) = conj(StabRes.v(j2,:,1));
    StabRes.w(j,:,1) = conj(StabRes.w(j2,:,1));
    StabRes.p(j,:,1) = conj(StabRes.p(j2,:,1));
    StabRes.alpha(j,1) = -conj(StabRes.alpha(j,1)); % growth rates should be equal between conjugate modes so -conj
    StabRes.phi(j,:,1) = conj(StabRes.phi(j2,:,1));
    end                                                                            
    end 
    
    if max(abs(Stab.bcw(2:end,:)))==0
    fprintf('Initialising with ILST solution at inflow. \n')
    else
    fprintf('Initialising with ILST solution at inflow and non-homogeneous b.c. at wall. \n')
    end

elseif Stab.IC=="ZERO"
    fprintf('Initialising only with non-homogeneous b.c. at wall. \n')

    % Determine which modes will be nonzero due to wall forcing
    for j = 1:nf
        if any(Stab.bcw(:,:,j),'all')
        nonzeromodeswall = [nonzeromodeswall j];
        end
    end
else

fprintf('Error: no valid inflow condition. \n')
return

end

nonzeromodes = unique([nonzeromodesic, nonzeromodeswall]);
RunJ = nonzeromodes; % Used in loops to assess only active modes

% Define phiIC
phiIC = squeeze(StabRes.phi(:,:,1));


%% create outflow buffer region

% Set up outflow buffer
ibuf=round(Opt.xb*Grid.nx/100); % x-station where buffer starts
bufc=ones(1,Grid.nx); % Pre-allocate buffer vector bufc

% Create attenuation function
for i=ibuf:Grid.nx
bufc(i)=0.5*(1+tanh(Opt.kappa*(1-2*(i-ibuf)/(Grid.nx-ibuf))));  % Buffer definition adapted from Joslin NASA TR3205
end

% Scale attenuation function for smooth introduction
bufc(ibuf:Grid.nx)=bufc(ibuf:Grid.nx)./bufc(ibuf); % Stretch results to be smooth at ibuf

% Set up outflow forcing buffer
ibuff=round((Opt.nltbufxb)*Grid.nx/100); % x-station where forcing buffer starts
bufcf=ones(1,Grid.nx); % Pre-allocate buffer vector bufcf

% Create attenuation function
for i=ibuff:Grid.nx
bufcf(i)=0.5*(1+tanh(Opt.kappa*(1-2*(i-ibuff)/(Grid.nx-ibuf))));  % from Joslin NASA TR3205
end

% Scale attenuation function for smooth introduction
bufcf(ibuff:Grid.nx)=bufcf(ibuff:Grid.nx)./bufcf(ibuff);

fprintf(' Buffer region from %.3g to 100 percent of the domain. \n',Opt.xb)

%% Set up step EBM for sharp surfaces 

if strcmpi(Grid.StepType,"FFS")
    % Forward-facing step
    bufs =1-(repmat(StabGrid.xun,Grid.ny,1)>Grid.StepX-0.00001).*(repmat(StabGrid.yun,1,Grid.nx)<=Grid.StepH); %=0 velocity buffer under step and at step surface
    bufsp =1-(repmat(StabGrid.xun,Grid.ny,1)>Grid.StepX+0.00001).*(repmat(StabGrid.yun,1,Grid.nx)<Grid.StepH-0.001); %=0 pressure buffer under step

    [~,stepi] = min(abs(StabGrid.x(end,:)-Grid.StepX)); % streamwise station of the step face
    [~,yi] = min(abs(StabGrid.y(:,stepi)-Grid.StepH)); % collocation point of the step corner
    
    bufsp(end,stepi)= 0; % pressure buffer value =0 at the inner corner exactly

    fprintf(' Buffer region implemented under the surface post-step. \n')
elseif strcmpi(Grid.StepType,"BFS")  % Backward Facing Step
    % Not yet implemented

elseif strcmpi(Grid.StepType,"GAP")  % Gap
    % Not yet implemented

elseif strcmpi(Grid.StepType,"SHUMP")% Sharp hump
    % Not yet implemented

else %No sharp features
    if Grid.StepH==0 % If no step is present, all values should be 1
        bufs = ones(Grid.ny,Grid.nx);
        bufsp =ones(Grid.ny,Grid.nx);
        stepi = 0;
        yi = 1;
    end
end




%% initialisation of the matrices

I      = eye(Grid.ny);
Z      = zeros(Grid.ny);

ML1MFD = spalloc(4*Grid.ny,4*Grid.ny*Grid.nx,4*(3*Grid.ny^2+4*Grid.ny)+7*Grid.ny^2+4*Grid.ny);  % LHS preallocate sparse matrix per i station. !!! dependent on formulation of finite differences
ML1 =    spalloc(4*Grid.ny,4*Grid.ny*Grid.nx,4*(3*Grid.ny^2+4*Grid.ny)+7*Grid.ny^2+4*Grid.ny);  % LHS preallocate sparse matrix per i station. !!! dependent on formulation of finite differences
ML2 = spalloc(4*Grid.ny,4*Grid.ny*Grid.nx,3*Grid.ny);  % LHS preallocate sparse matrix per i station. !!! dependent on formulation of finite differences
ML3 = spalloc(4*Grid.ny,4*Grid.ny*Grid.nx,5*Grid.ny);  % LHS preallocate sparse matrix per i station. !!! dependent on formulation of finite differences
ML4 = spalloc(4*Grid.ny,4*Grid.ny*Grid.nx,3*Grid.ny);  % LHS preallocate sparse matrix per i station. !!! dependent on formulation of finite differences

R=zeros(4*Grid.nx*Grid.ny,1);         % RHS

iML1MFDtotal = []; iML1total = []; iML2total = []; iML3total = []; iML4total = [];
jML1MFDtotal = []; jML1total = []; jML2total = []; jML3total = []; jML4total = [];
sML1MFDtotal = []; sML1total = []; sML2total = []; sML3total = []; sML4total = [];

% construction of the matrix blocks

dxi  =   StabGrid.xiun(1,2) - StabGrid.xiun(1,1);

%% Initialize Nonlinear convergence loop

Aactive = zeros(size(L));       % Active mode registrer
Aactive(nonzeromodes) = 1;      % Set active modes to introduced modes
Aold = zeros(nf,Grid.nx);       % Amplitudes of previous iteration
q = zeros(nf,4*Grid.ny*Grid.nx);% Solution vector
TH = Opt.TH;                    % Nonlinear introduction threshold


%% Nonlinear convergence loop

 dal = ones(nf,1);   % Convergence measure
 dalsave(:,1) = dal; % Convergence measure register
 Conv = Opt.Conv;    % Convergence threshold
 mmry = 0;           % memory variable so that LHS is created only once
 k1 = 0;             % Iteration counter
 k2 = 0;             % Converged iteration counter
 AF = 1;             % Amplitude factor = 1 initially. damped later to values <1

 while max(abs(dal)) >= Conv || AF ~= 1% % Convergence criterium

dal = zeros(nf,1); %Reset convergence measure

for j = [RunJ(RunJ >= round(nf/2))]
    %Set up LHS
    omega = StabRes.omegavec(j);
    beta = StabRes.betavec(j);
    
if mmry ~=1 % Run LHS generation only in the first iteration
    mmry = 1; % Set memory to 1 so this part of the code is not repeated
    
    fprintf('Constructing LHS matrices for all modes ... \n')
    perc = 0; % Progress variable (percentage)
    fprintf('%.3g percent ',perc)
    tic
for i=1:Grid.nx

 
 %% Prepare LHS building

 % Create diagonal matrix form of base flow 
      U = diag(BF.Ur(:,i));
    dxU = diag(BF.dxUr(:,i));
    dyU = diag(BF.dyUr(:,i));
    
      V = diag(BF.Vr(:,i));
    dxV = diag(BF.dxVr(:,i)); 
    dyV = diag(BF.dyVr(:,i));
    
      W = diag(BF.Wr(:,i));
    dxW = diag(BF.dxWr(:,i));
    dyW = diag(BF.dyWr(:,i));

 % Local transformations
    xix  = StabGrid.xix(:,i);
    etax = StabGrid.etax(:,i);
    xiy  = StabGrid.xiy(:,i);
    etay = StabGrid.etay(:,i);
    xixx = StabGrid.xixx(:,i);
    xiyy = StabGrid.xiyy(:,i);
    etaxx= StabGrid.etaxx(:,i);
    etayy= StabGrid.etayy(:,i);
            
 % Create common terms (ct) for matrices A-F
    Act =  U.*etax*D1...
          +V.*etay*D1...
          -1/BF.Re.*(etax.^2+etay.^2).*D2...
          -1/BF.Re.*(etaxx+etayy).*D1;
        
    Ect =  U.*xix...
          +V.*xiy...
          -1/BF.Re*(xixx+xiyy).*I...
          -1/BF.Re*(2*etax.*xix+2*etay.*xiy).*D1;

%% construct coefficient matrices
    Ai=[(Act+dxU).*bufs(:,i)       dyU.*bufs(:,i)        Z    etax.*D1.*bufs(:,i)
            dxV.*bufs(:,i)        (Act+dyV).*bufs(:,i)   Z    etay.*D1.*bufs(:,i)
            dxW.*bufs(:,i)         dyW.*bufs(:,i)       Act.*bufs(:,i) Z  
           (etax.*D1).*bufsp(:,i) (etay.*D1).*bufsp(:,i) Z             Z];        

    Bi=[-iu*I.*bufs(:,i)  Z  Z  Z
         Z  -iu*I.*bufs(:,i) Z  Z
         Z  Z  -iu*I.*bufs(:,i) Z
         Z  Z                Z  Z];
     
    Ci=[iu*W.*bufs(:,i)    Z   Z    Z
         Z    iu*W.*bufs(:,i)  Z    Z
         Z    Z   iu*W.*bufs(:,i)   iu*I.*bufs(:,i)
         Z    Z   iu*I.*bufsp(:,i)  Z];
     
    Di=[I/BF.Re.*bufs(:,i)  Z  Z  Z
         Z  I/BF.Re.*bufs(:,i) Z  Z
         Z  Z  I/BF.Re.*bufs(:,i) Z
         Z                      Z  Z  Z];

    Ei=[Ect.*bufs(:,i)      Z        Z    I.*xix.*bufs(:,i)
        Z       Ect.*bufs(:,i)       Z    I.*xiy.*bufs(:,i)
        Z       Z        Ect.*bufs(:,i)            Z
        I.*xix.*bufsp(:,i)   I.*xiy.*bufsp(:,i) Z  Z];

    Fi=[-I/BF.Re.*(xix.^2+xiy.^2).*bufs(:,i)   Z       Z       Z
        Z       -I/BF.Re.*(xix.^2+xiy.^2).*bufs(:,i)   Z       Z
        Z       Z       -I/BF.Re.*(xix.^2+xiy.^2).*bufs(:,i)   Z
        Z       Z                                          Z       Z];
    
    
    % Ensure pressure is extrapolated downward from the surface
    Pres = I.*(1-bufsp(:,i)) + diag(-1*(1-bufsp(2:end,i)),-1);

    % Set up buffer and embedded boundary method
    Gi=[I.*(1-bufs(:,i))+I.*bufs(:,end)*(1-bufc(i)) Z Z Z
        Z I.*(1-bufs(:,i)) Z Z
        Z Z I.*(1-bufs(:,i))+I.*bufs(:,end)*(1-bufc(i))  Z 
        Z Z Z Pres];


 
    % construct the ML block for station i

    sten=4*Grid.ny*(i-1)+1:4*Grid.ny*i; % row coordinates of FD stencil
    sten1=1:4*Grid.ny;                  % stencil height

    if i==1
        % At i = 1 LHS is the identity matrix to force 1*q = R. The inflow 
        % BC is presented in R such that the solution is the inflow BC.

        ML1MFD(sten1,sten)=eye(4*Grid.ny); % identity matrix for enforcing initial condition
        ML1(sten1,sten)=eye(4*Grid.ny);    % identity matrix for enforcing initial condition
        ML2(sten1,sten)=0*eye(4*Grid.ny);  % zero, 1's already in ML1
        ML3(sten1,sten)=0*eye(4*Grid.ny);  % zero, 1's already in ML1
        ML4(sten1,sten)=0*eye(4*Grid.ny);  % zero, 1's already in ML1

    elseif (i==2)
    
        stenr=[sten-4*Grid.ny sten sten+4*Grid.ny sten+8*Grid.ny sten+12*Grid.ny];    % row stencil
        
        ML1MFD(sten1,stenr)=[-3*Ei/(12*dxi)+11*Fi/(12*dxi^2), ...
                         (Ai-10*Ei/(12*dxi)-20*Fi/(12*dxi^2))+Gi/bufc(i),...
                             18*Ei/(12*dxi)+ 6*Fi/(12*dxi^2), ...
                            - 6*Ei/(12*dxi)+ 4*Fi/(12*dxi^2), ...
                                Ei/(12*dxi)- 1*Fi/(12*dxi^2)]*bufc(i); 
                            
        ML1(sten1,stenr)=ML1MFD(sten1,stenr);
        ML2(sten1,sten)=Bi*bufc(i);
        ML3(sten1,sten)=Ci*bufc(i);
        ML4(sten1,sten)=Di*bufc(i);

        % apply boundary condition
        for k=1:3
        if k ~= 2 % MFD v should be free at top BC to allow for BL growth
        ML1MFD((k-1)*Grid.ny+1,:)=0;
        ML1MFD((k-1)*Grid.ny+1,4*Grid.ny*(i-1)+(k-1)*Grid.ny+1)=1;
        end
        ML1((k-1)*Grid.ny+1,:)=0;  ML2((k-1)*Grid.ny+1,:)=0; ML3((k-1)*Grid.ny+1,:)=0; ML4((k-1)*Grid.ny+1,:)=0;   
        ML1((k-1)*Grid.ny+1,4*Grid.ny*(i-1)+(k-1)*Grid.ny+1)=1;    
        
        ML1MFD((k)*Grid.ny,:)=0; ML1((k)*Grid.ny,:)=0;  ML2((k)*Grid.ny,:)=0; ML3((k)*Grid.ny,:)=0; ML4((k)*Grid.ny,:)=0;  
        ML1MFD((k)*Grid.ny,4*Grid.ny*(i-1)+(k)*Grid.ny)=1;    
        ML1((k)*Grid.ny,4*Grid.ny*(i-1)+(k)*Grid.ny)=1;  
        end

    elseif (i==Grid.nx-1)

        sten1=1:4*Grid.ny;                  % stencil height
        stenr=[sten-8*Grid.ny sten-4*Grid.ny sten sten+4*Grid.ny];    % row coordinates of FD stencil
        
        ML1MFD(sten1,stenr)=[Ei/(6*dxi), ...
                          -6*Ei/(6*dxi)+  Fi/(dxi^2),...
                       (Ai+3*Ei/(6*dxi)-2*Fi/(dxi^2))+Gi/bufc(i),...
                           2*Ei/(6*dxi)+  Fi/(dxi^2)]*bufc(i);
                       
        ML1(sten1,stenr)=ML1MFD(sten1,stenr);
        ML2(sten1,sten)=Bi*bufc(i);
        ML3(sten1,sten)=Ci*bufc(i);
        ML4(sten1,sten)=Di*bufc(i);

        % apply boundary condition
        for k=1:3
        if k ~= 2 % MFD v should be free at top BC to allow for BL growth
        ML1MFD((k-1)*Grid.ny+1,:)=0;
        ML1MFD((k-1)*Grid.ny+1,4*Grid.ny*(i-1)+(k-1)*Grid.ny+1)=1;
        end
        ML1((k-1)*Grid.ny+1,:)=0;  ML2((k-1)*Grid.ny+1,:)=0; ML3((k-1)*Grid.ny+1,:)=0; ML4((k-1)*Grid.ny+1,:)=0;   
        ML1((k-1)*Grid.ny+1,4*Grid.ny*(i-1)+(k-1)*Grid.ny+1)=1;    
        
        ML1MFD((k)*Grid.ny,:)=0; ML1((k)*Grid.ny,:)=0;  ML2((k)*Grid.ny,:)=0; ML3((k)*Grid.ny,:)=0; ML4((k)*Grid.ny,:)=0;  
        ML1MFD((k)*Grid.ny,4*Grid.ny*(i-1)+(k)*Grid.ny)=1;    
        ML1((k)*Grid.ny,4*Grid.ny*(i-1)+(k)*Grid.ny)=1;  
        end

    elseif (i==Grid.nx)
        sten1=1:4*Grid.ny;                  % stencil height
        stenr=[sten-8*Grid.ny sten-4*Grid.ny sten];   % row coordinates of FD stencil
        
        ML1MFD(sten1,stenr)=[Ei/(2*dxi)+  Fi/(dxi^2), ...
                          -4*Ei/(2*dxi)-2*Fi/(dxi^2), ...
                       (Ai+3*Ei/(2*dxi)+  Fi/(dxi^2))+Gi/bufc(i)]*bufc(i);
                  
        ML1(sten1,stenr)=ML1MFD(sten1,stenr);
        ML2(sten1,sten)=Bi*bufc(i);
        ML3(sten1,sten)=Ci*bufc(i);
        ML4(sten1,sten)=Di*bufc(i);

        % apply boundary condition
        for k=1:3
        if k ~= 2 % MFD v should be free at top BC to allow for BL growth
        ML1MFD((k-1)*Grid.ny+1,:)=0;
        ML1MFD((k-1)*Grid.ny+1,4*Grid.ny*(i-1)+(k-1)*Grid.ny+1)=1;
        end
        % ML
        ML1((k-1)*Grid.ny+1,:)=0;  ML2((k-1)*Grid.ny+1,:)=0; ML3((k-1)*Grid.ny+1,:)=0; ML4((k-1)*Grid.ny+1,:)=0;   
        ML1((k-1)*Grid.ny+1,4*Grid.ny*(i-1)+(k-1)*Grid.ny+1)=1;    
        
        ML1MFD((k)*Grid.ny,:)=0; ML1((k)*Grid.ny,:)=0;  ML2((k)*Grid.ny,:)=0; ML3((k)*Grid.ny,:)=0; ML4((k)*Grid.ny,:)=0;  
        ML1MFD((k)*Grid.ny,4*Grid.ny*(i-1)+(k)*Grid.ny)=1;    
        ML1((k)*Grid.ny,4*Grid.ny*(i-1)+(k)*Grid.ny)=1;  
        end
    elseif i == stepi-1 && Grid.StepH ~= 0
        
        % Values under step height
        sten1=[0*Grid.ny+((yi+1):Grid.ny) 1*Grid.ny+((yi+1):Grid.ny) 2*Grid.ny+((yi+1):Grid.ny) 3*Grid.ny+((yi+1):Grid.ny) ];
        stenr=[sten-12*Grid.ny sten-8*Grid.ny sten-4*Grid.ny sten sten+4*Grid.ny];    % row coordinates of FD stencil
        
        ML1MFD(sten1,stenr)=[- 1*Ei(sten1,:)/(12*dxi)-   Fi(sten1,:)/(12*dxi^2),...
                               6*Ei(sten1,:)/(12*dxi)+ 4*Fi(sten1,:)/(12*dxi^2),...
                             -18*Ei(sten1,:)/(12*dxi)+ 6*Fi(sten1,:)/(12*dxi^2),...
                  Ai(sten1,:)+10*Ei(sten1,:)/(12*dxi)-20*Fi(sten1,:)/(12*dxi^2)+Gi(sten1,:)/bufc(i),...
                               3*Ei(sten1,:)/(12*dxi)+11*Fi(sten1,:)/(12*dxi^2)]*bufc(i);
        
        ML1(sten1,stenr)=ML1MFD(sten1,stenr);
        ML2(sten1,sten)=Bi(sten1,:)*bufc(i);
        ML3(sten1,sten)=Ci(sten1,:)*bufc(i);
        ML4(sten1,sten)=Di(sten1,:)*bufc(i);
        
        % Values above step
        sten1=[0*Grid.ny+(1:yi) 1*Grid.ny+(1:yi) 2*Grid.ny+(1:yi) 3*Grid.ny+(1:yi) ];
        stenr=[sten-8*Grid.ny sten-4*Grid.ny sten sten+4*Grid.ny sten+8*Grid.ny];    % row stencil
        
        ML1MFD(sten1,stenr)=[   Ei(sten1,:)/(12*dxi)-   Fi(sten1,:)/(12*dxi^2),... 
                             -8*Ei(sten1,:)/(12*dxi)+16*Fi(sten1,:)/(12*dxi^2),...
                (Ai(sten1,:)                        -30*Fi(sten1,:)/(12*dxi^2))+Gi(sten1,:)/bufc(i),... 
                              8*Ei(sten1,:)/(12*dxi)+16*Fi(sten1,:)/(12*dxi^2),...
                             -  Ei(sten1,:)/(12*dxi)-   Fi(sten1,:)/(12*dxi^2)]*bufc(i);
                         
        ML1(sten1,stenr)=ML1MFD(sten1,stenr);
        ML2(sten1,sten)=Bi(sten1,:)*bufc(i);
        ML3(sten1,sten)=Ci(sten1,:)*bufc(i);
        ML4(sten1,sten)=Di(sten1,:)*bufc(i);    

        % apply boundary condition
        for k=1:3
        if k ~= 2 % MFD v should be free at top BC to allow for BL growth
        ML1MFD((k-1)*Grid.ny+1,:)=0;
        ML1MFD((k-1)*Grid.ny+1,4*Grid.ny*(i-1)+(k-1)*Grid.ny+1)=1;
        end
        ML1((k-1)*Grid.ny+1,:)=0;  ML2((k-1)*Grid.ny+1,:)=0; ML3((k-1)*Grid.ny+1,:)=0; ML4((k-1)*Grid.ny+1,:)=0;   
        ML1((k-1)*Grid.ny+1,4*Grid.ny*(i-1)+(k-1)*Grid.ny+1)=1;    
        
        ML1MFD((k)*Grid.ny,:)=0; ML1((k)*Grid.ny,:)=0;  ML2((k)*Grid.ny,:)=0; ML3((k)*Grid.ny,:)=0; ML4((k)*Grid.ny,:)=0;  
        ML1MFD((k)*Grid.ny,4*Grid.ny*(i-1)+(k)*Grid.ny)=1;    
        ML1((k)*Grid.ny,4*Grid.ny*(i-1)+(k)*Grid.ny)=1;  
        end
        
    elseif i == stepi && Grid.StepH ~= 0

        % Values under step height
        sten1=[0*Grid.ny+((yi+1):Grid.ny) 1*Grid.ny+((yi+1):Grid.ny) 2*Grid.ny+((yi+1):Grid.ny) 3*Grid.ny+((yi+1):Grid.ny) ];
        stenr=[sten-16*Grid.ny sten-12*Grid.ny sten-8*Grid.ny sten-4*Grid.ny sten];   % row coordinates of FD stencil
        
        ML1MFD(sten1,stenr)=[+ 3*Ei(sten1,:)/(12*dxi)+ 11*Fi(sten1,:)/(12*dxi^2),...
                             -16*Ei(sten1,:)/(12*dxi)- 56*Fi(sten1,:)/(12*dxi^2),...
                             +36*Ei(sten1,:)/(12*dxi)+114*Fi(sten1,:)/(12*dxi^2),...
                             -48*Ei(sten1,:)/(12*dxi)-104*Fi(sten1,:)/(12*dxi^2),... 
                 (Ai(sten1,:)+25*Ei(sten1,:)/(12*dxi)+ 35*Fi(sten1,:)/(12*dxi^2))+Gi(sten1,:)/bufc(i)]*bufc(i);
            
        ML1(sten1,stenr)=ML1MFD(sten1,stenr);
        ML2(sten1,sten)=Bi(sten1,:)*bufc(i);
        ML3(sten1,sten)=Ci(sten1,:)*bufc(i);
        ML4(sten1,sten)=Di(sten1,:)*bufc(i);
 
        % Values above step
        sten1=[0*Grid.ny+(1:yi) 1*Grid.ny+(1:yi) 2*Grid.ny+(1:yi) 3*Grid.ny+(1:yi) ];
        stenr=[sten-8*Grid.ny sten-4*Grid.ny sten sten+4*Grid.ny sten+8*Grid.ny];   % row coordinates of FD stencil
        
        ML1MFD(sten1,stenr)=[Ei(sten1,:)/(12*dxi)- 1*Fi(sten1,:)/(12*dxi^2), ...
                          -8*Ei(sten1,:)/(12*dxi)+16*Fi(sten1,:)/(12*dxi^2), ...
             (Ai(sten1,:)                        -30*Fi(sten1,:)/(12*dxi^2))+Gi(sten1,:)/bufc(i), ...
                           8*Ei(sten1,:)/(12*dxi)+16*Fi(sten1,:)/(12*dxi^2), ...
                          -1*Ei(sten1,:)/(12*dxi)- 1*Fi(sten1,:)/(12*dxi^2)]*bufc(i);
                      
        ML1(sten1,stenr)=ML1MFD(sten1,stenr);
        ML2(sten1,sten)=Bi(sten1,:)*bufc(i);
        ML3(sten1,sten)=Ci(sten1,:)*bufc(i);
        ML4(sten1,sten)=Di(sten1,:)*bufc(i); 

        % apply boundary condition
        for k=1:3
        if k ~= 2 % MFD v should be free at top BC to allow for BL growth
        ML1MFD((k-1)*Grid.ny+1,:)=0;
        ML1MFD((k-1)*Grid.ny+1,4*Grid.ny*(i-1)+(k-1)*Grid.ny+1)=1;
        end
        ML1((k-1)*Grid.ny+1,:)=0;  ML2((k-1)*Grid.ny+1,:)=0; ML3((k-1)*Grid.ny+1,:)=0; ML4((k-1)*Grid.ny+1,:)=0;   
        ML1((k-1)*Grid.ny+1,4*Grid.ny*(i-1)+(k-1)*Grid.ny+1)=1;    
        
        ML1MFD((k)*Grid.ny,:)=0; ML1((k)*Grid.ny,:)=0;  ML2((k)*Grid.ny,:)=0; ML3((k)*Grid.ny,:)=0; ML4((k)*Grid.ny,:)=0;  
        ML1MFD((k)*Grid.ny,4*Grid.ny*(i-1)+(k)*Grid.ny)=1;    
        ML1((k)*Grid.ny,4*Grid.ny*(i-1)+(k)*Grid.ny)=1;  
        end
    else
        sten1=1:4*Grid.ny;                  % stencil height
        stenr=[sten-8*Grid.ny sten-4*Grid.ny sten sten+4*Grid.ny sten+8*Grid.ny];    % row coordinates of FD stencil

        ML1MFD(sten1,stenr)=[Ei/(12*dxi)- 1*Fi/(12*dxi^2), ...
                          -8*Ei/(12*dxi)+16*Fi/(12*dxi^2), ...
                        (Ai             -30*Fi/(12*dxi^2))+Gi/bufc(i), ...
                           8*Ei/(12*dxi)+16*Fi/(12*dxi^2), ...
                            -Ei/(12*dxi)- 1*Fi/(12*dxi^2)]*bufc(i);
         
        ML1(sten1,stenr)=ML1MFD(sten1,stenr);
        ML2(sten1,sten)=Bi*bufc(i);
        ML3(sten1,sten)=Ci*bufc(i);
        ML4(sten1,sten)=Di*bufc(i);    

        % apply boundary condition
        for k=1:3
        if k ~= 2 % MFD v should be free at top BC to allow for BL growth
        ML1MFD((k-1)*Grid.ny+1,:)=0;
        ML1MFD((k-1)*Grid.ny+1,4*Grid.ny*(i-1)+(k-1)*Grid.ny+1)=1;
        end
        ML1((k-1)*Grid.ny+1,:)=0;  ML2((k-1)*Grid.ny+1,:)=0; ML3((k-1)*Grid.ny+1,:)=0; ML4((k-1)*Grid.ny+1,:)=0;   
        ML1((k-1)*Grid.ny+1,4*Grid.ny*(i-1)+(k-1)*Grid.ny+1)=1;    
        
        ML1MFD((k)*Grid.ny,:)=0; ML1((k)*Grid.ny,:)=0;  ML2((k)*Grid.ny,:)=0; ML3((k)*Grid.ny,:)=0; ML4((k)*Grid.ny,:)=0;  
        ML1MFD((k)*Grid.ny,4*Grid.ny*(i-1)+(k)*Grid.ny)=1;    
        ML1((k)*Grid.ny,4*Grid.ny*(i-1)+(k)*Grid.ny)=1;  
        end
    end
 
    % Update progress percentage
    perc=round(i*100/Grid.nx);

    if perc<10
        fprintf('\b\b\b\b\b\b\b\b\b\b%.3g percent ',perc)
    elseif perc<100
        fprintf('\b\b\b\b\b\b\b\b\b\b\b%.3g percent ',perc)
    else
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b%.3g percent ',perc)
    end

% Store coordinate (iML,jML) and value (s) data for LHS
[iML1MFD,jML1MFD,sML1MFD] = find(ML1MFD); iML1MFD = iML1MFD+4*Grid.ny*(i-1);
[iML1,jML1,sML1] = find(ML1); iML1= iML1+4*Grid.ny*(i-1);
[iML2,jML2,sML2] = find(ML2); iML2= iML2+4*Grid.ny*(i-1);
[iML3,jML3,sML3] = find(ML3); iML3= iML3+4*Grid.ny*(i-1);
[iML4,jML4,sML4] = find(ML4); iML4= iML4+4*Grid.ny*(i-1);

% iML locations
iML1MFDtotal = [iML1MFDtotal; iML1MFD];
iML1total = [iML1total; iML1];
iML2total = [iML2total; iML2];
iML3total = [iML3total; iML3];
iML4total = [iML4total; iML4];

% jML locations
jML1MFDtotal = [jML1MFDtotal; jML1MFD];
jML1total = [jML1total; jML1];
jML2total = [jML2total; jML2];
jML3total = [jML3total; jML3];
jML4total = [jML4total; jML4];

% ML values
sML1MFDtotal = [sML1MFDtotal; sML1MFD];
sML1total = [sML1total; sML1];
sML2total = [sML2total; sML2];
sML3total = [sML3total; sML3];
sML4total = [sML4total; sML4];
    
% Reset building blocks
ML1MFD=0*ML1MFD; ML1=ML1*0; ML2=ML2*0; ML3=ML3*0; ML4=ML4*0;

end
clear ML1MFD ML1 ML2 ML3 ML4

% Create LHS from coordinate data
MLAMFD = sparse(iML1MFDtotal,jML1MFDtotal,sML1MFDtotal,4*Grid.ny*Grid.nx,4*Grid.ny*Grid.nx);
MLA1 = sparse(iML1total,jML1total,sML1total,4*Grid.ny*Grid.nx,4*Grid.ny*Grid.nx);
if any(StabRes.omegavec)
MLA2 = sparse(iML2total,jML2total,sML2total,4*Grid.ny*Grid.nx,4*Grid.ny*Grid.nx);
end
if any(StabRes.betavec)
MLA3 = sparse(iML3total,jML3total,sML3total,4*Grid.ny*Grid.nx,4*Grid.ny*Grid.nx);
MLA4 = sparse(iML4total,jML4total,sML4total,4*Grid.ny*Grid.nx,4*Grid.ny*Grid.nx);
end

% Clear up workspace
clear iML1MFDtotal jML1MFDtotal sML1MFDtotal iML1total jML1total sML1total ...
    iML2total jML2total sML2total iML3total jML3total sML3total iML3total ...
    jML3total sML3total iML1 iML1MFD iML2 iML3 iML4 jML1 jML2 jML3 jML4 sML1 ...
    sML2 sML3 sML4 jML4total jML1MFD iML4total jML1MFD jML4total sML4total ...
    sML1MFD A1 A2 A3 A4 Ai Bi Ci Di Ei Fi dxU dxBL.U dxBL.V dxW dxBL.W Gi ...
    stenr sten sten1 U BL.U V BL.V W BL.W 

fprintf(', finished in %.2f seconds.\n',toc)
end

% Construct ML 
if omega == 0 && beta == 0  % MFD ignore beta and omega additions
 ML = MLAMFD; 
elseif omega == 0           % ignore omega additions
 ML = MLA1+MLA3*beta+MLA4*beta^2;
elseif beta == 0            % ignore beta additions
 ML = MLA1+MLA2*omega;      
else                        % add all  
 ML = MLA1+MLA2*omega+MLA3*beta+MLA4*beta^2;   
end

%Construct R
R = 0*R; % Reset RHS
for i = 1:Grid.nx
  if i==1
        if ismember(j,nonzeromodes)
            R(4*Grid.ny*(i-1)+1:4*Grid.ny*i,1)=AF*phiIC(j,:); % Inflow BC is inflow shape corrected for ramping amplitude factor AF
        else   
            R(4*Grid.ny*(i-1)+1:4*Grid.ny*i,1)=0*phi(j,:,1);  % Force inflow of other modes to be 0 
        end
  else
    for k = 1:3
        R(4*Grid.ny*(i-1)+(k-1)*Grid.ny+1) = bct(k,i,j); % top bc_k (k1->u,k2->v,k3->w)
        R(4*Grid.ny*(i-1)+(k)*Grid.ny) = bcw(k,i,j);     % wall bc_k (k1->u,k2->v,k3->w) 
    end 
  end
end

%%
tic
fprintf([' Solving for mode (m,n) = (' num2str(Mvec(j)) ',' num2str(Nvec(j)) ')...'])

% Ensure inhomogeneous BC's are not affected by nonlinear terms
f(j,R~=0)= 0; % f has to be 0 whenever R is not

% Solve the problem via backslash operator (LU operator employed)
q(j,:)=ML\(R+f(j,:).');

% Solution of the conjugate is the conjugate of the solution
q(nf+1-j,:)=conj(q(j,:));

clear ML

fprintf('\b\b\b, finished in %.2f seconds.\n',toc)


end %j-loop
%% damp inflow amplitude

if nnz(Stab.A0)~=1 % Does not reduce inflow amplitude in linear problems
if k1 == 0 % Only perform this after the first iteration

    % Reshape q into phi
    for j = RunJ
    for i = 1:Grid.nx
    StabRes.phi(j,:,i)=q(j,(i-1)*4*Grid.ny+(1:4*Grid.ny));
    end
    end
    
    % Calculate the maximum amplitude reached linearly
    [AVAL,~] = max(max(abs(squeeze(StabRes.phi(round(nf/2)+1,1:3*Grid.ny,:)))));
    
    % Find Amplitude Factor (AF) such that evolution is linear
    if AVAL >= Opt.AMAX % If maximum linear amplitude > threshold
        AF = Opt.AMAX/AVAL; % AF is the amplitude factor such that qmax = AMAX
    end
    
    q = AF*q; % Multiply linear result by AF before calculation of NLT

end 
end

%% Update u,v,w from solution q
for j = RunJ
   for i = 1:Grid.nx
       StabRes.u(j,:,i) = q(j,           (i-1)*4*Grid.ny+(1:Grid.ny));
       StabRes.v(j,:,i) = q(j,Grid.ny+   (i-1)*4*Grid.ny+(1:Grid.ny));
       StabRes.w(j,:,i) = q(j,2*Grid.ny+ (i-1)*4*Grid.ny+(1:Grid.ny));
       StabRes.p(j,:,i) = q(j,3*Grid.ny+ (i-1)*4*Grid.ny+(1:Grid.ny));  
   end
end

%% Determine active modes and allow select new ones to enter

%Find Amplitudes for all modes over xi
% Update active modes
for i = 1:Grid.nx
   for j = RunJ
      StabRes.A(j,i) = max(abs((q(j,(i-1)*Grid.ny*4 + (1:Grid.ny) ))));
   end
end

% Find amplitudes of RHS mode interactions over entire domain
for i = 1:Grid.nx
Amat(:,:,i) = StabRes.A(:,i)*StabRes.A(:,i)';
end

% Find maximum in xi
Amat = max(Amat,[],3);

% Find location where the maximum amplitude of the primary mode occurs
[~,index] = max(StabRes.A(:,:),[],2);
index = max(index); % Need only one stage, take the most downstream one


% Determine which modes are allowed to enter
if Opt.Sweep == 1 % If intermedate results are requested, allow all modes to enter the system while ramping.
    [nzA, ] = find((StabRes.A(:,index))); %find nonzero Amplitudes
    nfm = max(nzA)+(1);     %find index of highest harmonic+1
    nfl = min(nzA)-(1);     %find index of conjugate of highest harmonic-1 
else
    % If amplitude ramping finished (AF=1) and converged introduce higher harmonics
    if AF>=1 && max(abs(dalsave(:,max(k1,1))))<=Conv*100 && ~isnan(dalsave(max(RunJ),k1))
    [nzA, ] = find((StabRes.A(:,index))); %find nonzero Amplitudes
    nfm = max(nzA)+(1);     %find index of highest harmonic+1
    nfl = min(nzA)-(1);     %find index of conjugate of highest harmonic-1  

    else %AF <1, Still ramping up, do not let more than 2 modes enter the system
        if any(Stab.A0) % Modes were introduced at the inflow (LOAD or ILST)
        [nzA, ] = find(Stab.A0(:)); %find nonzero Amplitudes
        nfm = max(nzA)+(1);     %find index of highest harmonic+1
        nfl = min(nzA)-(1);     %find index of conjugate of highest harmonic-1  
    
        else %Blowing and suction was active for mode introduction (WALL)
        [nzA, ] = [round(nf/2)-1 round(nf/2)+1]; %find nonzero Amplitudes
        nfm = max(nzA)+(1);     %find index of highest harmonic+1
        nfl = min(nzA)-(1);     %find index of conjugate of highest harmonic-1    
        end
    end
end

% Reset Amode
Amode = zeros(size(StabRes.A(:,1)));

% Predict amplitude of to be introduced modes
for j=1:nf
    if j > nfm || j<nfl
       Amode(j) = 0; 
    else
       Amode(j)  = sum(sum(Amat.*HB(:,:,j)));     
    end
end

% Set active modes as either currently or forced greater than threshold
% Save currently active modes
Aactive_old = Aactive; 

% Active modes are either already introduced or are expected to cross the
% threshold TH
for j = 1:nf
Aactive(j) = any((StabRes.A(j,:) + Amode(j,:)) >= TH); 
end

% Find which modes are new
NewModes = find(Aactive - Aactive_old > 0);

% Modes to run are currently active modes and new ones
RunJ = sort([RunJ NewModes]);



%% Calculate Nonlinear terms

[ f ] = NLT_HNS(StabGrid,RunJ,StabRes,D1,HB);

% Apply attenuation to forcing terms
for jj = RunJ 
for i =ibuff:Grid.nx
    f(jj,(i-1)*4*Grid.ny+(1:4*Grid.ny)) = f(jj,(i-1)*4*Grid.ny+(1:4*Grid.ny))*bufcf(i);
end
end

%% Intermediate plotting
for j = RunJ
for i = 1:Grid.nx
    phi(j,:,i)=q(j,(i-1)*4*Grid.ny+(1:4*Grid.ny));
    f_phi(j,:,i) = f(j,(i-1)*4*Grid.ny+(1:4*Grid.ny));
end
end

figure(473)
hold off
for j = RunJ(RunJ>=round(nf/2))
    semilogy(StabGrid.xiun(1:round(Grid.nx)),squeeze(max(abs(f_phi(j,1:Grid.ny,1:round(Grid.nx)))))) 
    hold on
end
ylim([1e-12 1])
title('forcing term')

figure(474)
subplot(4,1,1:3)
hold off
for j = RunJ(RunJ>=round(nf/2))
    semilogy(StabGrid.xiun(1:round(Grid.nx)),squeeze(max(abs(phi(j,1:Grid.ny,1:round(Grid.nx)))))) 
    hold on
end
plot(StabGrid.xiun(round(Grid.nx*Opt.xb/100))*ones(1,2),[1e-10 1])
xlim([StabGrid.xiun(1) StabGrid.xiun(end)])
ylim([1e-12 1])
title('Amplitude evolution')
set(gcf, 'Position', [650 400 500 400]);
ylabel('A')

subplot(4,1,4)
plot(StabGrid.xiun,bufc)
xlim([StabGrid.xiun(1) StabGrid.xiun(end)])
xlabel('x')



%% Convergence
% Calculate convergence measure, Change based on the paper
for j = RunJ
    dal(j) = sum(abs(StabRes.A(j,1:ibuf)-Aold(j,1:ibuf)))/max(StabRes.A(j,1:ibuf))/Grid.nx;
end

% End simulation if only 1 mode is present (Linear)
if length(RunJ) == 1 ; dal = zeros(size(dal)); end

% Save amplitudes
Aold = StabRes.A;

%Only increase initial amplitude if current amplitude is quasi converged
if k1 ~= 0
   % If Opt.Sweep = 1, intermediate results should be converged fully
   % Else, convergence requirement for intermediate results is reduced
if max(abs(dal))<=(1-Opt.Sweep)*Opt.Conv*Opt.ConvF+Opt.Sweep*Opt.Conv
   AF = AF*Opt.AFg; % Increase Amplitude at the inflow
   if ge(AF,1); AF = 1; end % Amplitude should not exceed desired amplitude

   % Initialize list of inflow amplitudes and nonlinear forcing
   if ~exist('AFlist','var')
       AFlist = [AF AF];
       flist = f;
       flist(:,:,size(flist,3)+1) = f;
   end

   % Append list of inflow amplitudes and forcing for extrapolation purposes
   AFlist = [AFlist AF];
   flist(:,:,size(flist,3)+1) = f;

   %Linearly (in A^2) extrapolate the forcing term
   if AFlist(end)~=AFlist(end-1) && AFlist(end-1)~=AFlist(end-2) % Extrapolation does not work otherwise
   f = flist(:,:,end)+(flist(:,:,end)-flist(:,:,end-1))/(AFlist(end-1)^2-AFlist(end-2)^2) * (AFlist(end)^2-AFlist(end-1)^2);
   end 
end
end

%% Intermediate results 
% If the sweep option is on, additional intermediate results are saved
% seperately for amplitude, alpha and shape functions denoted by the sweep
% suffix. The last entry (...,k1) is the iteration number. 

k1 = k1+1; 

if Opt.Sweep == 1 && max(abs(dal))<=Opt.Sweep*Opt.Conv
k2 = k2+1;
% Store amplitude factor AF
StabRes.AF(k2) = AF;

% Prepare output via shape function normalization

% Calculate du/dxi
for j = RunJ
    dux(j,:,:) = - FD1d4o( squeeze(StabRes.u(j,:,:)), dxi );
end

% Create y integration vector
    yinteg=linspace(0,max(StabGrid.etaun),4000);

for j = RunJ
for i = 1:Grid.nx    
    
    [~,yimax]=max(abs(interp1(StabGrid.etaun,StabRes.u(j,:,i),yinteg,'spline')));
    dumax=interp1(StabGrid.etaun,dux(j,:,i),yinteg(yimax),'spline');
    umax=interp1(StabGrid.etaun,StabRes.u(j,:,i),yinteg(yimax),'spline');

    % Calculate an alpha a posteriori for comparison purposes and save
    % current iteration
    StabRes.alphasweep(j,i,k2) = (dumax/umax)/iu;
    
    % Normalize shape function by local u_max and save current iteration
    StabRes.usweep(j,:,i,k2) = StabRes.u(j,:,i)./abs(umax);
    StabRes.vsweep(j,:,i,k2) = StabRes.v(j,:,i)./abs(umax);
    StabRes.wsweep(j,:,i,k2) = StabRes.w(j,:,i)./abs(umax);
    StabRes.psweep(j,:,i,k2) = StabRes.p(j,:,i)./abs(umax);

    % Save current iteration phi results
    StabRes.phisweep(j,:,i,k2)=StabRes.phi(j,:,i);
    
    % Store amplitudes based on u in StabRes.A
    if j == round(nf/2)
    StabRes.Asweep(j,i,k2) = abs(umax); 
    else
    StabRes.Asweep(j,i,k2) = 2*abs(umax); % Double amplitude as conjugates will be deleted
    end
end
end



end
%% Plot convergence

dalsave(:,k1) = abs(dal);

dalsave(dalsave <= 0) = NaN;
dalsave(dalsave == 1) = NaN;

figure(15)
hold off
semilogy(linspace(1,k1,k1),dalsave(RunJ,:))
hold on
semilogy(linspace(1,k1,k1),Conv.*ones(1,k1),'k--')
set(gcf, 'Position', [100 400 500 400]);
xlim([0 k1])
ylabel('Relative error')
xlabel('Iteration')
pause(0.05)
end %convergence loop

%% Prepare output via shape function normalization
% Calculate du/dxi
for j = RunJ
dux(j,:,:) = - FD1d4o( squeeze(StabRes.u(j,:,:)), dxi );
end

% Create y integration vector
yinteg=linspace(0,max(StabGrid.etaun),4000);

for j = RunJ
for i = 1:Grid.nx    
    
    [~,k]=max(abs(interp1(StabGrid.etaun,StabRes.u(j,:,i),yinteg,'spline')));
    dumax=interp1(StabGrid.etaun,dux(j,:,i),yinteg(k),'spline');
    umax=interp1(StabGrid.etaun,StabRes.u(j,:,i),yinteg(k),'spline');

    % Calculate an alpha a posteriori for comparison purposes
    StabRes.alpha(j,i) = (dumax/umax)/iu;
    
    %Normalize shape function by local u_max
    StabRes.u(j,:,i) = StabRes.u(j,:,i)./abs(umax);
    StabRes.v(j,:,i) = StabRes.v(j,:,i)./abs(umax);
    StabRes.w(j,:,i) = StabRes.w(j,:,i)./abs(umax);
    StabRes.p(j,:,i) = StabRes.p(j,:,i)./abs(umax);
    
    % Store amplitudes based on u in StabRes.A
    if j == round(nf/2)
    StabRes.A(j,i) = abs(umax); 
    else
    StabRes.A(j,i) = 2*abs(umax); % Double amplitude as conjugates will be deleted
    end
end
end

% Delete conjugates from results
StabRes.omegavec(:,1:round(nf/2)-1) = [];
StabRes.betavec(:,1:round(nf/2)-1)  = [];
StabRes.A(1:round(nf/2)-1,:)        = [];
StabRes.alpha(1:round(nf/2)-1,:)    = [];
StabRes.u(1:round(nf/2)-1,:,:)      = [];
StabRes.v(1:round(nf/2)-1,:,:)      = [];
StabRes.w(1:round(nf/2)-1,:,:)      = [];
StabRes.p(1:round(nf/2)-1,:,:)      = [];
StabRes.phi(1:round(nf/2)-1,:,:)    = [];

if Opt.Sweep
% Delete conjugates from results
StabRes.Asweep(1:round(nf/2)-1,:,:)        = [];
StabRes.alphasweep(1:round(nf/2)-1,:,:)    = [];
StabRes.usweep(1:round(nf/2)-1,:,:,:)      = [];
StabRes.vsweep(1:round(nf/2)-1,:,:,:)      = [];
StabRes.wsweep(1:round(nf/2)-1,:,:,:)      = [];
StabRes.psweep(1:round(nf/2)-1,:,:,:)      = [];
StabRes.phisweep(1:round(nf/2)-1,:,:,:)    = [];
end
 
end