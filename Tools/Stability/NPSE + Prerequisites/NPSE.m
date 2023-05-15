function [StabRes,StabGrid] = NPSE(BF,Grid,Stab,NonDim,Opt)
% NPSE function: 
% Executes NPSE simulation from min(StabGrid.xiun) to StabGrid.xiun(:,nx), at each given x-
% station of the given mean flow (Ur,Vr,Wr). The fields are to be given 
% such as to go from the free-stream to the wall as the (row-)index 
% increases. Downstream stations correspond to an increasing column index. 

% NPSE Notes
% - Nonlinear terms are modelled as a source term.
% - Shapefunctions are not adjusted during the loop to have u_max = 1, rather
%   the Amplitude is saved and adjusted in the end to increase numerical
%   stability.
% - Since this concerns an NPSE solver, omega and beta are now arrays whose 
%   combinations make up the modes.
% - This solver will fail when trying to run an NPSE from a nonlinear
%   starting point. This will cause many modes to be introduced in one step
%   which are not yet converged. Initiate from linear region (upstream, lower amplitude).

% IF the solver fails try the following steps:
% Case: Solver fails early in the simulation after introducing one of the
% first harmonics
% Suggestion: Adjust H so the EV filter can pick the correct eigenfunction.

% IF the MFD pressure blows up, try increasing the division (suppression) of
% dp/dx in the systems for the MFDI (introduction)[path: NPSE -> IC -> MFDI]
% and apply an equal suppresion in the MFD system of the ILPSE. 

% IF the higher harmonics cause a failure to converge (Error at the NLT
% caller for alpha and dal calculation.) when using ILST forced mode
% introduction.
% Try to decrease H (not too much that earlier harmonics cannot be
% introduced correctly)
% Then try to have y_i closer to the wall and increase ny to remain a small
% dy for the earlier modes.

% Sven Westerbeek, MSc 
% Thesis "Development of a Nonlinear Parabolized Stability Equation (NPSE) 
% Solver for Transition Prediction in Boundary Layers" 2020

% Loop legend: i = marching (x), j = modes

%% Load default settings

% Define imaginary unit
iu=sqrt(-1);
warning('off', 'MATLAB:nearlySingularMatrix') % Required to suppress warning

if ~isfield(Grid,'ft'), Grid.ft = 1; end
if ~isfield(Grid,'ytype'), Grid.ytype ="malik"; end
if ~isfield(Grid,'xtype'), Grid.xtype ="equidistant"; end

if ~isfield(Opt,'TH'); Opt.TH = 1e-11; end
if ~isfield(Opt,'Conv'); Opt.Conv = 1e-8 ; end
if ~isfield(Opt,'UR'); Opt.UR = 0; end 

%% Setup mode vector and harmonic balancing matrices.

% Calculate Mode INTeraction matrix normalized by the fundamental mode
[Nmat, Mmat, Modevec,Mvec,Nvec] = Mintv2(Stab.M,Stab.N);

% Perform Harmonic Balancing for each mode
[HB] = Hbalancing(Stab.M, Stab.N, Mmat, Nmat, Modevec ); 
%Rvec = [Nvec;Mvec];

% Create vectors of omegas and betas per mode
StabRes.omegavec = Stab.omega_0*Mvec;
StabRes.betavec = Stab.beta_0*Nvec;

L = 1:(2*Stab.N+1)*(2*Stab.M+1); % Mode counter
nf = max(L); % Number of modes

%% Grid Generation
fprintf('Generating grid. \n')
[StabGrid,D1,D2]=grid_gen(Grid);

%% Base Flow and BC Interpolation on numerical grid

% Interpolate base flow velocities onto numerical grid
fprintf('Interpolating base flow on the numerical grid. \n')
Ur = fillmissing(griddata(BF.X,BF.Y,BF.U,StabGrid.x,StabGrid.y,'cubic'),'pchip');
Vr = fillmissing(griddata(BF.X,BF.Y,BF.V,StabGrid.x,StabGrid.y,'cubic'),'pchip');
Wr = fillmissing(griddata(BF.X,BF.Y,BF.W,StabGrid.x,StabGrid.y,'cubic'),'pchip');

% Enforce no-slip condition
Ur(end,:)=0; Vr(end,:)=0; Wr(end,:)=0;

% Interpolate derivatives of the base flow velocities onto numerical grid
% Filling missing values for robustness, this can happen when the domain
% description does not exactly match the base flow domain 
dxUr = fillmissing(griddata(BF.X,BF.Y,BF.dxU,StabGrid.x,StabGrid.y,'cubic'),'pchip');
dxVr = fillmissing(griddata(BF.X,BF.Y,BF.dxV,StabGrid.x,StabGrid.y,'cubic'),'pchip');
dxWr = fillmissing(griddata(BF.X,BF.Y,BF.dxW,StabGrid.x,StabGrid.y,'cubic'),'pchip');

dyUr = fillmissing(griddata(BF.X,BF.Y,BF.dyU,StabGrid.x,StabGrid.y,'cubic'),'pchip');
dyVr = fillmissing(griddata(BF.X,BF.Y,BF.dyV,StabGrid.x,StabGrid.y,'cubic'),'pchip');
dyWr = fillmissing(griddata(BF.X,BF.Y,BF.dyW,StabGrid.x,StabGrid.y,'cubic'),'pchip');

%% Initial Condition Generation 
% Execute ILST analysis to initialize marching loop

% Find boundary layer thickness for filtering
id = find(Ur(:,1)/max(Ur(:,1))<0.999,1); % index where U starts decreasing
d99 = interp1(Ur(id:end,1)/max(Ur(id:end,1)),StabGrid.etaun(id:end),0.99,'spline');

% Initialize solution vector phi = [u v w p]'(x,y), amplitude and source f
StabRes.A       = zeros(size(L));
StabRes.phi     = zeros(nf,(4*(Grid.ny-2)+3),Grid.nx); 
f               = zeros(nf,(4*(Grid.ny-2)+3),Grid.nx);
StabRes.alpha   = zeros(nf,Grid.nx);

% Initializing u, v, w, p matrices with (size that includes BC)
StabRes.u       = zeros(nf,Grid.ny,Grid.nx);
StabRes.v       = zeros(nf,Grid.ny,Grid.nx);
StabRes.w       = zeros(nf,Grid.ny,Grid.nx);
StabRes.p       = zeros(nf,Grid.ny,Grid.nx);

[~,nonzeromodes]= find(Stab.A0>0); % Determines which modes require IC
RunJ = nonzeromodes; % RunJ vector sets which modes to run for


for j = flip(nonzeromodes) %Run initial condition for all nonzero modes 
if j >= round(nf/2)% Mode is physical
    [StabRes.u(j,:,1),StabRes.v(j,:,1),StabRes.w(j,:,1),StabRes.p(j,:,1),StabRes.alpha(j,1),StabRes.phi(j,:,1)]=...
    IC(j,1,NonDim.Re,D1,D2,Ur,Vr,Wr,dxUr,dxVr,dxWr,dyUr,dyVr,dyWr,d99,f,HB,RunJ,StabGrid,Grid,StabRes);       

else % Conjugate is mirror of physical
    j2 = round(nf/2)-(j-round(nf/2));% Defines conjugate mode
    
    StabRes.u(j,:,1)   = conj(StabRes.u(j2,:,1));
    StabRes.v(j,:,1)   = conj(StabRes.v(j2,:,1));
    StabRes.w(j,:,1)   = conj(StabRes.w(j2,:,1));
    StabRes.p(j,:,1)   = conj(StabRes.p(j2,:,1));
    StabRes.alpha(j,1) = -conj(StabRes.alpha(j2,1)); % growth rates should be equal between conjugate modes so -conj
    StabRes.phi(j,:,1) = conj(StabRes.phi(j2,:,1));
end                                                                            
end 

%% Initialize Nonlinear convergence loop
%Numerical initializations
count       = zeros(1,Grid.nx);  % Debugger value
skip        = 0;            % Will be used to wait one stage after mode has turned significant
NewModesRun = [];
RunSource   = []; %Initializes the source term one stage before mode initiation
Introi      = zeros(1,nf);
iref        = 0; 
Aactive     = zeros(size(StabRes.A));
Aactive(nonzeromodes) = 1;
A0save      = 2*Stab.A0; 
TH          = Opt.TH; %
%% Marching algorithm 
for i = 2:Grid.nx
    %% Extrapolate results to current stage as initial guess

    % Estimate shape functions and alpha to be the same initially
for j = 1:nf
    StabRes.u(j,:,i)   = StabRes.u(j,:,i-1);
    StabRes.v(j,:,i)   = StabRes.v(j,:,i-1);
    StabRes.w(j,:,i)   = StabRes.w(j,:,i-1);
    StabRes.p(j,:,i)   = StabRes.p(j,:,i-1);
    StabRes.alpha(:,i) = StabRes.alpha(:,i-1);   
end

%% Converging loop for phi, f and alpha 
 
 dal = ones(nf,1); % Convergence measure
 Conv = Opt.Conv;  % Convergence criterium
 while max(abs(dal)) >= Conv 
 count(i) = count(i)+1; %Counter, Debugging help
        
%% Initialization
dal = zeros(nf,1); %Reset convergence measure
tic

% Update Amplitude at stage i
for j = RunJ
  Nfac=trapz(StabGrid.xiun(1:i),-imag(StabRes.alpha(j,1:i)));
  StabRes.A(j)=Stab.A0(j)*exp(Nfac);
end

% Amplitudes of RHS mode interactions
Amat = StabRes.A'*StabRes.A; 

%% Check RHS A of EXPECTED active modes

mac = RunJ - round(nf/2);
RNmac = 0;
RMmac = 0;
for ii = 1:length(mac)
 Nmac = Nvec(RunJ); %N of active modes
 Mmac = Mvec(RunJ); %m of active modes
 RNmac = [RNmac Nmac(ii)+Nmac];
 RMmac = [RMmac Mmac(ii)+Mmac];
end
Rmac = unique([RNmac; RMmac]', 'rows')' ; %Mode combinations under consideration

% Find which mode numbers are given in Rmac
for ii = 1:length(Rmac)
    if ~isempty(ModeToModeNumber_v2(Rmac(1,ii),Rmac(2,ii),Stab.N,Stab.M))
    rmac(ii) = ModeToModeNumber_v2(Rmac(1,ii),Rmac(2,ii),Stab.N,Stab.M);
    end
end

% Sort and if necessary delete modes out of spectral domain
rmac = unique(sort([rmac RunJ]));
rmac(rmac<1)=[];
rmac(rmac>nf)=[];

Amode = zeros(size(StabRes.A));
%Limit Mode interactions to one harmonic over highest harmonic (n_max+1)
if StabRes.A(round(nf/2)+1) == 0 %CF case
[~, nzA] = find(StabRes.A); %find nonzero Amplitudes
nfm = max(nzA)+(Stab.M+Stab.N+1);     %find index of highest harmonic+1
nfl = min(nzA)-(Stab.M+Stab.N+1);     %find index of conjugate of highest harmonic-1
for j=rmac
    if j > nfm || j<nfl
       Amode(j) = 0; 
    else 
       Amode(j)  = sum(Amat(find(HB(:,:,j)))); 
    end
end
else %TS case
    [~, nzA] = find(StabRes.A); %find nonzero Amplitudes
nfm = max(nzA)+(1);     %find index of highest harmonic+1
nfl = min(nzA)-(1);     %find index of conjugate of highest harmonic-1
for j=rmac
    if j > nfm || j<nfl
       Amode(j) = 0; 
    else
       Amode(j)  = sum(Amat(find(HB(:,:,j))));     
    end
end
end

%% 
Aactive_old = Aactive; 
if iref ~= i
% Check for new modes compared to previous stage only at the start of the
% stage
Aactive = (StabRes.A + Amode) >= TH; 
end
NewModes = find(Aactive - Aactive_old > 0);
iref = i;
%% Run IC for New Modes 
for j = flip(NewModesRun) 
    % Flips order to first run phyiscal modes so conjugates can be mirrored
    fprintf([' \n Introducing mode ', num2str(j-round(nf/2)) '.'])
    pause(0.05)
    Introi(j) = i; %
    if j >= round(nf/2)% Mode is physical
    [StabRes.u(j,:,i),StabRes.v(j,:,i),StabRes.w(j,:,i),StabRes.p(j,:,i),StabRes.alpha(j,i),StabRes.phi(j,:,i),StabRes.A(j)]=...
    IC(j,i,NonDim.Re,D1,D2,Ur,Vr,Wr,dxUr,dxVr,dxWr,dyUr,dyVr,dyWr,d99,f,HB,RunJ,StabGrid,Grid,StabRes);  
    
    Stab.A0(j) = StabRes.A(j); %Set amplitude at introduction A0
    
    % For MFD set A = 1 as amplitude is maintained in the shape function
    if j == round(nf/2)
        Stab.A0(j) = 1;
        StabRes.A(j) = 1;
    end

    else %Mode is conjugate of physical mode, This forces an equal phase
    j2 = round(nf/2)-(j-round(nf/2));% Defines conjugate mode
    
    StabRes.u(j,:,i) = conj(StabRes.u(j2,:,i));
    StabRes.v(j,:,i) = conj(StabRes.v(j2,:,i));
    StabRes.w(j,:,i) = conj(StabRes.w(j2,:,i));
    StabRes.p(j,:,i) = conj(StabRes.p(j2,:,i));
    StabRes.alpha(j,i) = -conj(StabRes.alpha(j2,i)); % growth rates should be equal between conjugate modes so -conj
    StabRes.phi(j,:,i)  = conj(StabRes.phi(j2,:,i));
    StabRes.A(j)       = StabRes.A(j2);
    Stab.A0(j)      = Stab.A0(j2);
    end
end

NewModesRun = NewModes;
RunJ = nonzeros(Aactive .* L)';

% Delay mode introduction by one stage
if ~isempty(NewModesRun) || i == skip
    if i ~= skip
    Delete = []; % Reset Delete array for new entries
    [~,ct] = size(NewModesRun);
    for ii = 1:ct
    Delete(ii) = find(RunJ == NewModesRun(ii));
    end
    end
skip = i;
RunJ(Delete)=[];
RunSource = NewModesRun;
end
%% Source term update
if ge(count(i),2)
fold = squeeze(f(:,:,i));
end

f = NLT_NPSE(f,i,RunJ,RunSource,StabGrid,StabRes,Grid,Stab,D1,HB);

%Source term under relaxation following Zhao Lei et al.
if Opt.UR == 1
if ge(count(i),2)
dt = 0.2;
if i >= 179 %Only for the first stages
for j = RunJ
f(j,:,i) = ((1-dt).*fold(:,j)+dt.*f(j,:,i));
end
end
end
end

%% ILPSE Solver
% Calls the ILPSE solver which is the LPSE solver that is adjusted to include a source term.
% This solver does not march on its own since the inputs must be adjusted
% for source terms, growth and relevant modes every step.

for j =flip(RunJ)
if j >= round(nf/2)% Mode is physical
    
    [StabRes] = solver_ILPSE(i,j,NonDim.Re,Ur,Vr,Wr,...
                    dxUr,dxVr,dxWr,dyUr,dyVr,dyWr,...
                    f,D1,D2,StabRes,Grid,StabGrid);
        
else %Mode is conjugate of physical mode
    j2 = round(nf/2)-(j-round(nf/2));% Defines the (physical) conjugate mode

    StabRes.u(j,:,i) = conj(StabRes.u(j2,:,i));
    StabRes.v(j,:,i) = conj(StabRes.v(j2,:,i));
    StabRes.w(j,:,i) = conj(StabRes.w(j2,:,i));
    StabRes.p(j,:,i) = conj(StabRes.p(j2,:,i));
    StabRes.alpha(j,i) = -conj(StabRes.alpha(j2,i)); % growth rates should be equal between conjugate modes so -conj
    StabRes.phi(j,:,i) = conj(StabRes.phi(j2,:,i));

end
end

%% Calculate alpha
[StabRes,dal] = AlphaCalc(i,RunJ,Grid,StabGrid,StabRes);


    %% Convergence plotting
    if ~mod(i,10) %Only plot every 10 steps to reduce CPU load
%     plotdal(i) = max(abs(dal))*Re;
%     figure(1001)
%     set(gcf, 'Position',  [701, 1, 700, 400])
%     axis([0 Grid.nx 0 0.001])
%     hold on
%     plot(i,plotdal(i),'.r')
%     pause(0.05)
    
    %% Shapes plotting
    j = round(nf/2)+0*(2*Stab.N+1); %
    if ismember(j,RunJ)
    if j == round(nf/2)
        frame = [0 max(abs(StabRes.u(j,:,i))) 0 Grid.H];
    else
        frame = [0 1 0 y_max];
    end
    
    figure(1002)
    set(gcf, 'Position',  [1, 1, 700, 785])
    title(['Velocity and forcing term plots for omega = ' num2str(StabRes.omegavec(j))])
    subplot(2,4,1)
    plot(abs(StabRes.u(j,:,i)),StabGrid.etaun)
    if j == round(nf/2)
    title('MFD u')
    else
    title('u')
    end
    axis(frame)
    
    subplot(2,4,2)
    plot(abs(StabRes.v(j,:,i)),StabGrid.etaun)
    title('v')
    axis(frame)
    
    subplot(2,4,3)
    plot(abs(StabRes.w(j,:,i)),StabGrid.etaun)
    title('w')
    axis(frame)
    
    subplot(2,4,4)
    plot(abs(StabRes.p(j,:,i)),StabGrid.etaun)
    title('p')
    
    subplot(2,4,5)
    plot(abs(f(j,1:Grid.ny,i)),StabGrid.etaun)
    title('f_x')
    
    subplot(2,4,6)
    plot(abs(f(j,Grid.ny-2+1:2*(Grid.ny-2)+1,i)),StabGrid.etaun(1:end-1))
    title('f_y')
    
    subplot(2,4,7)
    plot(abs(f(j,2*(Grid.ny-2)+2:3*(Grid.ny-2)+1,i)),StabGrid.etaun(2:end-1))
    title('f_z')
    
    subplot(2,4,8)    
    title('f_{c}') %Always 0
    end
    end
    
end %End of while loop
    
    %Save amplitude of every mode at every stage
    A2(:,i) = 2*StabRes.A;
    mfdA(i) = max(abs(StabRes.u(round(nf/2),:,i)));
    
    figure(1003)
    set(gcf, 'Position',  [701, 400, 700, 450])
    title(['Amplitudes'])
    hold off
    semilogy(StabGrid.xiun(1:i),A2(RunJ(RunJ>round(nf/2)),1:i))
    hold on
    semilogy(StabGrid.xiun(1:i),A2(round(nf/2),1:i).*mfdA/2,'k.')
    
%     if i == 2
%     wb = waitbar(percentage/100, 'Please wait, the NPSE is running.'); 
%     set(gcf, 'Position',  [1, 1, 700, 785])
%     end
%     % Notify progress   
%     if round(100*i/Grid.nx) == percentage
%         %for j = 1:ceil(log10(percentage))+9; fprintf('\b'); end
%         %fprintf('%.0f percent ',percentage)
%         waitbar(percentage/100) 
%         percentage = percentage + 1;
%     end
    %A_MFD(i) = StabRes.A(round(nf/2));
end %End of marching loop

     %Normalize Shape functions and correct Amplitude
     %Also register MFD amplitude
     StabRes.A = A2; % Amplitude register 
     
     for i = 1:Grid.nx
     StabRes.A(round(nf/2),i) = max(abs(StabRes.phi(round(nf/2),1:Grid.ny,i)));
     for j = RunJ(RunJ>=round(nf/2))

     %    Calculate new maximum amplitude after spline interpolation
    y_mdInt             = linspace(StabGrid.etaun(1),StabGrid.etaun(end),4000);
    u1int               = interp1(StabGrid.etaun,abs(StabRes.u(j,:,i)),y_mdInt,'spline');
    [ampltd, ymax(1)]   = max(abs(u1int));

     
     %normalize shapes
     StabRes.u(j,:,i) = StabRes.u(j,:,i)./ampltd;
     StabRes.v(j,:,i) = StabRes.v(j,:,i)./ampltd;
     StabRes.w(j,:,i) = StabRes.w(j,:,i)./ampltd;
     StabRes.p(j,:,i) = StabRes.p(j,:,i)./ampltd;
     StabRes.phi(j,:,i) = StabRes.phi(j,:,i)./ampltd;
     
     %Correct Amplitude
     if j ~=round(nf/2)
     StabRes.A(j,i) = StabRes.A(j,i)*ampltd;
     j2 = nf-j+1;
     StabRes.A(j2,i) = StabRes.A(j,i)*ampltd;
     end
     end
     
     end
     
     % Use recorded amplitude 
       for j = 1:nf
           if j ~= round(nf/2)
       StabRes.A(j,1) = A0save(j);
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
       
end %End of function