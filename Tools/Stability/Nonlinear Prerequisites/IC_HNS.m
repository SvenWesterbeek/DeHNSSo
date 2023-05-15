function [StabRes] = IC_HNS(j,i,Re,Ur,Wr,dyUr,dyWr,D1,D2,StabRes,StabGrid)
%% License

%% Description
%   The IC_HNS code calls ILST (that finds a solution to the local 
%   Orr-Sommerfeld problem). The result is normalized and stored in
%   StabRes.

%   Note: In HNS, ILST should only be used to introduce modes at linear 
%   amplitudes at the inflow therefore no nonlinear introduction methods 
%   are included here. To initialize nonlinear shape functions, please use 
%   the Stab.IC = "LOAD" option.

%% Overview of inputs
% i streamwise station
% j mode to solve for
% Re Reynolds number (global)
% Ur Streamwise velocity profile of the base flow
% Wr Spanwise velocity profile of the base flow
% dyUr Wall-normal derivative of streamwise velocity profile of base flow
% dyWr Wall-normal derivative of spanwise velocity profile of base flow
% D1 First-order spectral derivative operator
% D2 Second-order spectral derivative operator
% StabRes.omegavec omega  per mode
% StabRes.betavec beta per mode

%% Overview of outputs
% StabRes.u streamwise perturbation velocity profile
% StabRes.v Wall-normal perturbation velocity profile
% StabRes.w Spanwise perturbation velocity profile
% StabRes.p Perturbation pressure profile
% StabRes.phi Perturbation vector of state variables (u,v,w,p)^T

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

% Visual Check for Mode Selection
% Eigenspectrum plot
figure
subplot(1,2,1)
plot(eig,'k.')
axis([-2 2 -2 2]) % make a function of the eig or alpha
hold on
plot(alpha,'ro')
subplot(1,2,2)
plot(abs( u(:,index)),StabGrid.etaun)

% Calculate new maximum amplitude after spline interpolation
u                 = u(:,index);
% Check if we can use ChebInt
y_mdInt           = linspace(StabGrid.etaun(1),StabGrid.etaun(end),4000);
u1int             = interp1(StabGrid.etaun,abs(u),y_mdInt,'spline');
[ampltd, ~]       = max(abs(u1int));
StabRes.u(j,:,i)  = u/ampltd;
StabRes.v(j,:,i)  = v(:,index)/ampltd;
StabRes.w(j,:,i)  = w(:,index)/ampltd;
StabRes.p(j,:,i)  = p(:,index)/ampltd;

% Create phi from eigenvectors
StabRes.phi(j,:,i) = ([ StabRes.u(j,:,i).';
                        StabRes.v(j,:,i).';
                        StabRes.w(j,:,i).';
                        StabRes.p(j,:,i).']);
end



