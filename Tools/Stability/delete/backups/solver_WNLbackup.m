function [alpha,u,v,w,p,y,phi] = solver_WNL(i,Re,Ur,Vr,Wr,dxUr,dxVr,dxWr,dyUr,dyVr,dyWr,y_md,x_md,omega,beta,N,y_max,y_i,alpha,alpha1)
% Executes Weakly Nonparallel Local stability simulation from min(x_md) to x_md(:,M), at each given x-
% station of the given mean flow (Ur,Vr,Wr), assuming it to be given on a 
% uniform grid. The fields are to be given such as to go from the free-
% stream to the wall as the (row-)index increases. Downstream stations
% correspond to an increasing column index. 

[~,nx] =size(x_md);
%% Spectral basis functions

% Set up Chebyshev grid and associated Pseudo-Spectral Differentiation
% Matrices (PSDMs) in the computational domain: i.e. eta in [-1,1]
[eta,DM] = chebdif(N,2);

% Transform the PSDMs to the physical domain: eta |-> y
[y,D1,D2] = MappingMalik(y_max,y_i,eta',DM(:,:,1),DM(:,:,2));
clear DM eta

% The PSDMs will generally pop up with the size (N-2)x(N-2), accordingly 
% define the identity and zero column and row
I      = eye(N-2);
Zc     = zeros(N-2,1);
Zr     = zeros(1,N-2);

%% Initial Condition Generation with ILST (expanded Orr-Sommerfeld)

% Execute ILST analysis for initialization
tic
fprintf(' Running ILST for inflow conditions...\n')
[EVOS,uOS,vOS,wOS,pOS,yOS] = solver_ILST(...
    Re,Ur,Wr,dyUr,dyWr,y_md,...
    omega,beta,N,y_max,y_i);
fprintf('\b\b\b\b, finished in %.2f seconds.\n',toc)

% Find boundary layer thickness for filtering
id = find(Ur(:,1)<0.999,1); % index where U starts decreasing
d99 = interp1(Ur(id:end,1),y_md(id:end),0.99,'spline');

% Filter eigenfunctions for correct eigenvalue
[OSeigval,index] = EVfilter(EVOS,vOS,yOS,omega,beta,d99);%Changed inputs wrt OG

% Initialize solution vector phi = [u v w p]'(x,y) and alpha = alpha(x)
phi    = zeros(4*(N-2)+2,nx);
alpha  = zeros(1,nx);

% Set ILST eigenfunction in first phi station
phi(:,1) = ([uOS(2:end-1,index);
            vOS(2:end-1,index);
            wOS(2:end-1,index);
            pOS(:,index)]);
        
phi(:,1) = phi(:,1)/max(abs(uOS(:,index))); %Normalized phi
%max (u,v,w,p) ~= 1 now
        
% Set ILST eigenvalue in first alpha station
alpha(1) = OSeigval;

%% 
iu=sqrt(-1);
 
    % Interpolate base flow components and create diagonal matrix form
      U = diag(interp1(y_md,  Ur(:,1),y(2:N-1),'spline'));
    dxU = diag(interp1(y_md,dxUr(:,1),y(2:N-1),'spline'));
    dyU = diag(interp1(y_md,dyUr(:,1),y(2:N-1),'spline'));
    
      V = diag(interp1(y_md,  Vr(:,1),y(2:N-1),'spline'));
    dxV = diag(interp1(y_md,dxVr(:,1),y(2:N-1),'spline')); % note that this term is order 1/Re^2
    dyV = diag(interp1(y_md,dyVr(:,1),y(2:N-1),'spline'));
    
      W = diag(interp1(y_md,  Wr(:,1),y(2:N-1),'spline'));
    dxW = diag(interp1(y_md,dxWr(:,1),y(2:N-1),'spline'));
    dyW = diag(interp1(y_md,dyWr(:,1),y(2:N-1),'spline'));


    
    % Reset convergence measure
    dal = 1;    
    while abs(dal)*Re >= 1e-7
        
        % Compile common convection-diffusion terms (parallel terms only)
        Del = -iu*omega*I + iu*alpha(i)*U + iu*beta*W + I*(alpha(i)^2 + beta^2)/Re+V*D1(2:end-1,2:end-1)-D2(2:end-1,2:end-1)/Re;
        
%        Define L0 in system of ODE'S; parallel part of L * q
        L0 = [   Del               dyU          0*I [Zc iu*alpha(i)*I Zc]
                  Zr       -D2(1,  2:end-1)/Re   Zr       D1(1,  :)
                 0*I           Del+dyV          0*I     D1(2:end-1,:)
                  Zr       -D2(end,2:end-1)/Re   Zr       D1(end,:)
                 0*I              dyW           Del   [Zc iu*beta*I Zc]
             I*alpha(i)*iu D1(2:end-1,2:end-1) I*beta*iu [Zc 0*I Zc]];
         
%       Define L1 in system of ODE'S; nonparallel part, d/dx terms of L * q
        L1 = [   dxU          0*I         0*I       [Zc 0*I Zc]
                  Zr          Zr          Zr        [0   Zr  0]
                 0*I          0*I         0*I       [Zc 0*I Zc]
                  Zr          Zr          Zr        [0   Zr  0]
                 dxW          0*I         0*I       [Zc 0*I Zc]
                 0*I          0*I         0*I       [Zc 0*I Zc]];

%       Define L2 in system of ODE'S; M * dq/dx
        L2 = [U-I*(2*alpha(i)*iu)/Re   I*0      I*0                 [Zc  I  Zc]
                  Zr                   Zr       Zr                  [0   Zr  0]
                 I*0  U-I*(2*alpha(i)*iu)/Re    I*0                 [Zc I*0 Zc] 
                  Zr                   Zr       Zr                  [0   Zr  0]
                 I*0                   I*0  U-I*(2*alpha(i)*iu)/Re  [Zc I*0 Zc]
                  I                    I*0      I*0                 [Zc I*0 Zc]];
              
%        Define L3 in system of ODE'S; N * da/dx q
        L3  = [-I*(iu/Re)      I*0      I*0      [Zc I*0 Zc] 
                  Zr           Zr       Zr       [0   Zr  0]
                 I*0        -I*(iu/Re)  I*0      [Zc I*0 Zc] 
                  Zr           Zr       Zr       [0   Zr  0]
                 I*0           I*0    -I*(iu/Re) [Zc I*0 Zc] 
                 I*0           I*0      I*0      [Zc I*0 Zc]];
             
%        Define L4 in system of ODE'S; Base flow expansion terms from L0
%        and L1
        L4 = [iu*alpha(i)*dxU+iu*beta*dxW   D1(2:end-1,2:end-1)*dxU     I*0 [Zc I*0 Zc]
             Zr                      Zr                                 Zr  [0  Zr   0]
             I*0                     iu*beta*dxW                        I*0 [Zc I*0 Zc] 
             Zr                      Zr                                 Zr  [0  Zr   0]
             I*0                     dxW      iu*alpha(i)*dxU+iu*beta*dxW   [Zc I*0 Zc]
             I*0                     I*0                                I*0 [Zc I*0 Zc]];

        % Set up system of equations
        A = [L0 + L1 + alpha1(i)*L3        L2 ;             %Herbert system
             L4 + alpha1(i)*L2             L0];
         
        A = [L0 + L1 + alpha1(i)*L3             L2 ;        %Derived system
             zeros(size(L0)) L0+L1+L4+alpha1(i)*L3];        
        [q, D] = eig(A); % WRONG!!!!

         
% Extract the unique eigenvector part
        q0 = q(1:end/2,:);
        q1 = q(end/2+1:end,:);
        
        s = 2*(4*(N-2)+2);
        v1  = [zeros(1,s);  q0(1*(N-2)+1:2*(N-2)  ,:); zeros(1,s)]; 
        %Correct eigenvalues and vectors must be selected
        [eigv,index2] = EVfilter(diag(D),v1,y_md,omega,beta,d99);
        figure (1)
        plot(real(20*v1(:,index2)),y_md, '-')
        hold on
        plot(real(vOS(:,index)),y_md, '.')
        
        
% Augment Dirichlet conditions to perturbation velocity components
phi  = [0;  q0(0*(N-2)+1:1*(N-2)  ,index2); 0;
        0;  q0(1*(N-2)+1:2*(N-2)  ,index2); 0;
        0;  q0(2*(N-2)+1:3*(N-2)  ,index2); 0;
            q0(3*(N-2)+1:4*(N-2)+2,index2)];
                          
% Augment Dirichlet conditions to perturbation velocity components
phix = [0; q1(0*(N-2)+1:1*(N-2)  ,index2); 0;
        0; q1(1*(N-2)+1:2*(N-2)  ,index2); 0;
        0; q1(2*(N-2)+1:3*(N-2)  ,index2); 0;
           q1(3*(N-2)+1:4*(N-2)+2,index2)];
                 
% Extract the individual eigenvectors
u(:,i)   = phi(0*N+1:1*N);
v(:,i)   = phi(1*N+1:2*N);
w(:,i)   = phi(2*N+1:3*N);
p(:,i)   = phi(3*N+1:4*N);

ampltd = max(abs(u));

% Alpha calculation eq (42) in Agard793 on Weakly nonparallel LST

dal = 2*iu/(x_md(i)-x_md(i-1)) *(trapz(y,conj(u(:,i)).*(u(:,i)-u(:,i-1)))...
                               /(trapz(y,conj(u(:,i) .* u(:,i)))))
alpha(i) = alpha(i)-dal;





% Calculate the residual growth in the eigenfunction, DO OUTSIDE OF FUNC!
%alphaphi = zeros(1,nx);

%alphaphi(:,i) = -sqrt(-1)*trapz(repmat(y,4,1),phix(:,i).*phi(:,i))/trapz(repmat(y,4,1),abs(phi(:,i)).^2);


% % Report whether stabilization was activated during any of the steps
% if stabon
%     fprintf('\n Stabilization was activated.\n')% at station i=%.0f (%.2f%%).\n',stabon,100*stabon/M);
% end 


end
