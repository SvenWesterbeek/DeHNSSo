function [EV,u,v,w,p,fn]=solver_ILST(Re,Ur,Wr,dUr,dWr,omega,beta,N,D1,D2)
% This function solves the Orr-Sommerfeld problem on domain y in [0,y_i,
% y_max], mapping half the nodes between y = 0 and y_i, using the same sca-
% ling as y_md. Clamped conditions are imposed by Dirichlet conditions at 
% the boundary nodes and using "adjusted" Lagrange basis (through cheb4c) 
% that identically satisfies the Neumann condition implicitly. The base 
% flow parameters are interpolated on the Chebyshev grid using the interpo-
% lation method for extrapolation when y_max > max(y_md).
%
% The most unstable mode, present in the search box, is selected from the
% 2*(4*N-6) eigenpairs, value eigval(i) and functions u, v, w and p, and 
% returned, where the eigenfunctions are defined on the N nodes stored in 
% y. fn times machine precision is the eigensolver error

% The transformed matrices represent the full operators. Removing the outer 
% elements yields the Dirichlet conditions
D1e = D1(1,  :);
D1w = D1(end,:);
D1x = D1(2:N-1,:);
D1  = D1(2:N-1,2:N-1);
D2e = D2(1,  :);
D2w = D2(end,:);
D2  = D2(2:N-1,2:N-1);
% Note that these matrices now have the size (N-2)x(N-2), accordingly de-
% fine the identity and zero matrix
I   = eye(N-2);
Z   = zeros(N-2);
Zc  = zeros(N-2,1);
Zr  = zeros(1,N-2);

% Transform Base Flow to matrix form
U  = diag( Ur(2:N-1));
DU = diag( dUr(2:N-1));
W  = diag( Wr(2:N-1));
DW = diag( dWr(2:N-1));

% Imaginary unit
iu=sqrt(-1);

% Diagonal convection-diffusion terms
CD = (beta^2*I - D2)/Re + iu*beta*W - iu*omega*I;

%% Define A in A*Xi = alpha*B*Xi + alpha^2*C*Xi

% Define the eigenvalue problem
A11 = CD;
A12 = DU;
%A13
%A14

%A21
A22 = CD;
PC2u2 = -D2e(2:end-1)/Re;
PC2l2 = -D2w(2:end-1)/Re;
%A23
A24 = D1x;
PC2u4 = D1e(1,:);
PC2l4 = D1w(end,:);

%A31
A32 = DW;
A33 = CD;
A34 = iu*beta*I;
clear CD

%A41
A42 = D1;
A43 = iu*beta*I;
%A43

%% Define B in A*Xi = alpha*B*Xi + alpha^2*C*Xi

B11 = -iu*U;
%B12
%B13
B14 = -iu*I;

%B21 
B22 = -iu*U;
%B23
%B24

%B31
%B32
B33 = -iu*U;
%B34

B41 = -iu*I;
%B42
%B43
%B44

%% Define C in A*Xi = alpha*B*Xi + alpha^2*C*Xi

C11 = -I/Re;
% C12 C13 C14

C22 = -I/Re;
% C22 C23 C24

C33 = -I/Re;
% C32 C33 C34

% C41 C42 C43 C44

%% Generate matrices

% LHS matrix
%     u    v    w   pe  p  pw
A = [A11  A12   Z  [Zc  Z  Zc]   % x-momentum
      Zr PC2u2  Zr    PC2u4      % y-compatibility (edge boundary)
      Z   A22   Z      A24       % y-momentum
      Zr PC2l2  Zr    PC2l4      % y-compatibility (wall boundary)
      Z   A32  A33 [Zc A34 Zc]   % z-momentum
      Z   A42  A43 [Zc  Z  Zc]]; % continuity

% RHS matrix (alpha)
B = [B11   Z    Z  [Zc B14 Zc] 
      Zr   Zr   Zr [0   Zr  0]
      Z   B22   Z  [Zc  Z  Zc] 
      Zr   Zr   Zr [0   Zr  0]
      Z    Z   B33 [Zc  Z  Zc] 
     B41   Z    Z  [Zc  Z  Zc] ];
 
% RHS matrix (alpha^2)
C = [C11   Z    Z  [Zc  Z  Zc] 
      Zr   Zr   Zr [0   Zr  0]
      Z   C22   Z  [Zc  Z  Zc] 
      Zr   Zr   Zr [0   Zr  0]
      Z    Z   C33 [Zc  Z  Zc] 
      Z    Z    Z  [Zc  Z  Zc] ];
  
%% Construct companion matrices 
% define A-sized zero and identity matrices 
Z = zeros(size(A));
I =   eye(size(A));
% build companian matrices
A = [A -B; Z I];
B = [Z  C; I Z];

% Determine Frobenius norm of A; error of eig = machine precision * ||A||_F
fn = norm(A,'fro');
%% Solve the eigenvalue problem
%tic
[V,D] = eig(A,B);
%toc

% Extract the eigenvalues
EV = diag(D);

% Extract the unique eigenvector part
V = V(1:end/2,:);

% Augment the Dirichlet boundary conditions to the functions
V = [zeros(1,2*(4*(N-2)+2)); V(0*(N-2)+1:1*(N-2)+0,:); zeros(1,2*(4*(N-2)+2));...
     zeros(1,2*(4*(N-2)+2)); V(1*(N-2)+1:2*(N-2)+0,:); zeros(1,2*(4*(N-2)+2));...
     zeros(1,2*(4*(N-2)+2)); V(2*(N-2)+1:3*(N-2)+0,:); zeros(1,2*(4*(N-2)+2));...
                             V(3*(N-2)+1:4*(N-2)+2,:)]; 

% Extract the eigenvectors
u   = V(0*N+1:1*N,:);
v   = V(1*N+1:2*N,:);
w   = V(2*N+1:3*N,:);
p   = V(3*N+1:4*N,:);

end

