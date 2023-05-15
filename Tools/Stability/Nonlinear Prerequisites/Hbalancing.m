function [ HB ] = Hbalancing(M, N, Mmat, Nmat, Modevec )
% Hbalancing creates a 3D matrix that will be used to calculate
% the source terms for the ILPSE solver. In each matrix HB(:,:,j), a 1 
% means that the modes corresponding to the row and column index of that 
% entry interact to force that mode j. For example, for N = 2, M = 0, 
% HB(:,:,3) corresponding to the MFD will look like:

%    n,m -2,0   -1,0    0,0     1,0     2,0        
%  n,m    __________________________________
% -2,0   | 0      0      0       0       1
% -1,0   | 0      0      0       1       0
%  0,0   | 0      0      1       0       0
% -1,0   | 0      1      0       0       0
% -2,0   | 1      0      0       0       0

% As the interactions (-2,0)x(2,0), (-1,0)x(1,0), (0,0)x(0,0), (1,0)x(-1,0
% and (2,0)x(-2,0) all force the MFD (0,0). In practice, the interaction
% (0,0)x(0x0) is excluded from forcing term calculations however.
%%  Inputs
%   M is the number of omega modes
%   N is the number of beta modes
%   Mmat contains the m part of the modes affected by the mode interactions
%   Nmat contains the n part of the modes affected by the mode interactions
%   Modevec contains the m and n part of the modecounter l

%%  Outputs
%   HB is a collection of [2M+1]x[2N+1] matrices of size [2M+1]x[2N+1]
%   stored as a 3D matrix where every matrix is a binary matrix that serves
%   as a registry for the source amplitudes that affect that respective mode.
%   (Harmonic Balancing)

L = (2*N+1)*(2*M+1);
for i = 1:L
% Read desired combination of m,n for mode i (mode l)
m=Modevec(1,i); 
n=Modevec(2,i);
A = Mmat == m;      %A tells us which entry has the correct mode for m
B = Nmat == n;      %B tells us which entry has the correct mode for n
HB(:,:,i) = A.*B;   %Simulates LOGIC AND, given that both are required.
end 

end

