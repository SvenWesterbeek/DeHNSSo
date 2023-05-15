function [mode] = ModeToModeNumber(n,m,N,M)
% Calculates the mode number for a given (n,m) out of a maximum of (N,M)
% M,m = omega modes
% N,n = beta modes

%% Error check
if n>N || m>M
    error('Input exceeds truncation (m,n) > (M,N)')
end

%% mint
L = (2*N+1)*(2*M+1); %Counter

%initialize basic building blocks (m,n)
M1 = -M:1:M; 
N1 = -N.*ones(1,2*N+1); 

% Set up omega vector
for i = 1:L 
 Mvec(i) = M1(mod(i-1,length(M1))+1);
end

% Set up beta vector
for k = 1:L/(2*M+1)
 Nvec([(k-1)*(2*M+1)+1:k*(2*M+1)]) = -N+(k-1)*ones(1,2*M+1);
end

% Find mode number
[~,mindex] = find(Mvec==m);
[~,nindex] = find(Nvec==n);
mode = intersect(mindex,nindex);


end

