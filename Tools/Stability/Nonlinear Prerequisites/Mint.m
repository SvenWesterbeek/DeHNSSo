function [Nmat, Mmat, Modevec,Mvec, Nvec] = Mint(M,N)

%% License

%% Description
% mint (Mode INTeractions) is a script that sets up a matrix containing the
% modes affected by the interactions of all currently represented modes

% Outputs Mmat and Nmat are mode numbers for beta and omega respectively.
% Output Modevec is the mode vector [m(1:2M+1);n(1:2M+1)] for every mode
% counter l_{m,n}.


%% mint
L = (2*N+1)*(2*M+1); %Counter

%initialize basic building block
M1 = -M:1:M; 

% Set up omega vector
for i = 1:L 
 Mvec(i) = M1(mod(i-1,length(M1))+1);
end

% Set up beta vector
for k = 1:L/(2*M+1)
 Nvec((k-1)*(2*M+1)+1:k*(2*M+1)) = -N+(k-1)*ones(1,2*M+1);
end
Modevec = [Mvec;Nvec];

% Find recipient mode for interaction of mode i with mode j
for i = 1:L %Multiplying powers is addition in power
    for j = 1:L
        Nmat(i,j) = Nvec(i)+Nvec(j); 
    end
end
for i = 1:L %Multiplying powers is addition in power
    for j = 1:L
        Mmat(i,j) = Mvec(i)+Mvec(j);
    end
end
end

