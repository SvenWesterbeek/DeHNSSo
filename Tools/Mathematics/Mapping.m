function [ ynd,D1y,D2y ] = Mapping(Grid,eta,D1eta,D2eta)


if Grid.ytype ~= "step"
%% Mapping Malik
% Transform Chebyshev to physical grid
ynd = (Grid.y_i*Grid.H*(1+eta)./(Grid.H - eta*(Grid.H - 2*Grid.y_i)))';

if nargout > 1
    % Compute scale factors w.r.t. differentiation
    detady   = (Grid.H - eta*(Grid.H - 2*Grid.y_i)).^2/(2*Grid.y_i*Grid.H*(Grid.H - Grid.y_i));
    d2etady2 = -2*detady.^2*(Grid.H - 2*Grid.y_i)./(Grid.H - eta*(Grid.H - 2*Grid.y_i));

    D1y = diag(detady)*D1eta;    
    D2y = diag(d2etady2)*D1eta + diag(detady.^2)*D2eta;
end

else
%% Mapping stretched
% This function maps the Chebyshev grid on the domain [-1,0,1] to
% [0,Grid.y_i,Grid.H]. Furthermore, it delivers the pseudo-spectral 
% differentiation matrices that correspond to this transformation

%Correction to Grid.y_i for stretching
Grid.y_i = 0.5*0.43655*Grid.y_i; 
n = length(eta);

% Transform Chebyshev to physical grid
%ynd = (Grid.y_i*Grid.H*(1+eta)./(Grid.H - eta*(Grid.H - 2*Grid.y_i)))';
%ynd = [ynd(1:floor(length(eta)/2)); ...
%        linspace(ynd(floor(length(eta)/2)+1),0,floor(length(eta)/2))'];


% Set Distribution of the first Constant part
dy1 = ones(1,floor(n/2)).*Grid.y_i/floor(n/2);

% Calculate the distance that is covered by the quadratic function
A3 = Grid.H-dy1(1)*n;

syms nn
alpha = A3/int((nn-floor(n/2))^2,floor(n/2),n);
alpha = double(alpha);

for j = 1:ceil(n/2)
    if j == 1
    dy2(j) = dy1(1);
    else
    dy2(j) = alpha*(j)^2+dy1(1);
    end
end

%ynd = flipud(cumtrapz([dy1 dy2])');

% Stretched original (works OKAY)
k = 2;
yndStretch = k*(Grid.y_i*Grid.H*(1+eta)./(Grid.H - eta*(Grid.H - 2*k*Grid.y_i)))';

% Stretched original (works OKAY)
k = 1;
yndOG = k*(Grid.y_i*Grid.H*(1+eta)./(Grid.H - eta*(Grid.H - 2*k*Grid.y_i)))';

% Stretching test
f= Grid.ystretch;

ynd = (1+eta')./(f+eta')*(f+1)/2*Grid.H;
ynd = ynd-ynd(end);
ynd = ynd/max(ynd)*Grid.H;


ynd = (ynd + 20*yndOG + yndStretch)/22;

figure(5)
hold on
plot(ynd,'.')

% Linear Distribution (equidistant) (For TESTING)
%ynd = flipud(linspace(0,Grid.H,n)');

% eta-based test distribution
%ynd = 1/sqrt(2)*((eta+1).^0.5)';

%% Calculate elements of the derivative matrix E
%Set up c (eq 3.22)
c = ones(1,n);
c(1) = 2;
c(end) = 2;


for j = 1:n
    for k = 1:n
       if k~=j
           E(j,k) = c(j)/c(k) * (-1)^(k+j) / (eta(j)-eta(k));
       else
           E(j,j) = - (eta(j) / (2*(1-eta(j)^2)));
       end
    end
end

E(1,1) = (2*n^2+1)/6; %
E(n,n) = -E(1,1);

% Calculate second derivative matrix
E2 = zeros(size(E));
for j = 1:n
    for k = 1:n
        for m = 1:n
E2(j,k) = E2(j,k) + E(j,m)*E(m,k);
        end
    end
end

%% Set up D1y and D2y
% Calculate the scaling factor for transformation physical -> comp
S = FD1d2o_uneven(ynd,eta')';

S2 = FD2d2o_uneven(ynd,eta')';


% Calculate first derivative matrix
% D1y = diag(detady)*D1eta;    
D1y = diag(S)*D1eta;
% D2y = diag(d2etady2)*D1eta + diag(detady.^2)*D2eta;
D2y = diag(S2)*D1eta + diag(S.^2)*D2eta;

end % grid mode
end % of function


