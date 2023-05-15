function [ ynd,D1y,D2y ] = MappingMalik( y_max,y_i,eta,D1eta,D2eta )
% This function maps the Chebyshev grid on the domain [-1,0,1] to
% [0,y_i,y_max]. Furthermore, it delivers the pseudo-spectral 
% differentiation matrices that correspond to this transformation

% Transform Chebyshev to physical grid
ynd = (y_i*y_max*(1+eta)./(y_max - eta*(y_max - 2*y_i)))';

% Compute scale factors w.r.t. differentiation
detady   = (y_max - eta*(y_max - 2*y_i)).^2/(2*y_i*y_max*(y_max - y_i));
d2etady2 = -2*detady.^2*(y_max - 2*y_i)./(y_max - eta*(y_max - 2*y_i));

D1y = diag(detady)*D1eta;    
D2y = diag(d2etady2)*D1eta + diag(detady.^2)*D2eta;

end

