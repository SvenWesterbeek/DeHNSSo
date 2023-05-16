function [BL,Y] = chebint_uniform(BL,y_max,y_i)
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  A MATLAB differentiation matrix suite
%              Copyright (C) 1998 J.A.C. Weideman, S.C. Reddy
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

%  Published as part of DeHNSSo with the written permission of J.A.C. Weideman.

%  J.A.C. Weideman, S.C. Reddy (1998)  A MATLAB differentiation matrix suite
%  Help notes modified by JACW, May 2003.

%%
% Interpolate the boundary solution matrices on a constant mapping defined
% by the length-scale "l" and with parameters "y_max" and "y_i".
%
% Mathematics used for interpolation are extracted from the 'chebint.m'
% function, which accelerates computations because the transformation
% changes per station while several matrices will remain constant.

fprintf(' Interpolating spectral BL on uniform grid...\n')
tic

% Get original grid dimensions
[ny,nx] = size(BL.u);
H = max(BL.Y);

% Define uniform grid
cheb = cos(linspace(0,pi,ny))'; % cosine grid on domain [1 0 -1]

% Interpolation matrix
chebM = repmat(cheb',ny,1);

% Weights for Chebyshev formula
W = ones(ny,1).*(-1).^(0:ny-1)';
W(1) = W(1)/2; 
W(ny) = W(ny)/2;

% Malik transformation:
eta = y_max*y_i*(1+cheb) ./ (y_max - cheb*(y_max-2*y_i)); % cosine grid on domain [ymax yi 0]

% Eta will now be constant with x, but is stored as matrix for
% convenience in the rest of the program.
BL.eta = repmat(eta,1,nx);

% Interpolate on a unique grid distribution per station
for i = 1:nx
    % Get local distribution from inversed Malik transformation 
    y_max_local = H/BL.l(i);
    y_i_local = BL.y_i/BL.l(i);

    % Mapping for interpolation function (inversed Malik transformation)
    eta_local = (y_max_local*eta-y_i_local*y_max_local) ./ ...
        (y_max_local*eta - 2*y_i_local*eta + y_i_local*y_max_local); % domain [1 0> -1>]

    % Interpolation matrix of new, scaled mapping
    etaM = repmat(eta_local,1,ny);

    % Interpolation
    D = etaM - chebM; % Compute quantities x-x(k)
    D = 1./(D+eps*(D==0));  % and their reciprocals.

    ChebInt = @(u) D*(W.*u)./(D*W); % Evaluate interpolant as matrix-vector products.

    % Perform interpolation
    BL.u(:,i) = ChebInt(BL.u(:,i));
    BL.w(:,i) = ChebInt(BL.w(:,i));
    
    if isfield(BL,'T') % Compressible BL has additional fields
        BL.T(:,i)   = ChebInt(BL.T  (:,i));
        BL.rho(:,i) = ChebInt(BL.rho(:,i));
        BL.mu(:,i)  = ChebInt(BL.mu (:,i));
        BL.k(:,i)   = ChebInt(BL.k  (:,i));
    end
end

if nargout==2
    Y = repmat(eta,1,nx).*repmat(BL.l,ny,1);
end

fprintf('\b\b\b\b, finished in %.2f seconds.\n',toc)
end