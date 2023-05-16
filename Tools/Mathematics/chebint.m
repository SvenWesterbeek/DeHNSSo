function p = chebint(fk, x)
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
%% Description
%  The function p = chebint(fk, x) computes the polynomial interpolant
%  of the data (xk, fk), where xk are the Chebyshev nodes.  
%  Two or more data points are assumed.
%
%  Input:
%  fk:  Vector of y-coordinates of data, at Chebyshev points 
%       x(k) = cos((k-1)*pi/(N-1)), k = 1...N.
%  x:   Vector of x-values where polynomial interpolant is to be evaluated.
%
%  Output:
%  p:    Vector of interpolated values.
%
%  The code implements the barycentric formula; see page 252 in
%  P. Henrici, Essentials of Numerical Analysis, Wiley, 1982.
%  (Note that if some fk > 1/eps, with eps the machine epsilon,
%  the value of eps in the code may have to be reduced.)



  fk = fk(:); x = x(:);                    % Make sure data are column vectors.

   N = length(fk); 
   M = length(x);
     
  xk = sin(pi*[N-1:-2:1-N]'/(2*(N-1)));    % Compute Chebyshev points.

   w = ones(N,1).*(-1).^[0:N-1]';          % w = weights for Chebyshev formula
   w(1) = w(1)/2; w(N) = w(N)/2;
 
   D = x(:,ones(1,N)) - xk(:,ones(1,M))';  % Compute quantities x-x(k)
   D = 1./(D+eps*(D==0));                  % and their reciprocals.
  
   p = D*(w.*fk)./(D*w);                   % Evaluate interpolant as
                                           % matrix-vector products.
