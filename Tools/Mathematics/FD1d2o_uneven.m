function [ D ] = FD1d2o_uneven(x,y)
% First derivative of numerical data "y" at "x", with respect to "x".
% Suited for uneven spacing of x, not fundamentally second order.
% http://websrv.cs.umt.edu/isis/index.php/Finite_differencing:_Introduction
% Data has to be provided in columns. If a matrix is inserted, derivatives are
% taken along each columns for all columns separately.

dx = diff(x);
D = y;

% Central differences
D(2:end-1,:) = (   dx(1:end-1,:).^2                  .* y(3:end  ,:) ...
                + (dx(2:end  ,:).^2-dx(1:end-1,:).^2).* y(2:end-1,:) ...
                -  dx(2:end  ,:).^2                  .* y(1:end-2,:)  ) ...
               ./ (dx(1:end-1,:).*dx(2:end,:).*(dx(1:end-1,:) + dx(2:end,:)));
          
% Forward differences
D(1,:)       = (-  dx(1,:).^2                        .* y(3,:) ...
                + (dx(1,:)    +    dx(2,:)).^2       .* y(2,:) ...
                - (dx(2,:).^2 + 2.*dx(1,:).*dx(2,:)) .* y(1,:)  ) ...
               ./ (dx(1,:).*dx(2,:).*(dx(1,:) + dx(2,:)));
    
% Backward differences
D(end,:)     = (  (dx(end-1,:).^2 + 2.*dx(end-1,:).*dx(end,:)) .* y(end,:) ...
                - (dx(end-1,:)    +    dx(end  ,:)).^2         .* y(end-1,:) ...
                +  dx(end  ,:).^2                              .* y(end-2,:)   ) ...
               ./ (dx(end-1,:).*dx(end,:).*(dx(end-1,:)+dx(end,:)));

end

