function [eigenvalue,index] = EVfilter(EV,V,y,omega,beta,delta99,Mach)
% Function which filters out incorrect eigenvalues from an eigenvalue
% spectrum by checking for an exponential fit
%
% The physical disturbance will decay exponentially along "y". Therefore, 
% each eigenfunction is compared to:
%       phi = exp( - i sqrt(alpha_r^2 + beta^2) y )
% where alpha is the eigenvalue corresponding to each eigenfunction 
% (real component only) and beta the spanwise wavenumber used as input 
% for the OS solver.
%
% It is required that the OS solver used a mapping that includes the full 
% exponential decay of the eigenfunction! Therefore, make sure that y_max 
% of the mapping is sufficiently high.
%
% M. Kotsonis, K. Groot & J.Y. Boersma - 2017

if nargin<7
    Mach = 0;
end

power1 = 6;     % Tweaking factor for the slope error
power2 = 2;     % Tweaking factor for the switch error
treshold = 14;  % logarithm of accuracy of the numerical eigenfunction

use = ~isinf(EV);
EV = EV(use);
V = V(:,use);

NC = length(y);     % length of Chebychev matrix
NE = length(EV);    % length of eigenspectrum

% Parse eigenfunctions V
[Vmax,~] = max(abs(V)); % Store maxima for later use
V = V./repmat(Vmax,NC,1); % Normalize
V(isnan(V))=0+0*sqrt(-1);% replace NaN with zeros

% Combined wavenumbers of alpha and beta
% wavenumber = sqrt(real(EV).^2 + beta^2 - (real(EV)-omega).^2 * Mach)';
wavenumber = real(sqrt(EV.^2 + beta^2 - (EV-omega).^2 * Mach^2))';
%wavenumber = real(sqrt(EV.^2 + beta^2))';
% the former result is derived by Balakumar & Malik (1992) - Discrete modes
% and continuous spectra in supersonic boundary layers, eqs. (3.7) & (2.19)
% note that it is assumed that Pr = 3/4 and Re >> 1, in other cases,
% implement equation (2.22)

% Loop over eigenvalues (and functions)
Error = zeros(1,NE);
SlopeError = Error;
Switches = Error;
y_top = Error;
y_bulge = Error;
for n = 1:NE
    %% Limits in y

    % Determine at which y the numerical accuracy of the eigenfunction will no
    % longer produce an exponential trendline
    y_top(n) = treshold/wavenumber(n);

    % Determine the point where the exponential decay starts and the bulge ends
    % (this point is not strictly defined, just an arbitrary height)
    y_bulge(n) = delta99 + 0.15*delta99/wavenumber(n);
    
    % Check the values
    if y_bulge(n) > y(1)/2
        y_bulge(n) = y(1)/2;
    end
    if y_top(n) < y_bulge(n)*2
        y_top(n) = inf;
    end


    %% Fit of the exponential trend
    % Exponential decay should be equal to the wavenumber

    % Selection "s1" selects the freestream
    s1 = max([find(y<=y_top(n),  1),        find(y<=max(y)*0.85,1)]) ...    % Always skip upper 15%
       : min([find(y>=y_bulge(n),1,'last'), NC-1]);                         % Always skip last node
    EF = log(abs(V(s1,n)));
   
    % Per-segment comparison
    SlopeDifference = wavenumber(n) + diff(EF)./diff(y(s1));
    SegmentLength = -diff(y(s1)); % Weighting factor
    SlopeError_Segment = mean(abs(SegmentLength.*SlopeDifference.^2));

    % Linear least-square fit to determine the slope
    EF_fit = polyfit(y(s1),EF,1);
    SlopeError_LS = abs(wavenumber(n)+EF_fit(1));
    
    % Combine segment comparison and least-square fit
    SlopeError(n) = sqrt(SlopeError_Segment.*SlopeError_LS).^(power1);

    %% Number of switches in gradient

    % Selection "s2" selects the bulge
    s2 = find(y>=y_bulge(n),1,'last'):NC;

    % Count switches, reduced by expected values for each component
    S_imag = countSwitches(imag(V(s2,n)),2);
    S_real = countSwitches(real(V(s2,n)),1);
    S_abs  = countSwitches( abs(V(s2,n)),1);

    % The sum of the number of switches should be zero
    Switches(n) = S_imag+S_real+S_abs;

    % Scale the slope error with the number of switches, using an exponent 
    % because of the logarithmic scaling.
    % Times "power" amplifies the error of multiple switches.
    Error(n) = SlopeError(n).*10.^(Switches(n)*power2);
    
    
    
    
%     %% Per EF plots
%     figure(3)
%     clf
%     
%     h1=subplot(1,2,1);
%     plot(wavenumber,'-x')
%     h1.YLim(1)=0;
%     
%     h2=subplot(1,2,2);
%     hold on
%     plot(log(abs(V(:,n))),y,'-x')
%     plot(-wavenumber(n)*y(s1),y(s1),'-o')
%     if y_top(n)<y(1)
%         plot(h2.XLim,y_top(n)*[1 1],'k--')
%     end
%     if y_bulge(n)<y(1)
%         plot(h2.XLim,y_bulge(n)*[1 1],'b--','linewidth',1.5)
%     end
%     hold off
%     grid on
%     title(sprintf('Top=%.2f, Bulge=%.2f',y_top(n),y_bulge(n)));
%     keyboard
    

end

%% Select eigenvalue

% Minimum error is the eigenvalue to keep
[~,keep] = min(Error);
eigenvalue = EV(keep);

% Select correct index of original EV list
remain = find(use);
index = remain(keep);


%% Plots


% figure(2)
% clf
% hold on
% plot(1:NE,log10(SlopeError))
% plot(1:NE,Switches)
% plot(1:NE,log10(Error))
% hold off
% grid on
% h=gca;
% h.YLim = [min(log10(SlopeError))-2, 10];
% xlabel('Index of eigenfunction')
% legend('Slopes of exponential (log)','Number of switches','Scaled with switches (log)')
% 
% figure(10)
% EF_wavenumber = -repmat(wavenumber,NC,1).*repmat(y,1,NE);
% [~,ordered] = sort(Error);
% % [~,ordered] = sort(SlopeError);
% % [~,ordered] = sort(switches);
% 
% for i = ordered(1)
%     
%     clf
%     subplot(1,2,1)
%     hold on
%     plot(abs(V(:,i)),y,real(V(:,i)),y,imag(V(:,i)),y)
% %     plot(GR(:,i)./max(abs(GR(:,i))),y(s3),'x-')
%     hold off
%     grid on
%     title(sprintf('Index: %.0f',remain(i)))
%     xlabel('$\phi$','interpreter','latex')
%     ylabel('$y$','interpreter','latex','rot',0,...
%         'horizontalalignment','right','verticalalignment','middle')
% 
%     h2=subplot(1,2,2);
%     hold on
%     plot(log(abs(V(:,i))),y,'-x')
%     plot(EF_wavenumber(:,i),y,'-o')
%     plot(h2.XLim,y_top(i)*[1 1],'k--')
%     plot(h2.XLim,y_bulge(i)*[1 1],'b--','linewidth',1.5)
%     hold off
%     title('Logarithmic fit')
%     h=legend('abs( $\phi$ )',...
%         'abs( $e^{-i \sqrt{\alpha_r^2+\beta^2} y}$ )',...
%         'Top limit',...
%         'Bulge limit',...
%         'location','best');
%     h.Interpreter = 'latex';
%     h.FontSize = 11;
%     grid on
%     xlabel('$ln(\phi)$','interpreter','latex')
%     ylabel('$y$','interpreter','latex','rot',0,...
%         'horizontalalignment','right','verticalalignment','middle')
%     pause(0.01)
%     keyboard
% end


end

function switches = countSwitches(data,expected)
% Count the number of 'switches' in the gradient of 'data', reduced by the
% number of 'expected' switches, up to a minimum of zero.

% Sign of difference of the data segments (hence, removes one element)
sn = sign(diff(data));
% Product of succeeding data segments (also removes one element)
count = sn(1:end-1,:).*sn(2:end,:); % -1 if switched, 1 if not
% Count number of switches
switches = sum(count+1==0);
% Reduce number of switches
switches = abs(switches-expected);
% % Apply minimum of zero
% switches(switches<0) = 0;
end
