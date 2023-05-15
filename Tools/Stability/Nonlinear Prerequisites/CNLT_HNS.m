function [ f ] = CNLT_HNS(StabGrid,RunJ,u,v,w,D1,HB,xi,betavec)
%% License 

%% Description

% Calculates the nonlinear forcing planes for all active modes presented in RunJ
% on the grid defined in StabGrid for perturbation velocities presented in
% u, v, and w. 

% Inputs
% StabGrid 

% Outputs

%% Find further necessary values from inputs

% Find nf, ny and nx
[nf,ny,nx] = size(u);

% Define imaginary unit
iu = sqrt(-1);

%reset f
f = zeros(nf,4*ny*nx);

% Rename transformation coefficients to improve readability
etay  = StabGrid.etay;
etax  = StabGrid.etax;
xiy   = StabGrid.xiy;
xix   = StabGrid.xix;

% Find step size in numerical domain 
dxi = xi(2)-xi(1);

%% Calculate derivative fields

% Calculate derivatives wrt xi
for j = RunJ
dudxi(j,:,:) = -FD1d4o(squeeze(u(j,:,:)),dxi);
dvdxi(j,:,:) = -FD1d4o(squeeze(v(j,:,:)),dxi);
dwdxi(j,:,:) = -FD1d4o(squeeze(w(j,:,:)),dxi);
end

% Calculate derivatives wrt eta
for j = RunJ
for i = 1:nx
   dudeta(j,:,i) =  (D1*u(j,:,i)')';
   dvdeta(j,:,i) =  (D1*v(j,:,i)')';
   dwdeta(j,:,i) =  (D1*w(j,:,i)')';
end
end

%% Calculate forcing

% Reset/initialize momentum sums
xmom = zeros(nf,ny,nx);
ymom = zeros(nf,ny,nx);
zmom = zeros(nf,ny,nx);

% Loop only over active modes, skip symmetries
for j = RunJ(RunJ>=round(nf/2))

    % Find modes that interact to force mode j from Harmonic Balancing (HB)
    [gharray, jkarray] = find(HB(:,:,j)); 
    
    % Delete modes gh and jk that are inactive
    gharray = ismember(gharray,RunJ).*gharray; 
    jkarray = ismember(jkarray,RunJ).*jkarray;
    [INT,~] = size(nonzeros(gharray.*jkarray)); %If both modes are nonzero loop
    
    % Delete interactions with 0-amplitude modes
    R = find(gharray.*jkarray ~=0); 
    gharray = gharray(R);
    jkarray = jkarray(R);

    for j2 = 1:INT 
        % Determine interacting modes
        gh = gharray(j2); 
        jk = jkarray(j2);
                
        % Set up x-momentum source terms
        ududxi    = u(gh,:,:) .* dudxi(jk,:,:);     
        ududeta   = u(gh,:,:) .* dudeta(jk,:,:);       
        wu        = w(gh,:,:) .* u(jk,:,:);      
        vdudeta   = v(gh,:,:) .* dudeta(jk,:,:);
        vdudxi    = v(gh,:,:) .* dudxi(jk,:,:);   
        
        % Set up y-momentum source terms
        udvdxi    = u(gh,:,:) .* dvdxi(jk,:,:);   
        vdvdeta   = v(gh,:,:) .* dvdeta(jk,:,:);   
        wv        = w(gh,:,:) .* v(jk,:,:);       
        udvdeta   = u(gh,:,:) .* dvdeta(jk,:,:);
        vdvdxi    = v(gh,:,:) .* dvdxi(jk,:,:);  
        
        % Set up z-momentum source terms
        udwdxi    = u(gh,:,:) .* dwdxi(jk,:,:);     
        vdwdeta   = v(gh,:,:) .* dwdeta(jk,:,:);   
        ww        = w(gh,:,:) .* w(jk,:,:);       
        udwdeta   = u(gh,:,:) .* dwdeta(jk,:,:);
        vdwdxi    = v(gh,:,:) .* dwdxi(jk,:,:);     
        
        % Define beta of the mode jk
        beta  = betavec(jk);
    
        % Calculate x, y, and z-momentum forcing
        xmom(j,:,:) = squeeze(xmom(j,:,:)) -squeeze(ududxi).*xix-squeeze(ududeta).*etax  -squeeze(vdudxi).*xiy-squeeze(vdudeta).*etay  -iu*beta.*squeeze(wu);%
        ymom(j,:,:) = squeeze(ymom(j,:,:)) -squeeze(udvdxi).*xix-squeeze(udvdeta).*etax  -squeeze(vdvdxi).*xiy-squeeze(vdvdeta).*etay  -iu*beta.*squeeze(wv);%
        zmom(j,:,:) = squeeze(zmom(j,:,:)) -squeeze(udwdxi).*xix-squeeze(udwdeta).*etax  -squeeze(vdwdxi).*xiy-squeeze(vdwdeta).*etay  -iu*beta.*squeeze(ww);%
    
     end % end of mode interaction loop
    
        % Reshape into f
        for i = 2:nx %no forcing at i = 1!
            f(j,(i-1)*4*ny+(1:4*ny)) = [squeeze(xmom(j,:,i)) squeeze(ymom(j,:,i)) squeeze(zmom(j,:,i)) zeros(1,ny)];
        end % end of i-loop

end % end of recipient mode loop

end % end of function



