function [ f ] = NLT_HNS(StabGrid,RunJ,StabRes,D1,HB)
%% License 

%% Description

% Calculates the nonlinear forcing planes for all active modes presented in RunJ
% on the grid defined in StabGrid for perturbation velocities presented in
% u, v, and w. 

% Inputs
% StabGrid  Structure containing grid information
% RunJ      (1, at most (2N+1)x(2M+1))  Modes to run for
% StabRes   Stability results
% D1        (ny,ny)  First-order wall-normal derivative operator
% HB        ((2N+1)x(2M+1),(2N+1)x(2M+1),(2N+1)x(2M+1))  Harmonic Balancing matrices

% Outputs
% f ((2N+1)x(2M+1), 4 ny nx) Forcing term

%% Find further necessary values from inputs

% Find nf, ny and nx
[nf,ny,nx] = size(StabRes.u);

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
dxi = StabGrid.xiun(2)-StabGrid.xiun(1);

%% Calculate derivative fields

% Calculate derivatives wrt xi
for j = RunJ
dudxi(j,:,:) = -FD1d4o(squeeze(StabRes.u(j,:,:)),dxi);
dvdxi(j,:,:) = -FD1d4o(squeeze(StabRes.v(j,:,:)),dxi);
dwdxi(j,:,:) = -FD1d4o(squeeze(StabRes.w(j,:,:)),dxi);
end

% Calculate derivatives wrt eta
for j = RunJ
for i = 1:nx
   dudeta(j,:,i) =  (D1*StabRes.u(j,:,i)')';
   dvdeta(j,:,i) =  (D1*StabRes.v(j,:,i)')';
   dwdeta(j,:,i) =  (D1*StabRes.w(j,:,i)')';
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
        ududxi    = StabRes.u(gh,:,:) .* dudxi(jk,:,:);     
        ududeta   = StabRes.u(gh,:,:) .* dudeta(jk,:,:);       
        wu        = StabRes.w(gh,:,:) .* StabRes.u(jk,:,:);      
        vdudeta   = StabRes.v(gh,:,:) .* dudeta(jk,:,:);
        vdudxi    = StabRes.v(gh,:,:) .* dudxi(jk,:,:);   
        
        % Set up y-momentum source terms
        udvdxi    = StabRes.u(gh,:,:) .* dvdxi(jk,:,:);   
        vdvdeta   = StabRes.v(gh,:,:) .* dvdeta(jk,:,:);   
        wv        = StabRes.w(gh,:,:) .* StabRes.v(jk,:,:);       
        udvdeta   = StabRes.u(gh,:,:) .* dvdeta(jk,:,:);
        vdvdxi    = StabRes.v(gh,:,:) .* dvdxi(jk,:,:);  
        
        % Set up z-momentum source terms
        udwdxi    = StabRes.u(gh,:,:) .* dwdxi(jk,:,:);     
        vdwdeta   = StabRes.v(gh,:,:) .* dwdeta(jk,:,:);   
        ww        = StabRes.w(gh,:,:) .* StabRes.w(jk,:,:);       
        udwdeta   = StabRes.u(gh,:,:) .* dwdeta(jk,:,:);
        vdwdxi    = StabRes.v(gh,:,:) .* dwdxi(jk,:,:);     
        
        % Define beta of the mode jk
        beta  = StabRes.betavec(jk);
    
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



