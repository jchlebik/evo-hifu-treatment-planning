function sensor_data = kspaceElastic2D(kgrid, medium, source, sensor, varargin)
%KSPACEELASTIC2D 2D time-domain simulation of elastic wave propagation.
%
% DESCRIPTION:
%     kspaceElastic2D implements the elastic k-space model described in B.
%     E. Treeby, B. T. Cox (2014): Modeling power law absorption and
%     dispersion in viscoelastic solids using a split-field and the
%     fractional Laplacian. In: J. Acoust. Soc. Am., 136 (4), pp.
%     1499-1510, 2014.
%
% INPUTS:
%     see pstdElastic2D for general inputs
% 
%     medium.sound_speed_compression
%     medium.sound_speed_shear
%     medium.sound_speed_ref_compression
%     medium.sound_speed_ref_shear

%     medium.density
%
%     medium.alpha_coeff_compression
%     medium.alpha_coeff_shear
%     medium.alpha_power_compression
%     medium.alpha_power_shear
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 16th January 2014
%     last update - 19th February 2017

% suppress mlint warnings that arise from using subscripts
%#ok<*NASGU>
%#ok<*COLND>
%#ok<*NODEF>
%#ok<*INUSL>

% check whether code should extract absorption, and replace any output data
% with this data
if any(strcmp('PlotAbsorption', varargin))
   
    % find the location of the optional input in the list
    index = find(strcmp('PlotAbsorption', varargin));
    
    % read the value of the parameter
    manually_check_absorption = varargin{index + 1};
        
    % delete the optional input from the input arguments list
    varargin(index:index+1) = [];
    
    % clean-up
    clear index;
    
else
    
    % set parameter to false
    manually_check_absorption = false;
    
end

% =========================================================================
% CHECK INPUT STRUCTURES AND OPTIONAL INPUTS
% =========================================================================

% start the timer and store the start time
start_time = clock;
tic;

% set the name of the simulation code
MFILE = mfilename;

% get the number of inputs and outputs (nargin and nargout can't be used in
% subscripts in MATLAB 2016b or later)
num_inputs  = nargin;
num_outputs = nargout;

% define literals
getkWaveDefaults;

% run subscript to check inputs
kspaceFirstOrder_inputChecking;

% assign short names for the sound speed
cp = medium.sound_speed_compression;
cs = medium.sound_speed_shear;

% assign the lame parameters, pre-multiplied by the time step
mu_dt     = ( cs.^2 .* rho0 ) .* kgrid.dt;
lambda_dt = ( cp.^2 .* rho0 ) .* kgrid.dt - 2*mu_dt;

% assign the viscosity coefficients
if kelvin_voigt_model

    % extract and convert the absorption parameters to the correct units
    alpha0_p_db = medium.alpha_coeff_compression;
    alpha0_s_db = medium.alpha_coeff_shear;
    y_p         = medium.alpha_power_compression;
    y_s         = medium.alpha_power_shear;

    alpha0_p = db2neper(alpha0_p_db, y_p);
    alpha0_s = db2neper(alpha0_s_db, y_s);
    
    % assign the viscosity coefficients, pre-multiplied by the time step
    eta_dt = ( - 2 .* rho0 .* cs .^3 .* alpha0_s ./ cos(pi*y_s/2) ) .* kgrid.dt;
    chi_dt = ( - 2 .* rho0 .* cp .^3 .* alpha0_p ./ cos(pi*y_p/2) ) .* kgrid.dt - 2 * eta_dt;
    
    % check value for chi isn't negative
    if min(chi_dt(:)) < 0
        disp('  WARNING: CHI is negative!');
    end
    
    if min(chi_dt(:) + 2*eta_dt(:)) < 0
        disp('  WARNING: CHI + 2*ETA is negative!');
    end
     
end
    
% =========================================================================
% CALCULATE MEDIUM PROPERTIES ON STAGGERED GRID
% =========================================================================

% rho0_sgx = rho0;
% rho0_sgy = rho0;
% mu_dt_sgxy = mu;
% cp_sgx = cp;
% cp_sgy = cp;
% cs_sgx = cs;
% cs_sgy = cs;   
% eta_sgxy = eta;

% calculate the values of the density at the staggered grid points
% using linear interpolation, where sgx  = (x + dx/2, y) and 
% sgy  = (x, y + dy/2)
if numDim(rho0) == 2 && use_sg

    % rho0 is heterogeneous and staggered grids are used
    rho0_sgx = interpn(kgrid.x, kgrid.y, rho0, kgrid.x + kgrid.dx/2, kgrid.y,              '*linear');
    rho0_sgy = interpn(kgrid.x, kgrid.y, rho0, kgrid.x,              kgrid.y + kgrid.dy/2, '*linear');
    
    % set values outside of the interpolation range to original values 
    rho0_sgx(isnan(rho0_sgx)) = rho0(isnan(rho0_sgx));
    rho0_sgy(isnan(rho0_sgy)) = rho0(isnan(rho0_sgy));
    
else
    % rho0 is homogeneous or staggered grids are not used
    rho0_sgx = rho0;
    rho0_sgy = rho0;
end

% calculate the values of mu at the staggered grid points using linear
% interpolation where sgxy = (x + dx/2, y + dy/2)
if numDim(mu_dt) == 2 && use_sg

    % mu is heterogeneous and staggered grids are used
    mu_dt_sgxy  = interpn(kgrid.x, kgrid.y, mu_dt, kgrid.x + kgrid.dx/2, kgrid.y + kgrid.dy/2, '*linear');
    
    % set values outside of the interpolation range to original values 
    mu_dt_sgxy(isnan(mu_dt_sgxy))  = mu_dt(isnan(mu_dt_sgxy));
    
else
    % mu is homogeneous or staggered grids are not used
    mu_dt_sgxy = mu_dt;
end

% calculate the values of cp at the staggered grid points using linear
% interpolation 
if numDim(cp) == 2
    
    % cp is heterogeneous
    cp_sgx = interpn(kgrid.x, kgrid.y, cp, kgrid.x + kgrid.dx/2, kgrid.y, '*linear');
    cp_sgy = interpn(kgrid.x, kgrid.y, cp, kgrid.x, kgrid.y + kgrid.dy/2, '*linear');
    
    % set values outside of the interpolation range to original values 
    cp_sgx(isnan(cp_sgx)) = cp(isnan(cp_sgx));
    cp_sgy(isnan(cp_sgy)) = cp(isnan(cp_sgy));
    
else
    
    % cp is homogeneous or staggered grids are not used
    cp_sgx = cp;
    cp_sgy = cp;
    
end

% calculate the values of cs at the staggered grid points using linear
% interpolation 
if numDim(cs) == 2
    
    % cs is heterogeneous
    cs_sgx = interpn(kgrid.x, kgrid.y, cs, kgrid.x + kgrid.dx/2, kgrid.y, '*linear');
    cs_sgy = interpn(kgrid.x, kgrid.y, cs, kgrid.x, kgrid.y + kgrid.dy/2, '*linear');
    
    % set values outside of the interpolation range to original values 
    cs_sgx(isnan(cs_sgx)) = cs(isnan(cs_sgx));
    cs_sgy(isnan(cs_sgy)) = cs(isnan(cs_sgy));
    
else
    
    % cs is homogeneous
    cs_sgx = cs;
    cs_sgy = cs;
    
end

% calculate the values of eta at the staggered grid points using linear
% interpolation 
if kelvin_voigt_model 
    if numDim(eta_dt) == 2

        % eta is heterogeneous and staggered grids are used
        eta_sgxy  = interpn(kgrid.x, kgrid.y, eta_dt, kgrid.x + kgrid.dx/2, kgrid.y + kgrid.dy/2, '*linear');

        % set values outside of the interpolation range to original values 
        eta_sgxy(isnan(eta_sgxy)) = eta_dt(isnan(eta_sgxy));

    else
        % eta is homogeneous or staggered grids are not used
        eta_sgxy = eta_dt;
    end
end

% =========================================================================
% PREPARE DERIVATIVE AND PML OPERATORS
% =========================================================================

% define the k-space derivative operators, multiply by the staggered
% grid shift operators, and then re-order using ifftshift (the option
% use_sg exists for debugging) 
if use_sg
    
    % define the shift derivative operators
    ddx_k_shift_pos = ifftshift( 1i*kgrid.kx_vec .* exp( 1i*kgrid.kx_vec*kgrid.dx/2) );
    ddx_k_shift_neg = ifftshift( 1i*kgrid.kx_vec .* exp(-1i*kgrid.kx_vec*kgrid.dx/2) );
    ddy_k_shift_pos = ifftshift( 1i*kgrid.ky_vec .* exp( 1i*kgrid.ky_vec*kgrid.dy/2) );
    ddy_k_shift_neg = ifftshift( 1i*kgrid.ky_vec .* exp(-1i*kgrid.ky_vec*kgrid.dy/2) );
    
    % define the shift operators neeed for the k-space dyadic (see eqn 39b
    % in Kamyar's paper)
    shift_neg_y_pos_x = ifftshift( exp(1i*kgrid.kx*kgrid.dx/2 - 1i*kgrid.ky*kgrid.dy/2) );
    shift_neg_x_pos_y = ifftshift( exp(1i*kgrid.ky*kgrid.dy/2 - 1i*kgrid.kx*kgrid.dx/2) );
    
else
    
    % define the derivative operators _without_ shifting
    ddx_k_shift_pos = ifftshift( 1i*kgrid.kx_vec );
    ddx_k_shift_neg = ifftshift( 1i*kgrid.kx_vec );
    ddy_k_shift_pos = ifftshift( 1i*kgrid.ky_vec );
    ddy_k_shift_neg = ifftshift( 1i*kgrid.ky_vec );
    
    % define the operators neeed for the k-space dyadic _without_ shifting
    shift_neg_y_pos_x = 1;
    shift_neg_x_pos_y = 1;

end

% force to be in the correct direction for use with BSXFUN
ddy_k_shift_pos = ddy_k_shift_pos.'; 
ddy_k_shift_neg = ddy_k_shift_neg.';

% define normalised wavenumber vectors
kx_norm               = kgrid.kx ./ kgrid.k;
kx_norm(kgrid.k == 0) = 0;
kx_norm               = ifftshift(kx_norm);

ky_norm               = kgrid.ky ./ kgrid.k;
ky_norm(kgrid.k == 0) = 0;
ky_norm               = ifftshift(ky_norm);

if kelvin_voigt_model

    % define the absorption derivative operators
    Lp1_k     =   ifftshift(kgrid.k.^(y_p-1));
    Lp1_x_sgx =   cp_sgx.^(y_p-1) .* sin(pi*y_p/2);
    Lp1_x_sgy =   cp_sgy.^(y_p-1) .* sin(pi*y_p/2);

    Lp2_k     =   ifftshift(kgrid.k.^(y_p-2));
    Lp2_x_sgx = - cp_sgx.^(y_p-2) .* cos(pi*y_p/2);
    Lp2_x_sgy = - cp_sgy.^(y_p-2) .* cos(pi*y_p/2);

    Ls1_k     =   ifftshift(kgrid.k.^(y_s-1));
    Ls1_x_sgx =   cs_sgx.^(y_s-1) .* sin(pi*y_s/2);
    Ls1_x_sgy =   cs_sgy.^(y_s-1) .* sin(pi*y_s/2);

    Ls2_k     =   ifftshift(kgrid.k.^(y_s-2));
    Ls2_x_sgx = - cs_sgx.^(y_s-2) .* cos(pi*y_s/2);
    Ls2_x_sgy = - cs_sgy.^(y_s-2) .* cos(pi*y_s/2);

    % remove singularity
    Lp1_k(isinf(Lp1_k)) = 0;
    Ls1_k(isinf(Ls1_k)) = 0;
    Lp2_k(isinf(Lp2_k)) = 0;
    Ls2_k(isinf(Ls2_k)) = 0;

end
    
% create k-space operator (the option use_kspace exists for debugging)
if use_kspace
    kappa_p = ifftshift(sinc(c_ref_compression*kgrid.dt*kgrid.k/2));
    kappa_s = ifftshift(sinc(c_ref_shear      *kgrid.dt*kgrid.k/2));
else
    kappa_p = 1;
    kappa_s = 1;
end

% clear the staggered grid sound speed variables
clear cp_sgx cp_sgy cs_sgx cs_sgy;

% =========================================================================
% PLOT THINGS
% =========================================================================

% figure
% hz_plots = 6;
% vt_plots = 4;
% plot_num = 0;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(lambda_dt);
% title('lambda_dt');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(mu_dt);
% title('mu_dt');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(mu_dt_sgxy);
% title('mu_dt_sgxy');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(rho0);
% title('rho0');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(rho0_sgx);
% title('rho0_sgx');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(rho0_sgy);
% title('rho0_sgy');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(eta_dt);
% title('eta_dt');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(eta_sgxy);
% title('eta_sgxy');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(chi_dt);
% title('chi_dt');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(Lp1_k);
% title('Lp1_k');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(Lp1_x_sgx);
% title('Lp1_x_sgx');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(Lp1_x_sgy);
% title('Lp1_x_sgy');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(Lp2_k);
% title('Lp2_k');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(Lp2_x_sgx );
% title('Lp2_x_sgx ');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(Lp2_x_sgy);
% title('Lp2_x_sgy');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(Ls1_k);
% title('Ls1_k');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(Ls1_x_sgx );
% title('Ls1_x_sgx ');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(Ls1_x_sgy);
% title('Ls1_x_sgy');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(Ls2_k );
% title('Ls2_k ');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(Ls2_x_sgx);
% title('Ls2_x_sgx');
% colorbar;
% 
% plot_num = plot_num + 1;
% subplot(vt_plots, hz_plots, plot_num);
% imagesc(Ls2_x_sgy);
% title('Ls2_x_sgy');
% colorbar;

% =========================================================================
% DATA CASTING
% =========================================================================

% preallocate the loop variables using the castZeros anonymous function
% (this creates a matrix of zeros in the data type specified by data_cast)
ux_sgx          = castZeros([kgrid.Nx, kgrid.Ny]);  % could be removed, see notes below 
uy_sgy          = castZeros([kgrid.Nx, kgrid.Ny]);  % could be removed, see notes below

ux_p_k_sgx      = castZeros([kgrid.Nx, kgrid.Ny]);  % could be removed, see notes below 
uy_p_k_sgy      = castZeros([kgrid.Nx, kgrid.Ny]);  % could be removed, see notes below

ux_s_k_sgx     	= castZeros([kgrid.Nx, kgrid.Ny]);  % could be removed, see notes below 
uy_s_k_sgy     	= castZeros([kgrid.Nx, kgrid.Ny]);  % could be removed, see notes below

sxx_p           = castZeros([kgrid.Nx, kgrid.Ny]);
syy_p           = castZeros([kgrid.Nx, kgrid.Ny]);
sxy_p_sgxy    	= castZeros([kgrid.Nx, kgrid.Ny]);

sxx_s           = castZeros([kgrid.Nx, kgrid.Ny]);
syy_s           = castZeros([kgrid.Nx, kgrid.Ny]);
sxy_s_sgxy   	= castZeros([kgrid.Nx, kgrid.Ny]);

duxdx_p_k       = castZeros([kgrid.Nx, kgrid.Ny]);  % could be removed, see notes below
duxdy_p_k_sgxy	= castZeros([kgrid.Nx, kgrid.Ny]);  % could be removed, see notes below
duydy_p_k       = castZeros([kgrid.Nx, kgrid.Ny]);  % could be removed, see notes below
duydx_p_k_sgxy  = castZeros([kgrid.Nx, kgrid.Ny]);  % could be removed, see notes below

duxdx_s_k       = castZeros([kgrid.Nx, kgrid.Ny]);  % could be removed, see notes below
duxdy_s_k_sgxy  = castZeros([kgrid.Nx, kgrid.Ny]);  % could be removed, see notes below
duydy_s_k       = castZeros([kgrid.Nx, kgrid.Ny]);  % could be removed, see notes below
duydx_s_k_sgxy  = castZeros([kgrid.Nx, kgrid.Ny]);  % could be removed, see notes below

% run subscript to cast the remaining loop variables to the data type
% specified by data_cast 
if ~strcmp(data_cast, 'off')
    
    % throw error for unsupported feature
    error('''DataCast'' input is not currently supported');
    
end

% =========================================================================
% CREATE INDEX VARIABLES
% =========================================================================

% setup the time index variable
if ~record.time_rev
    index_start = 1;
    index_step = 1;
    index_end = kgrid.Nt; 
else
    % throw error for unsupported feature
    error('Time reversal is not currently supported');
    
    % reverse the order of the input data
    sensor.time_reversal_boundary_data = fliplr(sensor.time_reversal_boundary_data); %#ok<UNRCH>
    index_start = 1;
    index_step = 1;
    
    % stop one time point before the end so the last points are not
    % propagated
    index_end = kgrid.Nt - 1;    
end

% =========================================================================
% PREPARE VISUALISATIONS
% =========================================================================

% pre-compute suitable axes scaling factor
if plot_layout || plot_sim
    [~, scale, prefix] = scaleSI(max([kgrid.x_vec; kgrid.y_vec])); 
end

% throw error for currently unsupported plot layout feature
if plot_layout
    error('''PlotLayout'' input is not currently supported');   
end

% initialise the figures used for animation if 'PlotSim' is set to 'true'
if plot_sim
    kspaceFirstOrder_initialiseFigureWindow;
end 

% initialise movie parameters if 'RecordMovie' is set to 'true'
if record_movie
    kspaceFirstOrder_initialiseMovieParameters;
end

% =========================================================================
% LOOP THROUGH TIME STEPS
% =========================================================================

% update command line status
disp(['  precomputation completed in ' scaleTime(toc)]);
disp('  starting time loop...');

% restart timing variables
loop_start_time = clock;
tic;

% start time loop
for t_index = index_start:index_step:index_end
        
    % calculate the velocity at the next time step using the components of
    % the split stress tensor at the current time step
    ux_sgx = ux_sgx + dt./rho0_sgx .* real(ifftn(...
        bsxfun(@times, ddx_k_shift_pos, kappa_p .* fftn(sxx_p)) + ...       % dsxxdx_p_sgx
        bsxfun(@times, ddy_k_shift_neg, kappa_p .* fftn(sxy_p_sgxy)) + ...  % dsxydy_p_sgx
        bsxfun(@times, ddx_k_shift_pos, kappa_s .* fftn(sxx_s)) + ...       % dsxxdx_s_sgx
        bsxfun(@times, ddy_k_shift_neg, kappa_s .* fftn(sxy_s_sgxy))...     % dsxydy_s_sgx
        ));
    
    uy_sgy = uy_sgy + dt./rho0_sgy .* real(ifftn(...
        bsxfun(@times, ddy_k_shift_pos, kappa_p .* fftn(syy_p)) + ...       % dsyydy_p_sgy
        bsxfun(@times, ddx_k_shift_neg, kappa_p .* fftn(sxy_p_sgxy)) + ...  % dsxydx_p_sgy        
        bsxfun(@times, ddy_k_shift_pos, kappa_s .* fftn(syy_s)) + ...       % dsyydy_s_sgy
        bsxfun(@times, ddx_k_shift_neg, kappa_s .* fftn(sxy_s_sgxy)) ...    % dsxydx_s_sgy
        ));
             
    % add in the velocity source terms
    if ux_source >= t_index
        if strcmp(source.u_mode, 'dirichlet')
            % enforce the source values as a dirichlet boundary condition
            ux_sgx(u_source_pos_index) = source.ux(:, t_index);
        else
            % add the source values to the existing field values 
            ux_sgx(u_source_pos_index) = ux_sgx(u_source_pos_index) + source.ux(:, t_index);
        end
    end
    if uy_source >= t_index
        if strcmp(source.u_mode, 'dirichlet')
            % enforce the source values as a dirichlet boundary condition        
            uy_sgy(u_source_pos_index) = source.uy(:, t_index);
        else
            % add the source values to the existing field values 
            uy_sgy(u_source_pos_index) = uy_sgy(u_source_pos_index) + source.uy(:, t_index);
        end
    end    
    
    % split the velocity into compressional (p) and shear (s) parts
    ux_p_k_sgx = kx_norm.^2       .* fftn(ux_sgx)   + kx_norm .* ky_norm .* shift_neg_y_pos_x .* fftn(uy_sgy);
    ux_s_k_sgx = (1 - kx_norm.^2) .* fftn(ux_sgx)   - kx_norm .* ky_norm .* shift_neg_y_pos_x .* fftn(uy_sgy);
    uy_p_k_sgy = ky_norm.^2       .* fftn(uy_sgy)   + ky_norm .* kx_norm .* shift_neg_x_pos_y .* fftn(ux_sgx);
    uy_s_k_sgy = (1 - ky_norm.^2) .* fftn(uy_sgy)   - ky_norm .* kx_norm .* shift_neg_x_pos_y .* fftn(ux_sgx);

    % manually store the split velocity at the components
    if record.use_sensor && manually_check_absorption
        
        if (sensor_mask_index(2) - sensor_mask_index(1)) > max(kgrid.Nx, kgrid.Ny)
        
            % y-direction source
            uy_p = real(ifftn(uy_p_k_sgy));
            ux_s = real(ifftn(ux_s_k_sgx));
            absorption_check_p(:,t_index) = uy_p(sensor_mask_index);
            absorption_check_s(:,t_index) = ux_s(sensor_mask_index);
            
        else
            
            % x-direction source
            ux_p = real(ifftn(ux_p_k_sgx));
            uy_s = real(ifftn(uy_s_k_sgy));
            absorption_check_p(:,t_index) = ux_p(sensor_mask_index);
            absorption_check_s(:,t_index) = uy_s(sensor_mask_index);

        end
    end
    
    % calculate the velocity gradients
    duxdx_p_k      = bsxfun(@times, ddx_k_shift_neg, kappa_p .* ux_p_k_sgx);   % no_sg
    duxdx_s_k      = bsxfun(@times, ddx_k_shift_neg, kappa_s .* ux_s_k_sgx);   % no_sg

    duydy_p_k      = bsxfun(@times, ddy_k_shift_neg, kappa_p .* uy_p_k_sgy);   % no_sg
    duydy_s_k      = bsxfun(@times, ddy_k_shift_neg, kappa_s .* uy_s_k_sgy);   % no_sg
    
    duxdy_p_k_sgxy = bsxfun(@times, ddy_k_shift_pos, kappa_p .* ux_p_k_sgx);   % sg_xy
    duxdy_s_k_sgxy = bsxfun(@times, ddy_k_shift_pos, kappa_s .* ux_s_k_sgx);   % sg_xy
    
    duydx_p_k_sgxy = bsxfun(@times, ddx_k_shift_pos, kappa_p .* uy_p_k_sgy);   % sg_xy
    duydx_s_k_sgxy = bsxfun(@times, ddx_k_shift_pos, kappa_s .* uy_s_k_sgy);   % sg_xy
    
    % update the normal components and shear components of stress tensor
    if kelvin_voigt_model
    
        % calculate the gradients of the split stress tensor
        dsxxdx_p_sgx = real(ifftn(bsxfun(@times, ddx_k_shift_pos, fftn(sxx_p)))); 
        dsxxdx_s_sgx = real(ifftn(bsxfun(@times, ddx_k_shift_pos, fftn(sxx_s)))); 

        dsyydy_p_sgy = real(ifftn(bsxfun(@times, ddy_k_shift_pos, fftn(syy_p)))); 
        dsyydy_s_sgy = real(ifftn(bsxfun(@times, ddy_k_shift_pos, fftn(syy_s)))); 

        dsxydx_p_sgy = real(ifftn(bsxfun(@times, ddx_k_shift_neg, fftn(sxy_p_sgxy))));
        dsxydx_s_sgy = real(ifftn(bsxfun(@times, ddx_k_shift_neg, fftn(sxy_s_sgxy))));

        dsxydy_p_sgx = real(ifftn(bsxfun(@times, ddy_k_shift_neg, fftn(sxy_p_sgxy))));
        dsxydy_s_sgx = real(ifftn(bsxfun(@times, ddy_k_shift_neg, fftn(sxy_s_sgxy))));
    
        % calculate fractional Laplacian power law absorption terms
        Lx_p_sgx = fftn(Lp1_x_sgx .* real(ifftn(Lp1_k .* ux_p_k_sgx)) + Lp2_x_sgx .* real(ifftn(Lp2_k .* fftn((dsxxdx_p_sgx + dsxydy_p_sgx) ./ rho0_sgx ))));
        Lx_s_sgx = fftn(Ls1_x_sgx .* real(ifftn(Ls1_k .* ux_s_k_sgx)) + Ls2_x_sgx .* real(ifftn(Ls2_k .* fftn((dsxxdx_s_sgx + dsxydy_s_sgx) ./ rho0_sgx ))));
        
        Ly_p_sgy = fftn(Lp1_x_sgy .* real(ifftn(Lp1_k .* uy_p_k_sgy)) + Lp2_x_sgy .* real(ifftn(Lp2_k .* fftn((dsyydy_p_sgy + dsxydx_p_sgy) ./ rho0_sgy ))));
        Ly_s_sgy = fftn(Ls1_x_sgy .* real(ifftn(Ls1_k .* uy_s_k_sgy)) + Ls2_x_sgy .* real(ifftn(Ls2_k .* fftn((dsyydy_s_sgy + dsxydx_s_sgy) ./ rho0_sgy ))));
        
        % calculate derivatives of fractional Laplacian terms
        dLxdx_p = bsxfun(@times, ddx_k_shift_neg, Lx_p_sgx); 
        dLxdx_s = bsxfun(@times, ddx_k_shift_neg, Lx_s_sgx);
        
        dLydy_p = bsxfun(@times, ddy_k_shift_neg, Ly_p_sgy); 
        dLydy_s = bsxfun(@times, ddy_k_shift_neg, Ly_s_sgy);

        dLxdy_p_sgxy = bsxfun(@times, ddy_k_shift_pos, Lx_p_sgx);
        dLxdy_s_sgxy = bsxfun(@times, ddy_k_shift_pos, Lx_s_sgx);
        
        dLydx_p_sgxy = bsxfun(@times, ddx_k_shift_pos, Ly_p_sgy); 
        dLydx_s_sgxy = bsxfun(@times, ddx_k_shift_pos, Ly_s_sgy);
    
        % update the normal components and shear components of stress tensor
        sxx_p = sxx_p + ...
            lambda_dt .* real(ifftn(   duxdx_p_k + duydy_p_k )) ...
            + mu_dt   .* real(ifftn( 2*duxdx_p_k )) ...
            + chi_dt  .* real(ifftn(   dLxdx_p   + dLydy_p )) ...
            + eta_dt  .* real(ifftn( 2*dLxdx_p ));

        sxx_s = sxx_s + ...
            lambda_dt .* real(ifftn(   duxdx_s_k + duydy_s_k )) ...
            + mu_dt   .* real(ifftn( 2*duxdx_s_k )) ...
            + chi_dt  .* real(ifftn(   dLxdx_s   + dLydy_s )) ...
            + eta_dt  .* real(ifftn( 2*dLxdx_s ));    

        syy_p = syy_p + ...
            lambda_dt .* real(ifftn(   duydy_p_k + duxdx_p_k )) ...
            + mu_dt   .* real(ifftn( 2*duydy_p_k)) ...
            + chi_dt  .* real(ifftn(   dLydy_p   + dLxdx_p )) ...
            + eta_dt  .* real(ifftn( 2*dLydy_p  ));

        syy_s = syy_s + ...
            lambda_dt .* real(ifftn(   duydy_s_k + duxdx_s_k )) ...
            + mu_dt   .* real(ifftn( 2*duydy_s_k )) ...
            + chi_dt  .* real(ifftn(   dLydy_s   + dLxdx_s )) ...
            + eta_dt  .* real(ifftn( 2*dLydy_s ));    

        sxy_p_sgxy = sxy_p_sgxy + ...
            mu_dt_sgxy .* real(ifftn( duydx_p_k_sgxy + duxdy_p_k_sgxy )) ...
            + eta_sgxy .* real(ifftn( dLydx_p_sgxy   + dLxdy_p_sgxy ));

        sxy_s_sgxy = sxy_s_sgxy + ...
            mu_dt_sgxy .* real(ifftn( duydx_s_k_sgxy + duxdy_s_k_sgxy )) ...
            + eta_sgxy .* real(ifftn( dLydx_s_sgxy   + dLxdy_s_sgxy ));    
        
    else
        
        % update the normal components and shear components of stress tensor
        sxx_p = sxx_p + ...
            lambda_dt .* real(ifftn(   duxdx_p_k + duydy_p_k )) ...
            + mu_dt   .* real(ifftn( 2*duxdx_p_k ));

        sxx_s = sxx_s + ...
            lambda_dt .* real(ifftn(   duxdx_s_k + duydy_s_k )) ...
            + mu_dt   .* real(ifftn( 2*duxdx_s_k ));    

        syy_p = syy_p + ...
            lambda_dt .* real(ifftn(   duydy_p_k + duxdx_p_k )) ...
            + mu_dt   .* real(ifftn( 2*duydy_p_k ));

        syy_s = syy_s + ...
            lambda_dt .* real(ifftn( duydy_s_k + duxdx_s_k )) ...
            + mu_dt   .* real(ifftn( 2*duydy_s_k ));    

        sxy_p_sgxy = sxy_p_sgxy + mu_dt_sgxy .* real(ifftn( duydx_p_k_sgxy + duxdy_p_k_sgxy ));

        sxy_s_sgxy = sxy_s_sgxy + mu_dt_sgxy .* real(ifftn( duydx_s_k_sgxy + duxdy_s_k_sgxy ));
        
    end
        
    
    % add in the pre-scaled stress source terms
    if sxx_source >= t_index
        if strcmp(source.s_mode, 'dirichlet')
            % enforce the source values as a dirichlet boundary condition
            sxx_p(s_source_pos_index) = source.sxx(:, t_index);
            sxx_s(s_source_pos_index) = source.sxx(:, t_index);
        else
            % add the source values to the existing field values 
            sxx_p(s_source_pos_index) = sxx_p(s_source_pos_index) + source.sxx(:, t_index);
            sxx_s(s_source_pos_index) = sxx_s(s_source_pos_index) + source.sxx(:, t_index);
        end
    end
    if syy_source >= t_index
        if strcmp(source.s_mode, 'dirichlet')
            % enforce the source values as a dirichlet boundary condition
            syy_p(s_source_pos_index) = source.syy(:, t_index);
            syy_s(s_source_pos_index) = source.syy(:, t_index);
        else
            % add the source values to the existing field values 
            syy_p(s_source_pos_index) = syy_p(s_source_pos_index) + source.syy(:, t_index);
            syy_s(s_source_pos_index) = syy_s(s_source_pos_index) + source.syy(:, t_index);
        end
    end
    if sxy_source >= t_index
        if strcmp(source.s_mode, 'dirichlet')
            % enforce the source values as a dirichlet boundary condition
            sxy_p_sgxy(s_source_pos_index) = source.sxy(:, t_index);
            sxy_s_sgxy(s_source_pos_index) = source.sxy(:, t_index);
        else
            % add the source values to the existing field values 
            sxy_p_sgxy(s_source_pos_index) = sxy_p_sgxy(s_source_pos_index) + source.sxy(:, t_index);
            sxy_s_sgxy(s_source_pos_index) = sxy_s_sgxy(s_source_pos_index) + source.sxy(:, t_index);
        end
    end    
    
    % compute pressue from normal components of the stress (used for
    % plotting and output)
    p   = -(sxx_p + sxx_s + syy_p + syy_s)./2;
    
    % extract required sensor data from the pressure and particle velocity
    % fields if the number of time steps elapsed is greater than
    % sensor.record_start_index (defaults to 1) 
    if record.use_sensor && t_index >= sensor.record_start_index
    
        % update index for data storage
        file_index = t_index - sensor.record_start_index + 1;
        
        % run sub-function to extract the required data from the acoustic
        % variables
        sensor_data = kspaceFirstOrder_extractSensorData(2, sensor_data, file_index, sensor_mask_index, record, p, ux_sgx, uy_sgy, []);
        
    end
    
    % estimate the time to run the simulation
    if t_index == ESTIMATE_SIM_TIME_STEPS
  
        % display estimated simulation time
        disp(['  estimated simulation time ' scaleTime(etime(clock, loop_start_time)*index_end/t_index) '...']);

        % check memory usage
        kspaceFirstOrder_checkMemoryUsage; 
                
    end    
    
    % plot data if required
    if plot_sim && (rem(t_index, plot_freq) == 0 || t_index == 1 || t_index == index_end) 

        % update progress bar
        waitbar(t_index/kgrid.Nt, pbar);
        drawnow;   

        % ensure p is cast as a CPU variable and remove the PML from the
        % plot if required
        if strcmp(data_cast, 'gpuArray')
            sii_plot = double(gather(p(x1:x2, y1:y2)));
            sij_plot = double(gather(sxy_p_sgxy(x1:x2, y1:y2) + sxy_s_sgxy(x1:x2, y1:y2)));
        else
            sii_plot = double(p(x1:x2, y1:y2));
            sij_plot = double(sxy_p_sgxy(x1:x2, y1:y2) + sxy_s_sgxy(x1:x2, y1:y2));
        end
        
        % update plot scale if set to automatic or log
        if plot_scale_auto || plot_scale_log
            kspaceFirstOrder_adjustPlotScale;
        end 
        
        % add display mask onto plot
        if strcmp(display_mask, 'default')
            sii_plot(double(sensor.mask(x1:x2, y1:y2)) == 1) = plot_scale(2);
            sij_plot(double(sensor.mask(x1:x2, y1:y2)) == 1) = plot_scale(end); 
        elseif ~strcmp(display_mask, 'off')
            sii_plot(display_mask(x1:x2, y1:y2) ~= 0) = plot_scale(2);
            sij_plot(display_mask(x1:x2, y1:y2) ~= 0) = plot_scale(end);
        end
        
        % update plot
        subplot(1, 2, 1);
        imagesc(kgrid.y_vec(y1:y2)*scale, kgrid.x_vec(x1:x2)*scale, sii_plot, plot_scale(1:2));
        colormap(COLOR_MAP);
        ylabel(['x-position [' prefix 'm]']);
        xlabel(['y-position [' prefix 'm]']);
        title('Normal Stress (\sigma_{ii}/2)')
        axis image;
        
        subplot(1, 2, 2);
        imagesc(kgrid.y_vec(y1:y2)*scale, kgrid.x_vec(x1:x2)*scale, sij_plot, plot_scale(end-1:end));
        colormap(COLOR_MAP);
        ylabel(['x-position [' prefix 'm]']);
        xlabel(['y-position [' prefix 'm]']);
        title('Shear Stress (\sigma_{xy})')
        axis image;        
        
        % force plot update
        drawnow;        

        % save movie frame if required
        if record_movie

            % set background color to white
            set(gcf, 'Color', [1 1 1]);

            % save the movie frame
            writeVideo(video_obj, getframe(gcf));
            
        end        
        
        % update variable used for timing variable to exclude the first
        % time step if plotting is enabled
        if t_index == 1
            loop_start_time = clock;
        end
    end
end

% update command line status
disp(['  simulation completed in ' scaleTime(toc)]);

% =========================================================================
% CLEAN UP
% =========================================================================

% save the movie frames to disk
if record_movie
    close(video_obj);
end

% clean up used figures
if plot_sim
    close(img);
    close(pbar);
end

% save the final pressure field if in time reversal mode
if record.time_rev
    record.p_final = true;
end

% save the final acoustic pressure if required
if record.p_final
    sensor_data.p_final = p(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside);
end

% save the final particle velocity if required
if record.u_final
    sensor_data.ux_final = ux_sgx(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside);
    sensor_data.uy_final = uy_sgy(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside);
end

% run subscript to cast variables back to double precision if required
if data_recast
    kspaceFirstOrder_dataRecast;
end

% reorder the sensor points if a binary sensor mask was used for Cartesian
% sensor mask nearest neighbour interpolation (this is performed after
% recasting as the GPU toolboxes do not all support this subscript)
if record.use_sensor && record.reorder_data
    kspaceFirstOrder_reorderCartData;
end

% filter the recorded time domain pressure signals if transducer filter
% parameters are given 
if record.use_sensor && ~record.time_rev && isfield(sensor, 'frequency_response')
    sensor_data.p = gaussianFilter(sensor_data.p, 1/kgrid.dt, sensor.frequency_response(1), sensor.frequency_response(2));
end

if ~record.use_sensor
    % if sensor is not used, return empty sensor data
    sensor_data = [];
elseif record.time_rev
    % if computing time reversal, reassign sensor_data.p_final to
    % sensor_data
    sensor_data = sensor_data.p_final;
elseif ~isfield(sensor, 'record')
    % if sensor.record is not given by the user, reassign sensor_data.p to
    % sensor_data
    sensor_data = sensor_data.p;
end

% if manually checking absorption, add additional output data
if manually_check_absorption
    sensor_data.absorption_check_p = absorption_check_p;
    sensor_data.absorption_check_s = absorption_check_s;
end

% update command line status
disp(['  total computation time ' scaleTime(etime(clock, start_time))]);

% switch off log
if create_log
    diary off;
end