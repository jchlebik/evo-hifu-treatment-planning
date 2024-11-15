function [amp_out, phase_out, field_t1, field_t2] = fastBeamPattern(amp_in, phase_in, dx, f0, c0, varargin)
%FASTBEAMPATTERN Calculate beam pattern from an arbitrary phased array.
%
% DESCRIPTION:
%     fastBeamPattern calculates the steady state beam pattern (amplitude
%     and phase) from an arbitrary phased array driven by a single
%     frequency continuous wave sinusoid in a homogeneous and lossless
%     medium. The phased array is defined in completely general form as a
%     matrix of amplitude and phase. This allows arrays of arbitrary
%     geometry and numbers of elements to be defined.
%
%     The solution is based on the Green's function for the homogeneous
%     wave equation expressed in the spatial frequency domain or k-space.
%     The temporal convolution integral is solved analytically, and the
%     remaining integrals are expressed in the form of the spatial Fourier
%     transform. This allows the acoustic pressure for all spatial
%     positions at any time t > 0 to be calculated in a single step without
%     numerical quadrature. To avoid wave wrapping, the domain size used
%     for calculation is automatically expanded to a suitable dimension
%     size with small prime factors.
%
% USAGE:
%     [amp_out, phase_out] = fastBeamPattern(amp_in, phase_in, dx, f, c0)
%     [amp_out, phase_out] = fastBeamPattern(amp_in, phase_in, dx, f, c0, ...)
%     [amp_out, phase_out, field_t1] = fastBeamPattern(amp_in, phase_in, dx, f0, c0, ...)
%     [amp_out, phase_out, field_t1, field_t2] = fastBeamPattern(amp_in, phase_in, dx, f0, c0, ...)
%
% INPUTS:
%     amp_in       - matrix of the source amplitude at each grid point [Pa]
%     phase_in     - matrix of the source phase at each grid point [rad]
%     dx           - grid spacing [m]
%     f0           - source frequency [Hz]
%     c0           - medium sound speed [m/s]
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%     
%     'GridExpansionFactor'
%                  - Option to specify the multiplicative factor used to
%                    calculate the minimum expanded grid size to avoid wave
%                    wrapping based on t1 (default = 1.1). Note, setting a
%                    value for the optional input 'ExpandedGridSize' will
%                    override this value.
%     'GridSearchRange'
%                  - Option to set the search range used to find the
%                    expanded grid size with the smallest prime factors
%                    (default = 50). Note, setting a value for the optional
%                    input 'ExpandedGridSize' will override this value.
%     'SaveToDisk' - String containing a filename (including pathname if
%                    required). If set, after the precomputation phase,
%                    the variables used in calculation are saved to the
%                    specified location in HDF5 format. The function then
%                    exits. The saved variables can be used to run
%                    simulations using the C++ code.
%     'Time'       - Option to specify the time t1 at which the wavefield
%                    is calculated to extract the amplitude and phase
%                    (default = 'auto').
%     'TimeExpansionFactor'
%                  - Option to specify the multiplicative factor used to
%                    calculate t1 based on the time taken to propagate
%                    across the longest grid diagonal (default = 1.5).
%                    Note, setting a value for the optional input 'Time'
%                    will override this value. 
%     'ExpandedGridSize' 
%                  - Option to specify the size of the grid after expansion
%                    used to avoid wave wrapping (default = 'auto').
%
% OUTPUTS:
%     amp_out      - matrix of the output amplitude at each grid point in
%                    steady state [Pa] 
%     phase_out    - matrix of the output phase at each grid point in
%                    steady state [rad]
%     field_t1     - pressure field at time t1
%     field_t2     - pressure field at time t2
%
% ABOUT:
%     author       - Bradley Treeby and Ben Cox
%     date         - 21st April 2016
%     last update  - 24th July 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2016-2017 Bradley Treeby and Ben Cox
%
% See also fastBeamPatternC, kspaceFirstOrder2D, kspaceFirstOrder3D

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

%#ok<*UNRCH>

% start the timer
start_time = clock;

% =========================================================================
% DEFINE LITERALS
% =========================================================================

% define general literals
CAST_TO_SINGLE              = false;
NUM_REQUIRED_INPUTS         = 5;
MACHINE_PRECISION           = eps * 10;

% define literals used for 'SaveToDisk' option
FBP_FILE_MAJOR_VERSION      = 1;
FBP_FILE_MINOR_VERSION      = 0;
HDF5_MATRIX_DATA_TYPE       = 'single';
HDF5_INTEGER_DATA_TYPE      = 'uint64';
HDF5_FILE_TYPE              = 'input';

% define attribute names used for 'SaveToDisk' option
FILE_MAJOR_VER_ATT_NAME     = 'major_version';
FILE_MINOR_VER_ATT_NAME     = 'minor_version';
FILE_DESCR_ATT_NAME         = 'file_description';
FILE_CREATION_DATE_ATT_NAME = 'creation_date';
CREATED_BY_ATT_NAME         = 'created_by';
FILE_TYPE_ATT_NAME          = 'file_type';

% define literals used for 'UseFirstOrder' option
FIRST_ORDER_PML_SIZE        = 20;
FIRST_ORDER_CFL             = 0.1;
FIRST_ORDER_DATA_CAST       = 'single';
FIRST_ORDER_PLOT_SIM        = false;
FIRST_ORDER_RAMP_LENGTH     = 0;

% define literals used as defaults
SAVE_TO_DISK_DEF            = false;
T1_DEF                      = 'auto';
EXPANDED_GRID_SZ_DEF        = 'auto';
USE_FIRST_ORDER_DEF         = false;
PLOT_PROPAGATORS_DEF        = false;
TIME_EXP_FACTOR_DEF         = 1.5;
GRID_EXP_FACTOR_DEF         = 1.1;
GRID_SEARCH_RANGE_DEF       = 50;

% =========================================================================
% OPTIONAL INPUT VARIABLES
% =========================================================================

% assign values with defaults
save_to_disk        = SAVE_TO_DISK_DEF;
t1                  = T1_DEF;
sz_ex               = EXPANDED_GRID_SZ_DEF;
use_first_order     = USE_FIRST_ORDER_DEF;
plot_propagators    = PLOT_PROPAGATORS_DEF;
time_exp_factor     = TIME_EXP_FACTOR_DEF;
grid_exp_factor     = GRID_EXP_FACTOR_DEF;
grid_search_range   = GRID_SEARCH_RANGE_DEF;

% replace with user defined values if provided
if nargin < NUM_REQUIRED_INPUTS
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'GridExpansionFactor'
                grid_exp_factor = varargin{input_index + 1};
            case 'GridSearchRange'
                grid_search_range = varargin{input_index + 1};                
            case 'PlotPropagators'
                plot_propagators = varargin{input_index + 1};
            case 'SaveToDisk'
                save_to_disk = varargin{input_index + 1};
            case 'Time'
                t1 = varargin{input_index + 1};
            case 'TimeExpansionFactor'
                time_exp_factor = varargin{input_index + 1};
            case 'ExpandedGridSize'
                sz_ex = varargin{input_index + 1}; 
            case 'UseFirstOrder'
                use_first_order = varargin{input_index + 1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

% =========================================================================
% PRECALCULATIONS
% =========================================================================

% get grid size
sz = size(amp_in);

% update command line status
disp('Calculating beam pattern...');
displayGridSize(sz, 'input grid size');
disp(['  maximum supported frequency: ' scaleSI( c0 / (2 * dx) ) 'Hz']);

% calculate angular frequency [rad/s]
w0 = 2*pi*f0;

% check for maximum supported frequency
if (w0 / c0) > (pi / dx)
    error(['Input frequency is higher than maximum supported frequency of ' scaleSI(c0/(2*dx)) 'Hz.' ]);
end

% calculate the length of time needed for the wave to propagate across the
% grid diagonal (using exactly the length of the grid diagonal results in
% some errors at the wave front, so an additional factor is included)
if strcmp(t1, 'auto')
    t1 = time_exp_factor * dx * sqrt(sum(sz.^2)) / c0;
end

% compute second pressure field offset by a quarter period (this means one
% wave is like a sin, one is like a cos, and the phase can be recovered
% using atan)
t2 = t1 + pi / (2 * w0);

% =========================================================================
% SIMULATION USING KSPACEFIRSTORDER
% =========================================================================

% undocumented option to calculate the beam pattern using 
% kspaceFirstOrderND for comparision
if use_first_order

    % create grid
    switch numDim(amp_in)
        case 2
            kgrid = kWaveGrid(sz(1), dx, sz(2), dx);
        case 3
            kgrid = kWaveGrid(sz(1), dx, sz(2), dx, sz(3), dx);
    end

    % define time step using an integer number of PPP
    PPW = c0 / (f0 * dx);
    PPP = floor(PPW / FIRST_ORDER_CFL);
    dt  = 1 / (f0 * PPP);
    Nt = round(t1 / dt);

    % create the time array, recording for an integer number of PPP if
    % returning amplitude and phase, or ensuring that t1 is divisible by dt
    % (accounting for the temporal staggered grid) if returning the field
    if nargout == 2
        Nt = PPP * round(Nt / PPP);
    else
        dt = t1 / Nt;
        t1 = Nt*dt + dt/2;
        dt = t1 / Nt;
    end
    
    kgrid.t_array = (0:(Nt - 1))*dt;

    % assign medium properties
    medium.sound_speed = c0;

    % create source
    source.p_mask = double(amp_in ~= 0);
    if numel(phase_in) == 1
        source.p = createCWSignals(kgrid.t_array, f0, amp_in(amp_in ~= 0), phase_in .* ones(sum(source.p_mask(:)), 1), FIRST_ORDER_RAMP_LENGTH);
    else
        source.p = createCWSignals(kgrid.t_array, f0, amp_in(amp_in ~= 0), phase_in(amp_in ~= 0), FIRST_ORDER_RAMP_LENGTH);
    end

    % create sensor
    sensor.mask = ones(sz);
    sensor.record = {'p', 'p_final'};

    % only record the last few periods in steady state
    num_periods = 3;
    sensor.record_start_index = Nt - num_periods * PPP + 1;

    % set input options
    input_args = {'PMLInside', false, 'PlotSim', FIRST_ORDER_PLOT_SIM, 'DataCast', FIRST_ORDER_DATA_CAST, 'PMLSize', FIRST_ORDER_PML_SIZE};
    
    % run simulation
    switch numDim(amp_in)
        case 2
            sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
        case 3
            sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
    end

    % extract amplitude and phase from recorded data
    [amp_out, phase_out] = extractAmpPhase(sensor_data.p, 1/kgrid.dt, f0, 'Dim', 2);

    % reshape to the size of the input
    amp_out = reshape(amp_out, size(amp_in));
    phase_out = reshape(phase_out, size(amp_in));

    % assign output pressure field if required
    if nargout == 3
        field_t1 = sensor_data.p_final;
    end
    if nargout == 4
        field_t1 = [];
    end
    
    % exit
    return
    
end

% =========================================================================
% GRID EXPANSION
% =========================================================================

% automatically calculate the grid expansion to avoid wave wrapping
if strcmp(sz_ex, 'auto')
    minimum_expansion = ceil(grid_exp_factor * t2 * c0 / dx);
    sz_ex = getOptimalGridSize(sz + minimum_expansion, grid_search_range);
end

% display expanded grid size
displayGridSize(sz_ex, 'expanded grid size');

% =========================================================================
% SAVE TO DISK
% =========================================================================

if save_to_disk

    % update command line status
    disp('  saving input files to disk...');
    
    % list of all the single precision variables
    variable_list = {'amp_in', 'phase_in', 'w0', 'c0', 't1', 't2', 'dx'};

    % list of all the integer variables
    integer_variable_list = {'sz_ex'}; 
    
    % check if HDF5 file exists, and delete if it does (the hdf5 library
    % will give an error if the file already exists)
    if exist(save_to_disk, 'file')
        delete(save_to_disk);
    end

    % change all the variables to be in single precision (float in C++),
    % then add to HDF5 File
    for cast_index = 1:length(variable_list)

        % cast matrix to single precision
        eval([variable_list{cast_index} ' = ' HDF5_MATRIX_DATA_TYPE '(' variable_list{cast_index} ');']);

        % write to HDF5 file
        writeMatrix(save_to_disk, eval(variable_list{cast_index}), variable_list{cast_index});

    end

    % change all the index variables to be in 64-bit unsigned integers
    % (long in C++), then add to HDF5 file
    for cast_index = 1:length(integer_variable_list)

        % cast matrix to 64-bit unsigned integer
        eval([integer_variable_list{cast_index} ' = ' HDF5_INTEGER_DATA_TYPE '(' integer_variable_list{cast_index} ');']);

        % write to HDF5 file
        writeMatrix(save_to_disk, eval(integer_variable_list{cast_index}), integer_variable_list{cast_index});
        
    end
    
    % get computer infor
    comp_info = getComputerInfo;

    % set file description if not provided by user
    file_description = ['Input data created by ' comp_info.user_name ' running MATLAB ' comp_info.matlab_version ' on ' comp_info.operating_system_type];

    % set additional file attributes
    h5writeatt(save_to_disk, '/', FILE_MAJOR_VER_ATT_NAME, FBP_FILE_MAJOR_VERSION);
    h5writeatt(save_to_disk, '/', FILE_MINOR_VER_ATT_NAME, FBP_FILE_MINOR_VERSION);
    h5writeatt(save_to_disk, '/', CREATED_BY_ATT_NAME, ['k-Wave ' comp_info.kwave_version]);
    h5writeatt(save_to_disk, '/', FILE_DESCR_ATT_NAME, file_description);
    h5writeatt(save_to_disk, '/', FILE_TYPE_ATT_NAME, HDF5_FILE_TYPE);
    h5writeatt(save_to_disk, '/', FILE_CREATION_DATE_ATT_NAME, getDateString);    
    
    % update command line status
    disp(['  computation completed in ' scaleTime(etime(clock, start_time))]);
    
    % stop evaluation if no outputs are required
    if nargout == 0
        return
    end

end

% =========================================================================
% GRID EXPANSION
% =========================================================================

% expand (zero pad) the amplitude matrix
amp_in_ex = zeros(sz_ex);
amp_in_ex(1:size(amp_in, 1), 1:size(amp_in, 2), 1:size(amp_in, 3)) = amp_in;

% expand (zero pad) the phase if given as a matrix, otherwise leave as a
% scalar
if numel(phase_in) > 1
    phase_in_ex = zeros(sz_ex);
    phase_in_ex(1:size(amp_in, 1), 1:size(amp_in, 2), 1:size(amp_in, 3)) = phase_in;
else
    phase_in_ex = phase_in;
end

% define spatial wavenumbers (these are already shifted)
k = getWavenumbers(numDim(amp_in), sz_ex);

% cast all variables to single precision
if CAST_TO_SINGLE
    amp_in_ex   = single(amp_in_ex); 
    phase_in_ex = single(phase_in_ex);
    k           = single(k);
    w0          = single(w0);
    c0          = single(c0);
    t1          = single(t1);
    t2          = single(t2);    
end

% =========================================================================
% CALCULATE WAVE FIELDS AND AMPLITUDE AND PHASE
% =========================================================================

% calculate forward FFTs
A_sin = fftn(amp_in_ex .* sin(phase_in_ex));
A_cos = fftn(amp_in_ex .* cos(phase_in_ex));

% compute pressure field at time t1
pressure_field1 = computePressureField(t1);

% compute pressure field at time t2
pressure_field2 = computePressureField(t2);

% compute amplitude and phase from the two beam patterns
amp_out   = sqrt(pressure_field1.^2 + pressure_field2.^2);
phase_out = atan2(pressure_field1, pressure_field2);

% trim matrices back to their original size    
amp_out   = amp_out  (1:size(amp_in, 1), 1:size(amp_in, 2), 1:size(amp_in, 3));
phase_out = phase_out(1:size(amp_in, 1), 1:size(amp_in, 2), 1:size(amp_in, 3));

% apply scaling factor to account for spatial spread of BLI
amp_out = amp_out * 2 * c0 / dx;

% assign the output pressure field if required, including the scaling
% factor, and trim back to the original grid size
if nargout >= 3
    field_t1 = pressure_field1(1:size(amp_in, 1), 1:size(amp_in, 2), 1:size(amp_in, 3)) * 2 * c0 / dx;
end
if nargout == 4
    field_t2 = pressure_field2(1:size(amp_in, 1), 1:size(amp_in, 2), 1:size(amp_in, 3)) * 2 * c0 / dx;
end

% update command line status
disp(['  computation completed in ' scaleTime(etime(clock, start_time))]);

% =========================================================================
% SUBFUNCTION TO CALCULATE PRESSURE FIELD
% =========================================================================

function pressure_field = computePressureField(time)
%COMPUTEPRESSUREFIELD Calculate pressure field at a given time.
%
% DESCRIPTION:
%     computePressureField is a nested function to compute the pressure
%     field at a given time. The remaining variables (e.g., frequency and
%     sound speed) are accessed from the workspace of the parent function.
%
% USAGE:
%     beam_pattern = computePressureField(time)
%
% INPUTS:
%     time            - time [s]
%
% OUTPUTS:
%     pressure_field  - pressure field at specified time [Pa]
%
% ABOUT:
%     author          - Bradley Treeby and Ben Cox
%     date            - 21st April 2016
%     last update     - 3rd April 2017

% define propagators
propagator1 = w0 * (          cos(w0 * time) -      cos(c0 * k * time) ) ./ ( (c0 * k).^2 -          w0.^2 );
propagator2 = w0 * ( c0 * k * sin(w0 * time) - w0 * sin(c0 * k * time) ) ./ ( (c0 * k).^3 - c0 * k * w0.^2 ) ...
    - sin(c0 * k * time) ./ (c0 * k);

% replace values where denominator goes to zero (k == w0/c0) with limit,
% using a small neighbourhood to account for numerical rounding
propagator1(abs(k - w0/c0) < max(k(:)) * MACHINE_PRECISION) = time * sin(w0 * time) / 2;
propagator2(abs(k - w0/c0) < max(k(:)) * MACHINE_PRECISION) = -(sin(w0 * time) + time * w0 * cos(w0 * time) ) / (2 * w0);

% replace values where denominator goes to zero (k == 0) with limit
propagator2(k == 0) = -sin(w0 * time) / w0;

% compute pressure field (amplitude is given in units of pressure not
% density, so the factor of c0^2 given in [1] is not included)
pressure_field = real(ifftn( propagator1 .* A_cos - ...
                             propagator2 .* A_sin ...
                          ));

% plot propagators
if plot_propagators
    
    figure;
    subplot(1, 2, 1);
    imagesc(fftshift(propagator1));
    axis image;
    title('I_1');
    
    subplot(1, 2, 2);
    imagesc(fftshift(propagator2));
    axis image;
    title('I_2');
    scaleFig(2, 1);
    
end
                      
end

% =========================================================================
% SUBFUNCTION TO DISPLAY GRID SIZE
% =========================================================================

function displayGridSize(mat_sz, disp_string)
%DISPLAYGRIDSIZE Display size of input matrix.
%
% USAGE:
%     displayGridSize(mat, disp_string)
%
% INPUTS:
%     mat             - matrix to measure
%     disp_string     - string to display
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 29th January 2017
%     last update     - 27th February 2017

% calculate axis scaling factor
[~, scale, prefix] = scaleSI(min(mat_sz .* dx));

switch length(mat_sz)
    case 2
        disp(['  ' disp_string ': ' num2str(mat_sz(1)) ' by ' num2str(mat_sz(2)) ' grid points (' num2str(scale * mat_sz(1) * dx) ' by ' num2str(scale * mat_sz(2) * dx) prefix 'm)']);
    case 3
        disp(['  ' disp_string ': ' num2str(mat_sz(1)) ' by ' num2str(mat_sz(2)) ' by ' num2str(mat_sz(3)) ' grid points (' num2str(scale * mat_sz(1) * dx) ' by ' num2str(scale * mat_sz(2) * dx) ' by ' num2str(scale * mat_sz(3) * dx) prefix 'm)']); 
end

end

% =========================================================================
% SUBFUNCTION TO CALCULATE EXPANDED GRID SIZE
% =========================================================================

function sz_opt = getOptimalGridSize(sz, test_range)
%GETOPTIMALGRIDSIZE Find grid size with smallest prime factors.
%
% DESCRIPTION:
%     getOptimalGridSize finds the grid size with the smallest prime
%     factors bound by sz on the left, and sz + test_range on the right. If
%     sz has more than one element, the optimal grid size for each
%     dimension is return. 
%
% USAGE:
%     sz_opt = getOptimalGridSize(sz, test_range)
%
% INPUTS:
%     sz              - minimum size
%     test_range      - maximum range of sizes to test
%
% OUTPUTS:
%     sz_opt          - the grid size in the range [sz, sz + test_range]
%                       with the smallest prime factors
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 29th January 2017
%     last update     - 29th January 2017

% extract factors for each dimension
facs = zeros(length(sz), test_range);
for dim = 1:length(sz)
    for index = 1:test_range + 1
        facs(dim, index) = max(factor(sz(dim) + index - 1));
    end
end

% get best dimension size
[~, ind_opt] = min(facs, [], 2);
sz_opt = (sz + ind_opt.' - 1);

end

% =========================================================================
% SUBFUNCTION TO CALCULATE WAVENUMBERS
% =========================================================================

function k_mat = getWavenumbers(ndims, sz)
%GETWAVENUMBERS Calculate scalar wavenumber matrix.
%
% DESCRIPTION:
%     getWavenumbers returns a matrix of the scalar wavenumbers following
%     the definition used in kWaveGrid.
%
% USAGE:
%     k = getWavenumbers(ndims, sz)
%
% INPUTS:
%     ndims           - number of dimensions
%     sz              - matrix size
%
% OUTPUTS:
%     k               - matrix of the scalar wavenumber
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 1st February 2017
%     last update     - 3rd April 2017

% preallocate output matrix
k_mat = zeros(sz);

% loop over the number of dimensions
for index = 1:ndims
    
    % create wavenumbers vector for current dimension
    N = sz(index);
    if rem(N, 2) == 0
        n_vec = ((-N/2:N/2-1)/N);
    else
        n_vec = ((-(N-1)/2:(N-1)/2)/N);
    end

    % force middle value to be zero in case 1/Nx is a recurring
    % number and the series doesn't give exactly zero
    n_vec(floor(N/2) + 1) = 0;

    % scale the wavenumber vector components
    k_vec = (2*pi/dx) .* n_vec;   

    % put in the correct direction for bsxfun
    switch index
        case 1
            k_vec = reshape(k_vec, [], 1, 1);
        case 2
            k_vec = reshape(k_vec, 1, [], 1);
        case 3
            k_vec = reshape(k_vec, 1, 1, []);
    end
    
    % add to the wavenumber matrix
    k_mat = bsxfun(@plus, k_vec.^2, k_mat);
    
end

% take the square root, and shift
k_mat = ifftshift(sqrt(k_mat));

end

% =========================================================================

end