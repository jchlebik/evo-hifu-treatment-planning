% Script to compare the C++ and MATLAB codes for a range of different input
% and output sizes.
%
% author: Bradley Treeby and Jakub Budisky
% date: 28th March 2017
% last update: 11th June 2017

clearvars;

% set options
plot_diff = true;

% set the test pass
amplitude_max_error = 1e-4;
phase_max_error = 0.05;

% input grid size
Nx = 128;       % [grid points]
Ny = 96;        % [grid points]
Nz = 50;        % [grid points]
dx = 0.1e-3;    % [m]

% expanded grid size
Nx_ex = [432 433 433 433 433 432 432 432];
Ny_ex = [384 384 385 384 385 385 385 384];
Nz_ex = [360 360 360 361 361 360 361 361];

% define the properties of the propagation medium
c0 = 1510;       % [m/s] 

% set the source frequency
f0 = 2e6;        % [Hz]

% set the amplitude
width = 15;
amp_in = zeros(Nx, Ny, Nz);
amp_in(1, Ny/2 - width:Ny/2 + width, Nz/2 - width:Nz/2 + width) = 1;

% set pass variable
all_tests_passed = true;

% loop through the propagators
for prop_index = 1:3
    
    % set the phase
    phase_in = zeros(Nx, Ny, Nz);
    switch prop_index
        case 1
            phase_in(amp_in == 1) = 0;
        case 2
            phase_in(amp_in == 1) = pi/2;
        case 3
            phase_in(amp_in == 1) = pi/4;
    end
            
    % loop through expanded sizes
    for index = 1:length(Nx_ex)

        % set the expanded grid size
        sz_ex = [Nx_ex(index), Ny_ex(index), Nz_ex(index)];

        % calculate beam pattern
        [amp_out, phase_out] = fastBeamPattern(amp_in, phase_in, dx, f0, c0, 'ExpandedGridSize', sz_ex);

        % run C++ code
        [amp_out_cpp, phase_out_cpp] = fastBeamPatternC(amp_in, phase_in, dx, f0, c0, 'ExpandedGridSize', sz_ex);

        % compute the am,plitude error
        L_inf_amp = max(abs(amp_out_cpp(:) - amp_out(:))) / max(abs(amp_out(:))); 

        % compute the phase error accounting for phase wrapping
        angle_diff = @(ref, cmp) (mod(ref - cmp + pi, 2*pi) - pi);
        L_inf_phase = max(abs(angle_diff(phase_out(:), phase_out_cpp(:)))) / (2*pi);

        % print errors
        disp(['Expanded grid size = ' num2str(sz_ex)]);
        disp(['Amplitude error    = ' num2str(L_inf_amp)]);
        disp(['Phase error        = ' num2str(L_inf_phase)]);

        if (L_inf_amp < amplitude_max_error) && (L_inf_phase < phase_max_error)
            disp('Test Passed');
        else
            disp('Test Failed');
            all_tests_passed = false;
        end

        if plot_diff

            % extract a plane
            xy_amp_matlab = squeeze(amp_out(:, :, Nz/2));
            xy_amp_cpp = squeeze(amp_out_cpp(:, :, Nz/2));
            xy_phase_matlab = squeeze(phase_out(:, :, Nz/2));
            xy_phase_cpp = squeeze(phase_out_cpp(:, :, Nz/2));

            % plot the field
            figure;
            subplot(2, 3, 1);
            imagesc(xy_amp_matlab);
            axis image;
            colorbar;
            title('MATLAB');

            subplot(2, 3, 2);
            imagesc(xy_amp_cpp);
            axis image;
            colorbar;
            title('C++');

            subplot(2, 3, 3);
            imagesc(abs(xy_amp_matlab - xy_amp_cpp) ./ max(abs(amp_out(:))));
            axis image;
            colorbar;
            title('Error');

            % plot the phase
            subplot(2, 3, 4);
            imagesc(xy_phase_matlab);
            axis image;
            colorbar;
            title('MATLAB');

            subplot(2, 3, 5);
            imagesc(xy_phase_cpp);
            axis image;
            colorbar;
            title('C++');

            subplot(2, 3, 6);
            imagesc(abs(xy_phase_matlab - xy_phase_cpp) ./ (2*pi));
            axis image;
            colorbar;
            title('Error');

            drawnow;

        end

    end
    
end

% display status
if all_tests_passed
    disp('ALL TESTS PASSED');
else
    disp('TEST FAIL');
end