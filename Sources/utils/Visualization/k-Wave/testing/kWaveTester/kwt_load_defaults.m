% Script to set the default options for kWaveTester. This will run all test
% options (including C++ comparisons) in 3D, using the C++ binary in the
% default location.
%
% author: Bradley Treeby
% date: 23rd February 2017
% last update: 21st August 2017

% general options
options.use_nonuniform_grid                     = false;
options.force_plot_off                          = false;
options.output_folder                           = '';
options.save_test_log                           = true;
options.split_log_after_n_tests                 = 50;
options.continue_past_errors                    = true;

% C++ options
options.run_cpp_comparison_tests                = true;
options.plot_cpp_errors                         = true;
options.save_cpp_comparison_plots_to_disk       = true;
options.image_foldername                        = '';
options.use_gpu_code                            = false;
options.cuda_binary_path                        = getkWavePath('binaries');
options.cpp_binary_path                         = getkWavePath('binaries');
options.cpp_save_to_disk_only                   = false;

% name of cpp binary to test
if isunix
    options.cpp_binary_name                     = 'kspaceFirstOrder3D-OMP';
    options.cuda_binary_name                    = 'kspaceFirstOrder3D-CUDA';
else
    options.cpp_binary_name                     = 'kspaceFirstOrder3D-OMP.exe';
    options.cuda_binary_name                    = 'kspaceFirstOrder3D-CUDA.exe';
end

% test options
options.run_source_tests                        = true;
options.run_bin_sensor_tests                    = true;
options.run_cuboid_corner_tests                 = true;
options.run_cart_sensor_lin_interp_tests        = true;
options.run_cart_sensor_nearest_interp_tests    = true;
options.run_display_tests                       = true;
options.run_time_reversal_tests                 = true;

% grid options
options.test_dim                                = 3;
options.test_type                               = 1;
options.data_cast                               = 'off';  
options.test_case_start_index                   = 1;

% custom test
options.custom_test                             = false;