% Script to call kWavetester to test the functionality of
% kspaceFirstOrder3D-OMP compared with kspaceFirstOrder3D.
%
% author: Bradley Treeby
% date: 23rd February 2017
% last update: 3rd May 2017

clearvars;

% load defaults
kwt_load_defaults;

% no display
options.run_display_tests = false;
options.force_plot_off = true;

% set folder for images
options.image_foldername = ['omp_comparison_tests_' getDateString filesep];

% run tests
options.test_dim = 3;
kwt_loop_test_type;