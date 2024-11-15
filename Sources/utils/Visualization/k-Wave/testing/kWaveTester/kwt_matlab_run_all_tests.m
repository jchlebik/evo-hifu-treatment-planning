% Script to call kWavetester to test the functionality of
% kspaceFirstOrder1D, kspaceFirstOrder2D, and kspaceFirstOrder3D.
%
% author: Bradley Treeby
% date: 23rd February 2017
% last update: 23rd February 2017

clearvars;

% load defaults
kwt_load_defaults;

% no comparison with C++ code
options.run_cpp_comparison_tests = false;

% double precision
options.data_cast = 'off'; 
kwt_loop_dim;

% single precision
options.data_cast = 'single'; 
kwt_loop_dim;

% double precision
options.data_cast = 'gpuArray-single'; 
kwt_loop_dim;

% single precision
options.data_cast = 'gpuArray-double'; 
kwt_loop_dim;