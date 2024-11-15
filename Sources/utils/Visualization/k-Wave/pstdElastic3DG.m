function sensor_data = pstdElastic3DG(varargin)
%KSPACEFIRSTORDER3DG 3D time-domain simulation of elastic propagation on a GPU using C++ CUDA code. 
%
% DESCRIPTION:
%     pstdElastic3DG provides an blind interface to the native C++ CUDA
%     version of pstdElastic3D (called pstdElastic3D-CUDA). This function
%     is a wrapper for kspaceFirstOrder3DC and replaces the default binary
%     name with the name of the elastic GPU binary. 
%
%     This function requires the C++ binary/executable of
%     pstdElastic3D-CUDA to be downloaded from
%     http://www.k-wave.org/download.php and placed in the "binaries"
%     directory of the k-Wave toolbox. 
% 
%     Note, not all input options are currently supported, and all display
%     options are ignored (only command line outputs are given). See the
%     k-Wave user manual for more information. 
%
% USAGE:
%     see pstdElastic3D and kspaceFirstOrder3DC
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 1st May 2017
%     last update     - 1st May 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017 Bradley Treeby
%
% See also kspaceFirstOrder3D, kspaceFirstOrder3DC, kspaceFirstOrder3DG

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

% set the name of the k-Wave MATLAB function to use to generate the input
% file
kwave_function_name = 'pstdElastic3D';

% check for a custom binary name
if any(strcmp('BinaryName', varargin))
    
    % if the binary name is given, directly pass to kspaceFirstOrder3DC
    sensor_data = kspaceFirstOrder3DC(varargin{:}, ...
        'FunctionName', kwave_function_name);
    
else

    % if the binary name is not given, specify to use the GPU binary
    if isunix
        binary_name = 'pstdElastic3D-CUDA';
    else
        binary_name = 'pstdElastic3D-CUDA.exe';
    end
    
    % pass this and the original inputs to kspaceFirstOrder3DC
    sensor_data = kspaceFirstOrder3DC(varargin{:}, ...
        'BinaryName', binary_name, 'FunctionName', kwave_function_name);
    
end