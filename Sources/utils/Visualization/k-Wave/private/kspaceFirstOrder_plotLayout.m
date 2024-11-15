% DESCRIPTION:
%     Subscript for the first-order k-Wave simulation functions to plot the
%     simulation layout, including the initial pressure distribution,
%     source and sensor mask, and material properties. In 3D, only the
%     central planes are displayed.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 7th February 2012
%     last update - 8th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2017 Bradley Treeby

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

% set the colormap to use in 2D and 3D
plot_layout_cmap = flipud(bone(256));

% pre-compute suitable axes scaling factor
[~, scale, prefix] = scaleSI(max([kgrid.x_vec; kgrid.y_vec; kgrid.z_vec]));

% create plot axes
x_vec = kgrid.x_vec(x1:x2) * scale;
if kgrid.dim > 1
    y_vec = kgrid.y_vec(y1:y2) * scale;
end
if kgrid.dim > 2
    z_vec = kgrid.z_vec(z1:z2) * scale;
end

switch kgrid.dim
    case 1
        
        % set the spacing as a fraction of the plot range
        plot_spc = 0.05;
                
        % create figure window
        figure;
        
        % plot the initial pressure distribution or the source mask
        subplot(2, 2, 1);
        if isfield(source, 'p0')
            plot(x_vec, double(source.p0(x1:x2)));
            axis tight;
            plot_layout_lim = get(gca, 'YLim');
            plot_rng = plot_layout_lim(2) - plot_layout_lim(1);
            set(gca, 'YLim', [plot_layout_lim(1) - plot_spc * plot_rng, plot_layout_lim(2) + plot_spc * plot_rng]);
            title('Initial Pressure');
        elseif p_source
            plot(x_vec, double(source.p_mask(x1:x2)));
            axis tight;
            plot_layout_lim = get(gca, 'YLim');
            plot_rng = plot_layout_lim(2) - plot_layout_lim(1);
            set(gca, 'YLim', [plot_layout_lim(1) - plot_spc * plot_rng, plot_layout_lim(2) + plot_spc * plot_rng]);            
            title('Source Mask');
        elseif u_source
            plot(x_vec, double(source.u_mask(x1:x2)));
            axis tight;
            plot_layout_lim = get(gca, 'YLim');
            plot_rng = plot_layout_lim(2) - plot_layout_lim(1);
            set(gca, 'YLim', [plot_layout_lim(1) - plot_spc * plot_rng, plot_layout_lim(2) + plot_spc * plot_rng]);            
            title('Source Mask');            
        end
        
        % plot the sensor mask
        if record.use_sensor
            subplot(2, 2, 2);
            bar(x_vec, double(sensor.mask(x1:x2)), 'b');
            axis tight;
            plot_layout_lim = get(gca, 'YLim');
            plot_rng = plot_layout_lim(2) - plot_layout_lim(1);
            set(gca, 'YLim', [plot_layout_lim(1) - plot_spc * plot_rng, plot_layout_lim(2) + plot_spc * plot_rng]);
            title('Sensor Mask');
        end
        
        % plot the sound speed distribution
        subplot(2, 2, 3);
        if numel(c) == 1
            plot(x_vec, double(c));
        else
            plot(x_vec, double(c(x1:x2)));
        end
        axis tight;
        plot_layout_lim = get(gca, 'YLim');
        plot_rng = plot_layout_lim(2) - plot_layout_lim(1);
        set(gca, 'YLim', [plot_layout_lim(1) - plot_spc * plot_rng, plot_layout_lim(2) + plot_spc * plot_rng]);
        title('Sound Speed');
        
        % plot the density distribution
        subplot(2, 2, 4);
        if numel(rho0) == 1
            plot(x_vec, double(rho0));
        else
            plot(x_vec, double(rho0(x1:x2)));
        end
        axis tight;
        plot_layout_lim = get(gca, 'YLim');
        plot_rng = plot_layout_lim(2) - plot_layout_lim(1);
        set(gca, 'YLim', [plot_layout_lim(1) - plot_spc * plot_rng, plot_layout_lim(2) + plot_spc * plot_rng]);
        title('Density');
        
        % add axis label
        xlabel(['(All horizontal axes in ' prefix 'm)']);
               
    case 2
                
        % create figure window
        figure;        

        % plot the initial pressure distribution or the source mask
        if p0_source || ux_source || uy_source
            subplot(2, 2, 1);
            if p0_source
                imagesc(kgrid.y_vec(y1:y2) * scale, x_vec, double(source.p0(x1:x2, y1:y2)));
                title('Initial Pressure');
            elseif p_source
                imagesc(kgrid.y_vec(y1:y2) * scale, x_vec, double(source.p_mask(x1:x2, y1:y2)));
                title('Source Mask');
            elseif ux_source || uy_source
                imagesc(kgrid.y_vec(y1:y2) * scale, x_vec, double(source.u_mask(x1:x2, y1:y2)));
                title('Source Mask');
            end
            axis image;
            colormap(plot_layout_cmap);    
            h = colorbar;
            if p0_source
                ylabel(h, '[Pa]');
            end
        end

        % plot the sensor mask
        if record.use_sensor
            subplot(2, 2, 2);
            imagesc(kgrid.y_vec(y1:y2) * scale, x_vec, double(sensor.mask(x1:x2, y1:y2)));
            axis image;
            colormap(plot_layout_cmap);
            title('Sensor Mask');
            colorbar;
        end

        % plot the sound speed distribution
        if numel(c) == 1
            subplot(2, 2, 3);
            imagesc(kgrid.y_vec(y1:y2) * scale, x_vec, ...
                double(c) * ones(size(p(x1:x2, y1:y2))));
        else
            subplot(2, 2, 3);
            imagesc(kgrid.y_vec(y1:y2) * scale, x_vec, ...
                double(c(x1:x2, y1:y2)));
        end
        axis image;
        colormap(plot_layout_cmap);
        title('Sound Speed');
        h = colorbar;
        ylabel(h, '[m/s]');

        % plot the density distribution
        if numel(rho0) == 1
            subplot(2, 2, 4);
            imagesc(kgrid.y_vec(y1:y2) * scale, x_vec, ...
                double(rho0) * ones(size(p(x1:x2, y1:y2))));
        else
            subplot(2, 2, 4);
            imagesc(kgrid.y_vec(y1:y2) * scale, x_vec, ...
                double(rho0(x1:x2, y1:y2)));
        end
        axis image;
        colormap(plot_layout_cmap);
        title('Density');
        h = colorbar;
        ylabel(h, '[kg/m^3]');
        
        % add axis label
        xlabel(['(All axes in ' prefix 'm)']);        
    
    case 3
        
        % plot the initial pressure distribution
        if isfield(source, 'p0')
            
            figure;
            plot_layout_data = double(source.p0(x1:x2, y1:y2, z1:z2));
            plot_layout_lim = double([min(plot_layout_data(:)), max(plot_layout_data(:))]);            
            
            subplot(2, 2, 1);
            imagesc(y_vec, x_vec, squeeze(plot_layout_data(:, :, round(end/2))), plot_layout_lim);
            title('Initial Pressure: x-y plane');
            axis image;
            h = colorbar;
            ylabel(h, '[Pa]');

            subplot(2, 2, 2);
            imagesc(z_vec, x_vec, squeeze(plot_layout_data(:, round(end/2), :)), plot_layout_lim);
            title('x-z plane');
            axis image;
            xlabel(['(All axes in ' prefix 'm)']);
            h = colorbar;
            ylabel(h, '[Pa]');

            subplot(2, 2, 3);
            imagesc(z_vec, y_vec, squeeze(plot_layout_data(round(end/2), :, :)), plot_layout_lim);
            title('y-z plane');
            axis image;
            h = colorbar;
            ylabel(h, '[Pa]');
            colormap(plot_layout_cmap);
            
        end

        % plot c if heterogeneous
        if numDim(c) == 3
            
            figure;
            plot_layout_data = double(c(x1:x2, y1:y2, z1:z2));
            plot_layout_lim = double([min(plot_layout_data(:)), max(plot_layout_data(:))]);
            
            subplot(2, 2, 1);
            imagesc(y_vec, x_vec, squeeze(plot_layout_data(:, :, round(end/2))), plot_layout_lim);
            title('Sound Speed: x-y plane');
            axis image;
            h = colorbar;
            ylabel(h, '[m/s]');

            subplot(2, 2, 2);
            imagesc(z_vec, x_vec, squeeze(plot_layout_data(:, round(end/2), :)), plot_layout_lim);
            title('x-z plane');
            axis image;
            xlabel(['(All axes in ' prefix 'm)']);
            h = colorbar;
            ylabel(h, '[m/s]');

            subplot(2, 2, 3);
            imagesc(z_vec, y_vec, squeeze(plot_layout_data(round(end/2), :, :)), plot_layout_lim);
            title('y-z plane');
            axis image;
            h = colorbar;
            ylabel(h, '[m/s]');
            colormap(plot_layout_cmap);
            
        end

        % plot rho0 if heterogeneous
        if numDim(rho0) == 3    
            
            figure;
            plot_layout_data = double(rho0(x1:x2, y1:y2, z1:z2));
            plot_layout_lim = double([min(plot_layout_data(:)), max(plot_layout_data(:))]);
            
            subplot(2, 2, 1);
            imagesc(y_vec, x_vec, squeeze(plot_layout_data(:, :, round(end/2))), plot_layout_lim);
            title('Density: x-y plane');
            axis image;
            h = colorbar;
            ylabel(h, '[kg/m^3]');

            subplot(2, 2, 2);
            imagesc(z_vec, x_vec, squeeze(plot_layout_data(:, round(end/2), :)), plot_layout_lim);
            title('x-z plane');
            axis image;
            xlabel(['(All axes in ' prefix 'm)']);
            h = colorbar;
            ylabel(h, '[kg/m^3]');

            subplot(2, 2, 3);
            imagesc(z_vec, y_vec, squeeze(plot_layout_data(round(end/2), :, :)), plot_layout_lim);
            title('y-z plane');
            axis image;
            h = colorbar;
            ylabel(h, '[kg/m^3]');
            colormap(plot_layout_cmap);
            
        end
        
end

% clean up unused variables
clear plot_layout_data plot_layout_lim plot_rng plot_spc plot_layout_cmap;