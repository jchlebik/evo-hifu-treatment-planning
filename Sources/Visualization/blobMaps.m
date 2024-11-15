%Set domain sizes
Nx = 495;
Ny = Nx;
dx = 0.1e-3;
dy = dx;

load('austin_woman_cmaes');

% medium properties
medium.diffusion_coeff = 0.14e-6;
medium.boundary_condition = 'periodic';

% k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

%Target map init
targetMap = zeros(Nx, Ny);
targetMap = targetMap + makeDisc(Nx, Ny, Nx/2+25, Ny/2, 15);
targetMap = targetMap + makeDisc(Nx, Ny, Nx/2-25, Ny/2, 15);
targetMap = targetMap + makeDisc(Nx, Ny, Nx/2, Ny/2+25, 15);
targetMap = targetMap + makeDisc(Nx, Ny, Nx/2, Ny/2-25, 15);
targetMap = targetMap + makeDisc(Nx, Ny, Nx/2, Ny/2, 25);

targetMap = targetMap + makeDisc(Nx, Ny, Nx/2+15, Ny/2+15, 15);
targetMap = targetMap + makeDisc(Nx, Ny, Nx/2+15, Ny/2-15, 15);
targetMap = targetMap + makeDisc(Nx, Ny, Nx/2-15, Ny/2+15, 15);
targetMap = targetMap + makeDisc(Nx, Ny, Nx/2-15, Ny/2-15, 15);

targetMap = (targetMap >= 1);

hole = (makeDisc(Nx, Ny, Nx/2, Ny/2, 10));
targetMap = targetMap - hole;

%Penalize map
penalizeMap = zeros(Nx, Ny);
%penalizeMap(15:Nx-15, 15:Ny-15) = 0;
%penalizeMap = (makeDisc(Nx, Ny, Nx/2, Ny/2, 60) == 0);
%penalizeMap = penalizeMap + makeDisc(Nx, Ny, Nx/2, Ny/2, 10);
penalizeMap = (makeDisc(Nx, Ny, Nx/2, Ny/2, 100));
penalizeMap = penalizeMap - makeDisc(Nx, Ny, Nx/2, Ny/2, 60)
penalizeMap = 1*penalizeMap + 3*hole;

subplot(1, 2, 1); 
imagesc(double(segments) + (targetMap*1000));
axis image;
colorbar;
title('Cíl');  

subplot(1, 2, 2); 
imagesc(double(segments) + (penalizeMap*1000));
axis image;
colorbar;
title('Penalizační maska');

%subplot(1, 3, 3);
%imagesc(double(segments))
%axis image;
%colorbar;
%title('Medium');  
