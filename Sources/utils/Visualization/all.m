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
imagesc(double(segments) + (targetMap*-350) + (penalizeMap*200));
axis image;
colorbar;
title('Flower map');  
%% Target area
%Target map init
centerBasic.Nx = Nx/2+44;
centerBasic.Ny = Ny/2+5;


targetMapBasic = zeros(Nx, Ny);
% Basic region
targetMapBasic = targetMapBasic + makeDisc(Nx, Ny, centerBasic.Nx, centerBasic.Ny, 20);
targetMapBasic = targetMapBasic + makeDisc(Nx, Ny, centerBasic.Nx+20, centerBasic.Ny+10, 20);
targetMapBasic = targetMapBasic + makeDisc(Nx, Ny, centerBasic.Nx+30, centerBasic.Ny+20, 20);
targetMapBasic = targetMapBasic + makeDisc(Nx, Ny, centerBasic.Nx, centerBasic.Ny+20, 20);
targetMapBasic = targetMapBasic + makeDisc(Nx, Ny, centerBasic.Nx+20, centerBasic.Ny+20, 20);
targetMapBasic = targetMapBasic + makeDisc(Nx, Ny, centerBasic.Nx+15, centerBasic.Ny-5, 10);
targetMapBasic = targetMapBasic + makeDisc(Nx, Ny, centerBasic.Nx+10, centerBasic.Ny+30, 10);
targetMapBasic = targetMapBasic + makeDisc(Nx, Ny, centerBasic.Nx-10, centerBasic.Ny+10, 10);
targetMapBasic= (targetMapBasic >= 1);

% Critical region
centerCritical.Nx = Nx/2+53;
centerCritical.Ny = Ny/2+15;

targetMapCritical = zeros(Nx, Ny);
targetMapCritical = targetMapCritical + makeDisc(Nx, Ny, centerCritical.Nx, centerCritical.Ny, 15);
targetMapCritical = targetMapCritical + makeDisc(Nx, Ny, centerCritical.Nx+12, centerCritical.Ny+7, 15);
targetMapCritical = targetMapCritical + makeDisc(Nx, Ny, centerCritical.Nx+1, centerCritical.Ny+7, 15);
targetMapCritical = (targetMapCritical >= 1);

targetMap = targetMapBasic + targetMapCritical * 2;

%% Penalize area
%Penalize map, Area 1 - bone
centerPenalizeA1.Nx = Nx/2+60;
centerPenalizeA1.Ny = Ny/2-30;

penalizeMapA1 = zeros(Nx, Ny);
penalizeMapA1 =  penalizeMapA1 + makeDisc(Nx, Ny, centerPenalizeA1.Nx, centerPenalizeA1.Ny, 8);
penalizeMapA1 =  penalizeMapA1 + makeDisc(Nx, Ny, centerPenalizeA1.Nx+10, centerPenalizeA1.Ny+2, 8);
penalizeMapA1 =  penalizeMapA1 + makeDisc(Nx, Ny, centerPenalizeA1.Nx+20, centerPenalizeA1.Ny+4, 8);
penalizeMapA1 =  penalizeMapA1 + makeDisc(Nx, Ny, centerPenalizeA1.Nx+30, centerPenalizeA1.Ny+6, 8);
penalizeMapA1 =  penalizeMapA1 + makeDisc(Nx, Ny, centerPenalizeA1.Nx+40, centerPenalizeA1.Ny+8, 8);

%Penalize map, Area 2 - close to bone
centerPenalizeA2.Nx = Nx/2+60;
centerPenalizeA2.Ny = Ny/2-30;

penalizeMapA2 = zeros(Nx, Ny);
penalizeMapA2 =  penalizeMapA2 + makeDisc(Nx, Ny, centerPenalizeA2.Nx, centerPenalizeA2.Ny, 17);
penalizeMapA2 =  penalizeMapA2 + makeDisc(Nx, Ny, centerPenalizeA2.Nx+15, centerPenalizeA2.Ny+5, 17);
penalizeMapA2 =  penalizeMapA2 + makeDisc(Nx, Ny, centerPenalizeA2.Nx+30, centerPenalizeA2.Ny+10, 17);
penalizeMapA2 =  penalizeMapA2 + makeDisc(Nx, Ny, centerPenalizeA2.Nx+45, centerPenalizeA2.Ny+15, 17);

%Penalize map, Area 3 - soft cover
centerPenalizeA3.Nx = Nx/2+60;
centerPenalizeA3.Ny = Ny/2-30;

penalizeMapA3 = zeros(Nx, Ny);

penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx+45, centerPenalizeA3.Ny+30, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx+55, centerPenalizeA3.Ny+40, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx+55, centerPenalizeA3.Ny+40, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx+55, centerPenalizeA3.Ny+60, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx+50, centerPenalizeA3.Ny+70, 17);

penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx+45, centerPenalizeA3.Ny+82, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx+30, centerPenalizeA3.Ny+100, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx+20, centerPenalizeA3.Ny+100, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx+10, centerPenalizeA3.Ny+100, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx+0, centerPenalizeA3.Ny+100, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-10, centerPenalizeA3.Ny+100, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-20, centerPenalizeA3.Ny+100, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-30, centerPenalizeA3.Ny+100, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-40, centerPenalizeA3.Ny+100, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-50, centerPenalizeA3.Ny+100, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-60, centerPenalizeA3.Ny+100, 17);

penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-60, centerPenalizeA3.Ny+100, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-60, centerPenalizeA3.Ny+90, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-60, centerPenalizeA3.Ny+80, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-60, centerPenalizeA3.Ny+70, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-60, centerPenalizeA3.Ny+60, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-60, centerPenalizeA3.Ny+50, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-60, centerPenalizeA3.Ny+40, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-60, centerPenalizeA3.Ny+30, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-60, centerPenalizeA3.Ny+20, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-60, centerPenalizeA3.Ny+10, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-60, centerPenalizeA3.Ny+0, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-60, centerPenalizeA3.Ny-10, 17);

penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-50, centerPenalizeA3.Ny-10, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-40, centerPenalizeA3.Ny-10, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-30, centerPenalizeA3.Ny-10, 17);
penalizeMapA3 =  penalizeMapA3 + makeDisc(Nx, Ny, centerPenalizeA3.Nx-20, centerPenalizeA3.Ny-5, 17);

%Penalize map, Area 4 - hard cover
centerPenalizeA4.Nx = Nx/2+68;
centerPenalizeA4.Ny = Ny/2-30;

penalizeMapA4 = zeros(Nx, Ny);

penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx+45, centerPenalizeA4.Ny+30, 15);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx+55, centerPenalizeA4.Ny+40, 15);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx+55, centerPenalizeA4.Ny+40, 15);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx+55, centerPenalizeA4.Ny+60, 15);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx+50, centerPenalizeA4.Ny+70, 15);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx+45, centerPenalizeA4.Ny+82, 15);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx+35, centerPenalizeA4.Ny+100, 15);

centerPenalizeA4.Ny = Ny/2-20;
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx+20, centerPenalizeA4.Ny+95, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx+10, centerPenalizeA4.Ny+95, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx+0, centerPenalizeA4.Ny+95, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-10, centerPenalizeA4.Ny+95, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-20, centerPenalizeA4.Ny+95, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-30, centerPenalizeA4.Ny+95, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-40, centerPenalizeA4.Ny+95, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-50, centerPenalizeA4.Ny+95, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-60, centerPenalizeA4.Ny+95, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-75, centerPenalizeA4.Ny+90, 13);

centerPenalizeA4.Nx = Nx/2+57;
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-60, centerPenalizeA4.Ny+95, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-60, centerPenalizeA4.Ny+90, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-60, centerPenalizeA4.Ny+80, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-60, centerPenalizeA4.Ny+70, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-60, centerPenalizeA4.Ny+60, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-60, centerPenalizeA4.Ny+50, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-60, centerPenalizeA4.Ny+40, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-60, centerPenalizeA4.Ny+30, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-60, centerPenalizeA4.Ny+20, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-60, centerPenalizeA4.Ny+10, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-60, centerPenalizeA4.Ny+0, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-60, centerPenalizeA4.Ny-10, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-60, centerPenalizeA4.Ny-20, 13);

centerPenalizeA4.Nx = Nx/2+57;
centerPenalizeA4.Ny = Ny/2-35;
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-50, centerPenalizeA4.Ny-10, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-40, centerPenalizeA4.Ny-10, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-30, centerPenalizeA4.Ny-10, 13);
penalizeMapA4 =  penalizeMapA4 + makeDisc(Nx, Ny, centerPenalizeA4.Nx-20, centerPenalizeA4.Ny-5, 13);

penalizeMap =  5 * (penalizeMapA1 >= 1) + (penalizeMapA2 >= 1) + (penalizeMapA3 >= 1) + 2 * (penalizeMapA4 >= 1);
%%
subplot(1, 2, 2); 
imagesc(double(segments) + (targetMap*-350) + (penalizeMap*200));
axis image;
colorbar;
title('Blob map');  