% Penetrable Wedge test case

global units
units = 'km';

%% GeoHat beams
bellhop3d pwedge3dRayCen     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol( 'pwedge3dRayCen.shd', 0, -20, 100 )
caxisrev( [ 80 100 ] )
axis image
axis( [ 0 70 -50 0 ] )

%% GeoHat beams, Cart. coord.
bellhop3d pwedge3dCart     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol( 'pwedge3dCart.shd', 0, -20, 100 )
caxisrev( [ 80 100 ] )
axis image
axis( [ 0 70 -50 0 ] )

%% GeoHat beams 2D
bellhop3d pwedge2d     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol( 'pwedge2d.shd', 0, -20, 100 )
caxisrev( [ 80 100 ] )
axis image
axis( [ 0 70 -50 0 ] )
%% GeoHat beams 2D rotated
bellhop3d pwedge2d_rot     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol( 'pwedge2d_rot.shd', -20, 0, 100 )
caxisrev( [ 80 100 ] )
axis image
axis( [ -50 0 -70 0 ] )
%% GeoHat beams 3D rotated
bellhop3d pwedge3d_rot     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol( 'pwedge3d_rot.shd', -20, 0, 100 )
caxisrev( [ 80 100 ] )
axis image
axis( [ -50 0 -70 0 ] )

%% Side view
bellhop3d slice2d     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshd( 'slice2d.shd' )
caxisrev( [ 60 100 ] )
% axis( [ -50 0 -70 0 ] )
