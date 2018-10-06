% run the seamount test case
% p517 computational acoustics

global units
units = 'km';
%%
makebty              % make the bathymetry
figure
plotbty3d Seamount
shading flat

%% ray trace
copyfile( 'Seamount.bty', 'Seamount_ray.bty' )   % copy over the bathymetry file
bellhop3d Seamount_ray

hold on
plotray3d Seamount_ray.ray

%%
copyfile( 'Seamount.bty', 'Seamount2D.bty' )   % copy over the bathymetry file

bellhop3d Seamount2D

% polar plot of the TL
figure
plotshdpol( 'Seamount2D.shd', 3, 0, 100 )
caxisrev( [ 40 100 ] )
axis( [ -6 3 0 9 ] )

%%
% 3d run (GeoHat Cartesian)

copyfile( 'Seamount.bty', 'Seamount3DHatcart.bty' )   % copy over the bathymetry file

bellhop3d Seamount3DHatcart

% polar plot of the TL
figure
plotshdpol( 'Seamount3DHatcart.shd', 3, 0, 100 )
caxisrev( [ 40 100 ] )
axis( [ -6 3 0 9 ] )

%%
% 3d run (GeoHat Ray centered)

copyfile( 'Seamount.bty', 'Seamount3DHatRaycen.bty' )   % copy over the bathymetry file

bellhop3d Seamount3DHatRaycen

% polar plot of the TL
figure
plotshdpol( 'Seamount3DHatRaycen.shd', 3, 0, 100 )
caxisrev( [ 40 100 ] )
axis( [ -6 3 0 9 ] )

%%
% 3d run (GeoGaussian)

copyfile( 'Seamount.bty', 'Seamount3DGaussian.bty' )   % copy over the bathymetry file

bellhop3d Seamount3DGaussian

% polar plot of the TL
figure
plotshdpol( 'Seamount3DGaussian.shd', 3, 0, 100 )
caxisrev( [ 40 100 ] )
axis( [ -6 3 0 9 ] )

%print -depsc2 Seamount3DGaussian
%print -djpeg Seamount3DGaussian

