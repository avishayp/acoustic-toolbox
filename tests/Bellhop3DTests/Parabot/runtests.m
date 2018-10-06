% run the ParaBot test case
% p517 computational acoustics

global units
units = 'km';
%%
makebty              % make the bathymetry
figure
plotbty3d ParaBot
shading flat

% ray trace
copyfile( 'ParaBot.bty', 'ParaBot_ray.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot_ray.ati' )   % copy over the bathymetry file

bellhop3d ParaBot_ray     % run BELLHOP3D on the wedge3d test case

hold on
plotray3d ParaBot_ray.ray

delete( 'ParaBot_ray.bty' )
delete( 'ParaBot_ray.ati' )

%%
copyfile( 'ParaBot.bty', 'ParaBot2D.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot2D.ati' )   % copy over the bathymetry file

bellhop3d ParaBot2D     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshd( 'ParaBot2D.shd' )
caxisrev( [ 60 100 ] )

delete( 'ParaBot2D.bty' )
delete( 'ParaBot2D.ati' )

%%
% 3d run (GeoHat Cartesian)

copyfile( 'ParaBot.bty', 'ParaBot3DHatCart.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot3DHatCart.ati' )   % copy over the bathymetry file

bellhop3d ParaBot3DHatCart     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshd( 'ParaBot3DHatCart.shd' )
caxisrev( [ 60 100 ] )

delete( 'ParaBot3DHatCart.bty' )
delete( 'ParaBot3DHatCart.ati' )

%%
% 3d run (GeoHat Ray centered)

copyfile( 'ParaBot.bty', 'ParaBot3DHatRaycen.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot3DHatRaycen.ati' )   % copy over the bathymetry file

bellhop3d ParaBot3DHatRaycen     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshd( 'ParaBot3DHatRaycen.shd' )
caxisrev( [ 60 100 ] )

delete( 'ParaBot3DHatRaycen.bty' )
delete( 'ParaBot3DHatRaycen.ati' )


%%
% 3d run (GeoGaussian)

copyfile( 'ParaBot.bty', 'ParaBot3DGaussian.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot3DGaussian.ati' )   % copy over the bathymetry file

bellhop3d ParaBot3DGaussian     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshd( 'ParaBot3DGaussian.shd' )
caxisrev( [ 60 100 ] )

delete( 'ParaBot3DGaussian.bty' )
delete( 'ParaBot3DGaussian.ati' )

delete( 'ParaBot.bty' )
delete( 'ParaBot.ati' )
%print -depsc2 ParaBot3DGaussian
%print -djpeg ParaBot3DGaussian

