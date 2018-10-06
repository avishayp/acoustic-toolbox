% Parabolic Bottom test
% mbp

% linear boundary interpolation

make_bdry( 'L' )

% the rays:
bellhopM ParaBot
%figure; plotray RAYFIL
%plotbty 'ParaBot'   % superimpose a bathymetry plot

% TL: Geometric ray theory
copyfile( 'ParaBot.bty', 'ParaBotTLGeom.bty' );
bellhopM ParaBotTLGeom
plotshd( 'ParaBotTLGeom.mat', 2, 1, 1 )
caxisrev( [ 60 100 ] )

plotbty 'ParaBot'   % superimpose a bathymetry plot

% TL: Gaussian beams
copyfile( 'ParaBot.bty', 'ParaBotTLGB.bty' );
bellhopM ParaBotTLGB
plotshd( 'ParaBotTLGB.mat', 2, 1, 2 )
caxisrev( [ 60 100 ] )

plotbty 'ParaBot'   % superimpose a bathymetry plot

% curvilinear boundary interpolation

make_bdry( 'C' )

% the rays:
bellhopM ParaBot
%figure; plotray RAYFIL
%plotbty 'ParaBot'   % superimpose a bathymetry plot

% TL: Geometric ray theory
copyfile( 'ParaBot.bty', 'ParaBotTLGeom.bty' );
bellhopM ParaBotTLGeom
plotshd( 'ParaBotTLGeom.mat', 2, 1, 1 )
caxisrev( [ 60 100 ] )

plotbty 'ParaBot'   % superimpose a bathymetry plot

% TL: Gaussian beams
copyfile( 'ParaBot.bty', 'ParaBotTLGB.bty' );
bellhopM ParaBotTLGB
plotshd( 'ParaBotTLGB.mat', 2, 1, 2 )
caxisrev( [ 60 100 ] )

plotbty 'ParaBot'   % superimpose a bathymetry plot

