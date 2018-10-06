% Run tests to verify the point vs. line source options

global units
units = 'm';

% BELLHOP point source cases
bellhopM 'freePointB'
plotshd( 'freePointB.shd.mat', 3, 1, 1 );
axis equal
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

bellhopM 'freePoint_gbtB'
plotshd( 'freePoint_gbtB.shd.mat', 3, 1, 2 );
axis equal
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

% SCOOTER

scooterM 'freeSPoint'
plotshd( 'freeSPoint.shd.mat', 3, 1, 3 );
axis equal
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

%%
% incoherent TL runs

% BELLHOP point source cases
bellhopM 'freePointB_Inc'
plotshd( 'freePointB_Inc.shd.mat', 3, 1, 1 );
axis equal
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

bellhopM 'freePoint_gbtB_Inc'
plotshd( 'freePoint_gbtB_Inc.shd.mat', 3, 1, 2 );
axis equal
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

plotshd( 'freeSPoint.shd.mat', 3, 1, 3 );
axis equal
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )
