% run the Taiwan 3D test cases
% mbp, July 2011

% the polar and sideview plots should be done in separate BELLHOP runs for
% speed ...

global units
units = 'km';

%Bathy = LoadBathymetry( 'Taiwan.xyz' );   % load the xyz file from ETOPO1
%writebty3d( 'Taiwan.bty', Bathy );        % save it in my bty format
figure
plotbty3d Taiwan.bty                          % howzit look?

% makebty;             % make the bathymetry
%%

% 2D Transmission Loss runs

fileroot = 'TaiwanNx2D'
   
copyfile( 'Taiwan.bty', [ fileroot '.bty' ] )
copyfile( 'Taiwan.ssp', [ fileroot '.ssp' ] )

bellhop3d( fileroot )

% polar plot of the TL
figure
plotshdpol( [ fileroot '.shd' ], [ 175 300 ], [ 200 300 400 ], 500 )
caxisrev( [ 60 100 ] )

% sideview of the TL
figure
plotshd(    [ fileroot '.shd' ] )
caxisrev( [ 60 100 ] )

delete( [ fileroot '.bty' ] )
delete( [ fileroot '.ssp' ] )

%%

% 3D Transmission Loss runs

fileroot = 'Taiwan3D'
   
copyfile( 'Taiwan.bty', [ fileroot '.bty' ] )
copyfile( 'Taiwan.ssp', [ fileroot '.ssp' ] )

bellhop3d( fileroot )

% polar plot of the TL
figure
plotshdpol( [ fileroot '.shd' ], [ 175 300 ], [ 200 300 400 ], 500 )
caxisrev( [ 60 100 ] )

% sideview of the TL
figure
plotshd(    [ fileroot '.shd' ] )
caxisrev( [ 60 100 ] )

delete( [ fileroot '.bty' ] )
delete( [ fileroot '.ssp' ] )


%%

% ray trace runs

fileroot = 'TaiwanNx2D_ray'

copyfile( 'Taiwan.bty', [ fileroot '.bty' ] )
copyfile( 'Taiwan.ssp', [ fileroot '.ssp' ] )

bellhop3d( fileroot )

figure
plotbty3d( [ fileroot '.bty' ] )
hold on

plotray3d( [ fileroot '.ray' ] )

delete( [ fileroot '.bty' ] )
delete( [ fileroot '.ssp' ] )

%%


fileroot = 'Taiwan3D_ray'

copyfile( 'Taiwan.bty', [ fileroot '.bty' ] )
copyfile( 'Taiwan.ssp', [ fileroot '.ssp' ] )

bellhop3d( fileroot )

figure
plotbty3d( [ fileroot '.bty' ] )
hold on

plotray3d( [ fileroot '.ray' ] )

delete( [ fileroot '.bty' ] )
delete( [ fileroot '.ssp' ] )
