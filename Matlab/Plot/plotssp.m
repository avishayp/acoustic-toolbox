function plotssp( envfil )
% plotssp.m
% Plots the sound speed profile
% mike porter 3/2009

[ ~, ~, SSP, ~, fid ] = read_env( envfil, 'KRAKEN' ); % read in the environmental file

%figure
hold on

for medium = 1 : SSP.NMedia
   hh = plot( real( SSP.raw( medium ).alphaR ), SSP.raw( medium ).z );
   set( hh, 'LineWidth', 2 );
   if ( any( SSP.raw( medium ).betaR ) )
      hh = plot( real( SSP.raw( medium ).betaR ), SSP.raw( medium ).z, 'r:' );
      set( hh, 'LineWidth', 2 );
   end
end

set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
%axis IJ

xlabel( 'Sound Speed (m/s)' )
ylabel( 'Depth (m)' )

