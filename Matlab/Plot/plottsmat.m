function plotts( filename )
% plot a time series from the given file

% open the file
load( [ filename '.rts.mat' ] )
rd  = Pos.r.depth;
nrd = length( rd );
t   = tout;
RTS = RTS';

% plot
figure
orient tall
title( PlotTitle )
set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature

% scale all time series so that max is unity
scale  = max( max( RTS( :, : ) ) );
RTS    = RTS / scale * rd( nrd )/nrd;
offset = linspace( rd( 1 ), rd( nrd ), nrd );

hold on
threshold = -1e30;

for ird = 1:nrd
   ii = find( RTS( :, ird ) >  threshold );
   jj = find( RTS( :, ird ) <= threshold );
   h=area( t( ii ), RTS( ii, ird ) + offset( ird ) ); % pos. part shading under line
   %ylabel( [ 'Rd = ', num2str( rd( ird ) ) ] );
   set( h, 'BaseValue', offset( ird ) );
   plot( t( jj ), RTS( jj, ird ) + offset( ird ) ); % negative part just line
end

xlabel( 'Time (s)' );
%ylabel( [ 'Rd = ', num2str( rd( nrd ) ) ] );

%set(1,'PaperPosition', [ 0.25 0.00 5.5 7.0 ] )
%print -deps bellhop.ps
