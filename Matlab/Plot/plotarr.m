function plotarr( ARRFIL, irr, ird, isd )

% plot the arrivals calculated by BELLHOP
%
% usage:
% plotarr( filename, irr, ird, isd )
% where:
% irr = index of receiver range
% ird = index of receiver depth
% isd = index of source   depth
%
% mbp, April 2009

% read

Narrmx = 5000;
[ Arr, Pos ] = read_arrivals_asc( ARRFIL, Narrmx );
%[ Arr, Pos ] = read_arrivals_bin( ARRFIL, Narrmx );

% stem plot for a single receiver
figure
Narr = Arr.Narr( irr, ird, isd );
stem( Arr.delay( irr, 1:Narr, ird, isd ), abs( Arr.A( irr, 1:Narr, ird, isd ) ) )
xlabel( 'Time (s)' )
ylabel( 'Amplitude' )
title( [ 'Sd = ', num2str( Pos.s.depth( isd ) ), ...
   ' m    Rd = ', num2str( Pos.r.depth( ird ) ), ...
   ' m    Rr = ', num2str( Pos.r.range( irr ) ), ' m' ] )

% depth-time stem plot
figure
for ird1 = 1 : size( Arr.A, 3 )
   Narr = Arr.Narr( irr, ird1, isd );
   stem3( Arr.delay( irr, 1:Narr, ird1, isd ), Pos.r.depth( ird1 ) * ones( length( Arr.delay( irr, 1:Narr, ird1, isd ) ), 1 ), ...
       abs( Arr.A( irr, 1:Narr, ird1, isd ) ) )
hold on
end

xlabel( 'Time (s)' )
ylabel( 'Depth (m)' )
title( [ 'Sd = ', num2str( Pos.s.depth( isd ) ), ' m    Rr = ', num2str( Pos.r.range( irr ) ), ' m' ] )

% range-time stem plot
figure
for irr = 1 : size( Arr.A, 1 )
   Narr = Arr.Narr( irr, ird, isd );
   stem3( Arr.delay( irr, 1:Narr, ird, isd ), Pos.r.range( irr ) * ones( length( Arr.delay( irr, 1:Narr, ird, isd ) ), 1 ), ...
       abs( Arr.A( irr, 1:Narr, ird, isd ) ) )
hold on
end

xlabel( 'Time (s)' )
ylabel( 'Range (m)' )
title( [ 'Sd = ', num2str( Pos.s.depth( isd ) ), ' m    Rd = ', num2str( Pos.r.depth( ird ) ), ' m' ] )
