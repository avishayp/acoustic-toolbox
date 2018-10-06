% generate parabolic bottom bathymetry
% mbp Dec. 27, 2005
% example from:
% R. W. McGirr, D.B. King, J.A. Davis, J. Campbell, "An evaluation of
% range-dependent ray theory models", NORDA report 115.

b = 250000;
c = 250;

z = 0 : 5 : 5000;    % for TL comparison to virtual source method
x = c * ( ( z / 0.002 / b ).^2 - 1 ) / 1000;
y = x;
nx = length( x );
ny = length( y );

btyfil = 'ParaBot.bty';
atifil = 'ParaBot.ati';
%interp_type = 'R';
interp_type = 'C';

% xmin = -20;
% xmax = 20;
% ymin = -20;
% ymax = 20;
% 
% nx = 210;
% ny = 210;
% x = linspace( xmin, xmax, nx );
% y = linspace( ymin, ymax, ny );

[ X, Y ] = meshgrid( x, y );
R = 1000 * sqrt( X.^2 + Y.^2 );   % range in meters
z = 0.002 * b * sqrt( 1 + R / c );

% altimetry
fid = fopen( atifil, 'w' );
fprintf( fid, '''%c'' \n', interp_type );

fprintf( fid, '%i \r\n', nx );
fprintf( fid, '%f ', x );
%fprintf( fid, '%f %f /', xmin, xmax );
fprintf( fid, '\r\n');

fprintf( fid, '%i \r\n', ny );
fprintf( fid, '%f ', y );
%fprintf( fid, '%f %f /', ymin, ymax );
fprintf( fid, '\r\n');

for iy = 1 : ny
   fprintf( fid, '%f ', -z( iy, : ) );
   fprintf( fid, '\r\n');
end

fclose( fid );

% bathymetry

fid = fopen( btyfil, 'w' );
fprintf( fid, '''%c'' \n', interp_type );

fprintf( fid, '%i \r\n', nx );
fprintf( fid, '%f ', x );
%fprintf( fid, '%f %f /', xmin, xmax );
fprintf( fid, '\r\n');

fprintf( fid, '%i \r\n', ny );
fprintf( fid, '%f ', y );
%fprintf( fid, '%f %f /', ymin, ymax );
fprintf( fid, '\r\n');

for iy = 1 : ny
   fprintf( fid, '%f ', z( iy, : ) );
   fprintf( fid, '\r\n');
end

fclose( fid );
