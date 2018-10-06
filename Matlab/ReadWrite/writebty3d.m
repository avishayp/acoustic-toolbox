function writebty3d( btyfil, Bathy )

% write a 3d bathymetry file for use by BELLHOP3D
% useage:
% writebty3d( btyfil, Bathy)
%
% mbp July 2011

% Some of my code uses Bathy.Lonkm ad .Latkm but I think this is confusing
% since they are Eastings and Northings at that point, not lat/longs.

%Bathy.X = Bathy.Lonkm;
%Bathy.Y = Bathy.Latkm;

Bathy.depth( isnan( Bathy.depth ) ) = 0.0;   % remove NaNs

interp_type = 'R';   % rectilinear grid
nx = length( Bathy.X );
ny = length( Bathy.Y );

fid = fopen( btyfil, 'w' );
fprintf( fid, '''%c'' \n', interp_type );

fprintf( fid, '%i \r\n', nx );
%fprintf( fid, '%f %f /', Bathy.X( 1 ), Bathy.X( end ) );

fprintf( fid, '%f ', Bathy.X( 1 : end ) );
fprintf( fid, '\r\n');

fprintf( fid, '%i \r\n', ny );
%fprintf( fid, '%f %f /', Bathy.Y( 1 ), Bathy.Y( end ) );
fprintf( fid, '%f ', Bathy.Y( 1 : end ) );
fprintf( fid, '\r\n');

for iy = 1 : ny
   fprintf( fid, '%9.3f ', Bathy.depth( iy, : ) );
   fprintf( fid, '\r\n');
end

fclose( fid );
