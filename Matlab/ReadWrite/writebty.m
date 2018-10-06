function writebty( btyfil, interp_type, rngdep )
% Write a bathymetry from the workspace variables

Npts = length( rngdep( :, 1 ) );

fid = fopen( btyfil, 'w' );
fprintf( fid, '''%c'' \n', interp_type );
fprintf( fid, '%i \r\n', Npts );
fprintf( fid, '%f %f \r\n', rngdep' );

fclose( fid );

