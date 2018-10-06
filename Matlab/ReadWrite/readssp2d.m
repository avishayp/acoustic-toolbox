function [ cmat, rProf, NProf, NSSP ] = readssp2d( sspfil )

% Read the SSPFIL

% open the file

fid = fopen( sspfil, 'r' );

if ( fid == -1 )
    warndlg( 'No SSPFIL file with that name exists', 'Warning' );
end

% read
NProf = fscanf( fid, '%i', 1 );

rProf = fscanf( fid, '%f', NProf );
cmat  = fscanf( fid, '%f', [ NProf, inf ] );
cmat  = cmat';

NSSP = size( cmat, 1 );

fclose( fid );