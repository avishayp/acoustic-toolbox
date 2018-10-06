function [ TitleEnv, Opt, Comp, MLimit, NProf, rProf, R, Pos ] = read_flp( fileroot )

% Read the field.flp file
% Usage:
%    read_flp
% mbp 4/09

fid = fopen( [ fileroot '.flp' ], 'r' );    % open the field parameters file

if ( fid == -1 )
   fileroot
   error( 'flp file does not exist' )
end

% read Title
TitleEnv = fgetl( fid );
% Extract letters between the quotes
nchars = strfind( TitleEnv, '''' );   % find quotes
if ( ~isempty( nchars ) )
   TitleEnv   = TitleEnv( nchars( 1 ) + 1 : nchars( 2 ) - 1 );
end
disp( TitleEnv )

% read Opt
Opt = fgetl( fid );
% Extract letters between the quotes
nchars = strfind( Opt, '''' );   % find quotes
Opt    = Opt( nchars( 1 ) + 1 : nchars( 2 ) - 1 );
disp( Opt )

if ( length( Opt ) <= 2 )
   Opt( 3 : 3 ) = 'O';   % default beampattern is omni
end

% select the component
if ( length( Opt ) >= 3 )
   Comp = Opt( 3 : 3 );
else
   Comp = 'P';
end

% read MLimit
MLimit   = fscanf( fid, '%i', 1 );
fprintf( 'MLimit = %i \n\n', MLimit )
fgetl( fid );

% read profile info
NProf = fscanf( fid, '%i', 1 );
fprintf( '\nNumber of profiles   = %i \n', NProf )
fprintf( '\nProfile ranges (km)\n' );

fgetl( fid );

rProf = fscanf( fid, '%f', NProf );
fgetl( fid );

if ( NProf < 10 )
   fprintf( '%8.2f  \n', rProf )   % print all the depths
else
   fprintf( '%8.2f ... %8.2f \n', rProf( 1 ), rProf( end ) ) % print first, last depth
end

if NProf > 2
   rProf = linspace( rProf( 1 ), rProf( 2 ), NProf ); % generate vector of profile ranges
   warning( 'Producing profile ranges by interpolation between rProf(1) and rProf(2)' )
end

% read receiver ranges
R = readr( fid );
R = 1000.0 * R; % convert km to m

Pos = readsdrd( fid );
