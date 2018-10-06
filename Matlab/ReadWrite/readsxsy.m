function [ sx, sy, Nsx, Nsy ] = readsxsy( fid )

% Read source x and y coordinates
fprintf( '\n_______________________ \n' )

% x coordinates

[ sx, Nsx ] = readvector( fid );

%Nsx = fscanf( fid, '%i', 1 );
fprintf( '\nNumber of source x coordinates = %i \n', Nsx )
%fgetl( fid );

%sx = fscanf( fid, '%f', Nsx );

fprintf( '\nSource x coordinates (km) \n' )
fprintf( '%f ', sx )
fprintf( '\n' )

%if Nsx > 2
%   sx = linspace( sx( 1 ), sx( 2 ), Nsx )'; % generate vector of receiver ranges
%end

%fgetl( fid );

% y coordinates

[ sy, Nsy ] = readvector( fid );

%Nsy = fscanf( fid, '%i', 1 );
fprintf( '\nNumber of source y coordinates = %i \n', Nsy )
%fgetl( fid );

%sy = fscanf( fid, '%f', Nsy );

fprintf( '\nSource y coordinates (km) \n' )
fprintf( '%f ', sy )
fprintf( '\n' )

%if Nsy > 2
%   sy = linspace( sy( 1 ), sy( 2 ), Nsy )'; % generate vector of receiver ranges
%end

%fgetl( fid );

sx = 1000.0 * sx;   % convert km to m
sy = 1000.0 * sy;

