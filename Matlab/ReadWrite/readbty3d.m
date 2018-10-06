function [ xBot, yBot, zBot, NbtyPtsx, NbtyPtsy ] = readbty3d( btyfil )

% readbty3d

if ( strcmp( btyfil, 'BTYFIL' ) == 0 && isempty( strfind( btyfil, '.bty' ) ) )
    btyfil = [ btyfil '.bty' ]; % append extension
end

fid = fopen( btyfil, 'r' );
if ( fid == -1 )
    error( 'Bathymetry file does not exist' )
end

btyType = fgetl( fid );

% Extract option letter between the quotes
nchars = strfind( btyType, '''' );   % find quotes
btyType = [ btyType( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) blanks( 2 - ( nchars( 2 ) - nchars( 1 ) ) ) ];

switch ( btyType )
    case ( 'R' )
        disp( 'Piecewise-linear approximation to bathymetry' )
    case ( 'C' )
        disp( 'Curvilinear approximation to bathymetry' )
    otherwise
        fclose all;
        disp( btyType )
        error( 'Fatal error: Unknown option for bathymetry type' )
end

%% x values
[ xBot, NbtyPtsx ] = readvector( fid );

% NbtyPtsx = fscanf( fid, '%i', 1 );
fprintf( 'Number of bathymetry points in x = %i \n\n', NbtyPtsx )
fprintf( ' x (km) \n' )
% 
% xBot = zeros( 1, NbtyPtsx );
% 
for ii = 1 : NbtyPtsx
%     xBot( ii ) = fscanf( fid, '%f', 1 );
   if ( ii < 50 || ii == NbtyPtsx )   % echo up to 51 values
      fprintf( '%9.5g \n', xBot( ii ) );
   end     
end

%xlimits = fscanf( fid, '%f ', 2 )
%x = linspace( xlimits( 1 ), xlimits( 2 ), NbtyPtsx );
%fgetl( fid );

%% y values
[ yBot, NbtyPtsy ] = readvector( fid );

% NbtyPtsy = fscanf( fid, '%i', 1 );
fprintf( 'Number of bathymetry points in y = %i \n\n', NbtyPtsy )
fprintf( ' y (km) \n' )
% 
% yBot = zeros( 1, NbtyPtsy );
% 
for ii = 1 : NbtyPtsy
%     yBot( ii ) = fscanf( fid, '%f', 1 );
   if ( ii < 50 || ii == NbtyPtsy )   % echo up to 51 values
         fprintf( '%9.5g \n', yBot( ii ) );
   end
end

%ylimits = fscanf( fid, '%f ', 2 );
%y = linspace( ylimits( 1 ), ylimits( 2 ), NbtyPtsy );
%fgetl( fid );

zBot = fscanf( fid, '%f ', [ NbtyPtsx, NbtyPtsy ] );
zBot = zBot';

fclose( fid );  % close the bathymetry file
