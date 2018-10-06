function [ theta ] = readtheta( fid )

% Read receiver angles
fprintf( '\n_______________________ \n' )
[ theta, Ntheta ] = readvector( fid );

%Ntheta = fscanf( fid, '%i', 1 );
fprintf( '\nNumber of receiver ranges = %i \n', Ntheta )
%fgetl( fid );

%theta = fscanf( fid, '%f', Ntheta );

fprintf( '\nReceiver angles (degrees) \n' )
fprintf( '%f ', theta )
fprintf( '\n' )

%if Ntheta > 2
%  theta = linspace( theta( 1 ), theta( 2 ), Ntheta )'; % generate vector of receiver angles
%end

%fgetl( fid );
