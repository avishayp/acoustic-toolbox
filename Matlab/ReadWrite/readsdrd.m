function Pos = readsdrd( fid )

% Read source depths and receiver depths
%
% Variable 'Pos' is a structure:
% Pos.r.depth = vector of receiver depths
% Pos.Nrd     = number of receiver depths
% Pos.s.depth = vector of source depths
% Pos.Nsd     = number of source depths

%%
% source depths
fprintf( '\n_______________________ \n' )

[ Pos.s.depth, Pos.Nsd ] = readvector( fid );

%Pos.Nsd = fscanf( fid, '%i', 1 );
fprintf( '\nNumber of source depths   = %i \n', Pos.Nsd )
fprintf( '\nSource depths (m)\n' );

%fgetl( fid );
%Pos.s.depth = fscanf( fid, '%f', Pos.Nsd );

if ( Pos.Nsd < 10 )
   fprintf( '%8.2f  \n', Pos.s.depth )   % print all the depths
else
   fprintf( '%8.2f ... %8.2f \n', Pos.s.depth( 1 ), Pos.s.depth( end ) ) % print first, last depth
end
%fgetl( fid );

if Pos.Nsd > 2
  %Pos.s.depth = linspace( Pos.s.depth( 1 ), Pos.s.depth( 2 ), Pos.Nsd ); % generate vector of receiver depths
  disp( 'Producing source depths by interpolation between sd(1) and sd(2)' )
end

%%
% receiver depths
fprintf( '\n_______________________ \n' )

[ Pos.r.depth, Pos.Nrd ] = readvector( fid );

%Pos.Nrd = fscanf( fid, '%i', 1 );
fprintf( '\nNumber of receivers depths = %i \n', Pos.Nrd )
fprintf( '\nReceiver depths (m)\n' );

%fgetl( fid );
%Pos.r.depth = fscanf( fid, '%f', Pos.Nrd );

if ( Pos.Nrd < 10 )
   fprintf( '%8.2f  \n', Pos.r.depth )   % print all the depths
else
   fprintf( '%8.2f ... %8.2f \n', Pos.r.depth( 1 ), Pos.r.depth( end ) ) % print first, last depth
end

disp( '  ' )

if Pos.Nrd > 2
  % Pos.r.depth = linspace( Pos.r.depth( 1 ), Pos.r.depth( 2 ), Pos.Nrd ); % generate vector of receiver depths
  disp( 'Producing receiver depths by interpolation between rd(1) and rd(2)' )
end

%fgetl( fid );
fprintf( '\n' )