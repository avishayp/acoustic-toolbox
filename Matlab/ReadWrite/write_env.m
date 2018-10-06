function write_env( envfil, model, TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, varargin )

% Write an environmental file
% mbp 2009
% note: I had an unusual case with a parabolic mirror where round-off in
% the receiver ranges was a problem (+/- 1 m in range)
% Changed the output format from %6.2f to %6f to accommodate that.
% Similar issues may occur elsewhere in the code below ...

if ( strcmp( envfil, 'ENVFIL' ) == 0 && ~strcmp( envfil( end-3: end ), '.env' ) )
  envfil = [ envfil '.env' ]; % append extension
end

if ( size( varargin ) == 0 )
    fid = fopen( envfil, 'w' );   % create new envfil
else
    fid = fopen( envfil, 'a' );   % append to existing envfil
end

if ( fid == -1 )
    envfil
    error( 'Unable to create environmental file', 'write_env' );
end

model = upper( model );   % convert to uppercase

fprintf( fid, '''%s'' ! Title \r\n', TitleEnv );
fprintf( fid, '%8.2f  \t \t \t ! Frequency (Hz) \r\n', freq );
fprintf( fid, '%5i    \t \t \t ! NMedia \r\n', SSP.NMedia );
fprintf( fid, '''%s'' \t \t \t ! Top Option \r\n', Bdry.Top.Opt );

if ( Bdry.Top.Opt( 2:2 ) == 'A' )
    fprintf( fid, '    %6.2f %6.2f %6.2f %6.2g %6.2f %6.2f /  \t ! upper halfspace \r\n', SSP.depth( 1 ), ...
        Bdry.Top.HS.alphaR, Bdry.Top.HS.betaR, Bdry.Top.HS.rho, Bdry.Top.HS.alphaI, Bdry.Top.HS.betaI );
end

% SSP
for medium = 1 : SSP.NMedia
    
    fprintf( fid, '%5i %4.2f %6.2f \t ! N sigma depth \r\n', SSP.N( medium ), SSP.sigma( medium ), SSP.depth( medium+1 ) );
    for ii = 1 : length( SSP.raw( medium ).z )
        fprintf( fid, '\t %6.2f %6.2f %6.2f %6.2g %10.6f %6.2f / \t ! z c cs rho \r\n', ...
            [ SSP.raw( medium ).z( ii ) ...
              SSP.raw( medium ).alphaR( ii ) SSP.raw( medium ).betaR( ii ) SSP.raw( medium ).rho( ii ) ...
              SSP.raw( medium ).alphaI( ii ) SSP.raw( medium ).betaI( ii ) ] );
    end
end

% lower halfspace
fprintf( fid, '''%s'' %6.2f  \t \t ! Bottom Option, sigma\r\n', Bdry.Bot.Opt, 0.0 ); % SSP.sigma( 2 ) );

if ( Bdry.Bot.Opt( 1:1 ) == 'A' )
    fprintf( fid, '    %6.2f %6.2f %6.2f %6.2g %6.2f %6.2f /  \t ! lower halfspace \r\n', SSP.depth( SSP.NMedia+1 ), ...
        Bdry.Bot.HS.alphaR, Bdry.Bot.HS.betaR, Bdry.Bot.HS.rho, Bdry.Bot.HS.alphaI, Bdry.Bot.HS.betaI );
end

if( strmatch( model, strvcat( 'SCOOTER', 'KRAKEN', 'KRAKENC', 'SPARC' ), 'exact' ) )
    fprintf( fid, '%6.0f %6.0f \t \t ! cLow cHigh (m/s) \r\n', cInt.Low, cInt.High );   % phase speed limits
    fprintf( fid, '%8.2f \t \t \t ! RMax (km) \r\n', RMax );    % maximum range
end

% source depths

fprintf( fid, '%5i \t \t \t \t ! NSD', length( Pos.s.depth ) );

if ( length( Pos.s.depth ) >= 2 && equally_spaced( Pos.s.depth ) )
    fprintf( fid, '\r\n    %6f %6f ', Pos.s.depth( 1 ), Pos.s.depth( end ) );
else
    fprintf( fid, '\r\n    %6f  ', Pos.s.depth );
end

fprintf( fid, '/ \t ! SD(1)  ... (m) \r\n' );

% receiver depths

fprintf( fid, '%5i \t \t \t \t ! NRD', length( Pos.r.depth ) );

if ( length( Pos.r.depth ) >= 2 && equally_spaced( Pos.r.depth ) )
    fprintf( fid, '\r\n    %6f %6f ', Pos.r.depth( 1 ), Pos.r.depth( end ) );
else
    fprintf( fid, '\r\n    %6f  ', Pos.r.depth );
end

fprintf( fid, '/ \t ! RD(1)  ... (m) \r\n' );

% receiver ranges
if ( strcmp( model, 'BELLHOP' ) ||  strcmp( model, 'FirePE' ) )
    % receiver ranges
    fprintf( fid, '%5i \t \t \t \t ! NRR', length( Pos.r.range ) );
    
    if ( length( Pos.r.range ) >= 2 && equally_spaced( Pos.r.range ) )
        fprintf( fid, '\r\n    %6f %6f ', Pos.r.range( 1 ), Pos.r.range( end ) );
    else
        fprintf( fid, '\r\n    %6f  ', Pos.r.range );
    end
    fprintf( fid, '/ \t ! RR(1)  ... (km) \r\n' );
    write_bell( fid, Beam );
end

fclose( fid );
