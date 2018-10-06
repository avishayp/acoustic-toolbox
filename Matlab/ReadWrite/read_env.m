function  [ TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, fid ] = read_env( envfil, model )
% Read an environmental file

% should check here that model is one of the known ones

disp( '' )
disp( '' )
disp( '' )
disp( '__________________________________________________________________' );

if ( ~strcmp( envfil, 'ENVFIL' ) && ~strcmp( envfil( end-3: end ), '.env' ) )
    envfil = [ envfil '.env' ]; % append extension
end

model = upper( model );   % convert to uppercase

[ TitleEnv, freq, SSP, Bdry, fid ] = read_env_core( envfil );    % read in the environmental file

if( strmatch( model, strvcat( 'SCOOTER', 'KRAKEN', 'KRAKENC', 'KRAKEL', 'SPARC' ), 'exact' ) )
    cInt.Low   = fscanf( fid, '%f', 1 );   % lower phase speed limit
    cInt.High  = fscanf( fid, '%f', 1 );   % upper phase speed limit
    fprintf( '\n cLow = %8.1f m/s  cHigh = %8.1f m/s\n', cInt.Low, cInt.High )
    fgetl( fid );
    
    RMax  = fscanf( fid, '%f', 1 );   % read max range, Rmax
    fprintf( 'RMax = %f km \n', RMax )
    fgetl( fid );
    
else   % dummy values for BELLHOP
    cInt.Low  = 1500;
    cInt.High = 1e9;
    RMax      = 100;
end

% BELLHOP3D has x-y coordinates of sources as well
if ( strcmp( model, 'BELLHOP3D' ) )
    [ sx, sy, Nsx, Nsy ] = readsxsy( fid );
    Pos.s.x = sx;
    Pos.s.y = sy;
    Pos.Nsx = Nsx;
    Pos.Nsy = Nsy;
end

% !!! check: does this delete Pos.s.x, etc. when executed
Pos = readsdrd( fid );                           % read in the source and receiver depths

if ( strcmp( model, 'BELLHOP' ) )
    Pos.r.range = readr( fid );     % read in receiver ranges
    Pos.Nrr     = length( Pos.r.range );
    Beam        = read_bell( fid, Bdry, freq, Bdry.Bot.depth, Bdry.Top.depth, Pos.r.range( end ) );
elseif ( strcmp( model, 'BELLHOP3D' ) )
    Pos.r.range = readr( fid );     % read in receiver ranges
    Pos.Nrr     = length( Pos.r.range );
    Pos.theta   = readtheta( fid );
    Pos.Ntheta  = length( Pos.theta );
    Beam        = read_bell( fid, Bdry, freq, Bdry.Bot.depth, Bdry.Top.depth, Pos.r.range( end ) );
else   % dummy value for models that don't use Beam parameters
    Beam.RunType    = 'CG';
    Beam.Nbeams     = 0;
    Beam.alpha( 1 ) = -15;
    Beam.alpha( 2 ) = +15;
    Beam.Box.z      = 1.05 * max( Pos.r.depth );
    Beam.Box.r      = 1.05 * RMax;
    Beam.deltas     = 0;
    Pos.r.range     = linspace( 0, RMax, 501 );   % set up receiver range vector
end

