function field( FileRoot, modes )

% calculates the field using modes produced by KRAKEN
% parameters read from field.flp
%
% usage: field( filename, modes )
% mbp

% Comp     = 'P';   % select component (P, H, V) (Pressure, Horizontal, Vertical)

clear read_modes_bin % to force rewind to beginning of mode file

[ TitleEnv, Opt, Comp, MLimit, NProf, rProf, R, Pos ] = read_flp( FileRoot );

rd   = Pos.r.depth;   % receiver depths
Nrd  = length( rd );
Nr   = length( R );   % receiver ranges
freq = 0;            % assume single frequency run
% optionally read in a source beam pattern
SBP = Opt( 3 : 3 );
SrcBmPat = readpat( FileRoot, SBP )

if ( NProf == 1 )               % Range-independent case
   
   filename = [ FileRoot '.mod' ];
   
   if nargin == 1
      Modes = read_modes( filename, freq );
   else
      Modes = read_modes( filename, freq, modes );
   end

   freqVec = Modes.freqVec;
   MSrc = length( Modes.k );
   M    = min( MLimit, MSrc );        % Set number of propagating modes
   
   % weight modes by mode excitation
   zs  = Pos.s.depth;
   isd = find( Modes.z >= zs );    % index of source depth
   isd = isd( 1 );
   C   = Modes.phi( isd, : ).';    % column vector with modal weights
   
   % apply the source beam pattern
   if ( SBP == '*' )
      c = 1500;   % reference sound speed, should be speed at the source depth
      omega = 2 * pi * Modes.freqVec( 1 );
      kz2 = omega^2 / c^2 - Modes.k.^2;   % vertical wavenumber squared
      kz2( kz2 < 0 ) = 0;                 % remove negative values
      
      theta = atand( real( sqrt( kz2 ) ./ Modes.k ) );   % calculate the angle in degrees
      S = interp1( SrcBmPat( :, 1 ), SrcBmPat( : , 2 ), theta );
      C = C .* S;         % apply the shading
   end
    
   % calculate mode values at receiver depths
   irdvec = zeros( 1, length( Pos.r.depth ) );
   for ii = 1 : length( Pos.r.depth )
      zr  = Pos.r.depth( ii );
      ird = find( Modes.z >= zr );	% index of source depth
      if ( numel( ird ) == 0 )
          ird( 1 ) = length( Pos.r.depth );
      end
      irdvec( ii ) = ird( 1 );
   end
   phiR = Modes.phi( irdvec, : );
   
   pressure = evalri( C, phiR, R, Modes.k, Opt, Comp );
else
   if ( Opt(2:2) == 'C' )       % Range-dependent case
      % Coupled mode theory
      clear evalcm % to force re-initialization of profile ranges
      pressure = evalcm( FileRoot, Pos, rProf, NProf, rd, Nrd, R, Nr, MLimit, Opt );
   else
      % Adiabatic mode theory
      pressure = evalad( FileRoot, Pos, rProf, NProf, rd, Nrd, R, Nr, M, Opt );
   end
end

PlotTitle   = TitleEnv;
PlotType    = 'rectilin  ';
freq        = 0.0; % Modes.freq;
atten       = 0.0;
Pos.r.range = R;

save( [ FileRoot '.shd.mat' ], 'PlotTitle', 'PlotType', 'freqVec', 'atten', 'Pos', 'pressure' )
