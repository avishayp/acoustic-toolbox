function InfluenceGeoGaussian( zs, alpha, RunType, Dalpha, deltas  )

% Computes the beam influence, i.e.
% the contribution of a single beam to the complex pressure
% This version uses a beam representation in Cartesian coordinates

global ray Nsteps omega
global Pos
global U MxNarr

IBWin = 4;

% some ugly code to have RunType 'a' treated like 'A':
RunTypeE = RunType( 1 : 1 );
if ( RunTypeE == 'a' )
    RunTypeE = 'A';
end

DS      = sqrt( 2.0 ) * sin( omega * zs * ray( 1 ).Tray( 2 ) );   % Lloyd mirror pattern
q0      = ray( 1 ).c( 1 ) / Dalpha;   % Reference for J = Q0 / Q
caustic = 1;   % caustic is multiplied by i after each caustic
freq    = omega / ( 2 * pi );
lambda  = ray( 1 ).c( 1 ) / freq;

if ( RunType( 4:4 ) == 'R' )   % point source
    Rat1 = sqrt( abs( cos( alpha ) ) ) / sqrt( 2 * pi );
else
    Rat1 = 1 / sqrt( 2 * pi );
end

% *** step along beams ***

for is = 1 : max( Nsteps ) - 1
    if ( is > 1 )
        if ( ( ( ray( is ).q(  1 ) <= 0.0 ) && ( ray( is - 1 ).q(  1 ) > 0.0 ) ) || ...
             ( ( ray( is ).q(  1 ) >= 0.0 ) && ( ray( is - 1 ).q(  1 ) < 0.0 ) ) )
            caustic = 1i * caustic;   % apply the caustic phase change
        end
    end
    
    % *** find indices of bracketted receiver ranges ***
    % following assumes r is uniformly spaced
    deltar = ( Pos.r.range( end ) - Pos.r.range( 1 ) ) / ( length( Pos.r.range ) - 1 );
    
    if ( ray( is + 1 ).x( 1 ) >= ray( is     ).x( 1 ) )   % stepping to the right
        ir1 = ceil(  1 + ( ray( is     ).x( 1 ) - Pos.r.range( 1 ) ) / deltar + eps( ray( is ).x( 1 ) ) );
        ir2 = floor( 1 + ( ray( is + 1 ).x( 1 ) - Pos.r.range( 1 ) ) / deltar );
        ir1 = max( ir1, 1 );
        ir2 = min( ir2, length( Pos.r.range ) );
    else                                                  % stepping to the left
        ir1 = ceil(  1 + ( ray( is + 1 ).x( 1 ) - Pos.r.range( 1 ) ) / deltar + eps( ray( is ).x( 1 ) ) );
        ir2 = floor( 1 + ( ray( is     ).x( 1 ) - Pos.r.range( 1 ) ) / deltar );
        ir1 = max( ir1, 1 );
        ir2 = min( ir2, length( Pos.r.range ) );
    end
    
    % compute ray tangent and normal
    if ( ir2 >= ir1 )
        xray = ones( size( Pos.r.depth ) )' * ray( is ).x;
        %xray = ray( is ).x;
        %xray = xray( ones( length( Pos.r.depth ), 1 ), : );

        tray = ray( is + 1 ).x - ray( is ).x;
        rlen = norm( tray );
        
        %IF ( rlen < TINY( xV( 1, I ) ) ) CYCLE   % if duplicate point in ray, skip to next step along the ray
        tray        = tray / rlen;
        tray_scaled = tray / rlen;
        
        nray = [ -tray( 2 ) tray( 1 ) ];      % unit normal  to ray
    end

    for ir = ir1 : ir2   % loop across receiver ranges
        
        %xrcvr = [ Pos.r.range( ir ) * ones( size( Pos.r.depth ) )' Pos.r.depth' ];
        xrcvr = Pos.r.range( ir );
        xrcvr = [ xrcvr( ones( length( Pos.r.depth ), 1 ) ) Pos.r.depth' ];

        s     =      ( xrcvr - xray ) * tray_scaled';  % proportional distance along ray
        n     = abs( ( xrcvr - xray ) * nray' );       % normal distance to ray
        
        q     = ray( is ).q( 1 ) + s * ( ray( is+1 ).q( 1 ) -  ray( is ).q( 1 ) );  % interpolated amplitude
        
        sigma = abs( q / q0 );          	            % beam radius
        %ii =  sigma < min( 0.2 * is * deltas / lambda , pi * lambda ) ;
        %sigma( ii ) = min( 0.2 * is * deltas / lambda , pi * lambda );
        ii =  sigma < min( 0.2 * freq * ray( is ).tau / lambda , pi * lambda ) ;
        sigma( ii ) = min( 0.2 * freq * ray( is ).tau / lambda , pi * lambda );
        
        ird = find( n < IBWin * sigma );   % take subset of beams that contributes
        
        if ( ~isempty( ird ) )
            A        = abs( q0 ./ q( ird ) );
            delay    = ray( is ).tau + s( ird ) * ( ray( is + 1 ).tau - ray( is ).tau ); % interpolated delay
            
            % shift phase for rays that have passed through a caustic
            caust = caustic * ones( 1, length( ird ) );
            ii    = ( ( q( ird ) <= 0.0 ) & ( ray( is ).q( 1 ) > 0.0 ) ) | ...
                    ( ( q( ird ) >= 0.0 ) & ( ray( is ).q( 1 ) < 0.0 ) );
            caust( ii ) = 1i * caust( ii );   % apply the caustic phase change

            const = Rat1 * ray( is ).Rfa * realsqrt( ray( is ).c ./ abs( q( ird ) ) ) .* caust.';
            
            if ( RunTypeE == 'S' ) % Coherent, semi-coherent, or incoherent TL:
                const = DS .* const;
            end
            
            amp = const .* exp( -0.5 * ( n( ird ) ./ sigma( ird ) ).^2 ) ./ ( sigma( ird ) .* A );
            
            switch( RunTypeE )
                case ( 'E' )      % eigenrays
                    %WRTRAY( xv, Nsteps, DepthT, DEPTHB )
                    
                case ( 'A' )      % arrivals
                    angle  = 180 / pi * atan2( tray( 2 ), tray( 1 ) );
                    AddArr( omega, ird, ir, amp, delay, angle, MxNarr );
                    
                case ( 'C'  )     % coherent TL
                    contri       = amp .* exp( -1i * omega * delay );
                    U( ird, ir ) = U( ird, ir ) + contri;
                    
                otherwise         % incoherent/semi-coherent TL
                    W =  exp( -0.5 * ( n( ird ) ./ sigma( ird ) ).^2 ) ./ ( 2. * sigma( ird ) .* A );   % Gaussian decay
                    contri       = ( abs( amp ) .* exp( omega * imag( delay ) ) ./ W ).^2 .* W;
                    U( ird, ir ) = U( ird, ir ) + contri;
                    
            end  % close switch on RunTypeE
        end   % next receiver depth
    end   % next receiver range
end   % next step along the ray