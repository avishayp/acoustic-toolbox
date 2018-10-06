function U = scalep( Dalpha, c, r, U, RunType, TopOpt, freq )

% Scale the pressure field

% Compute scale factor
switch ( RunType(2:2) )
    case ( 'C' )
        const = -Dalpha * sqrt( freq ) / c;
    case ( 'R' )
        const = -Dalpha * sqrt( freq ) / c;
    otherwise
        const = -1.0;
end

% Thorp attenuation?
if ( TopOpt(4:4) == 'T' )
    f2 = ( freq / 1000.0 )^2;
    
    % Original Thorp (1967) formula
%     alpha = 40.0 * f2 / ( 4100.0 + f2 ) + 0.1 * f2 / ( 1.0 + f2 );
%     alpha = alpha / 914.4;     % dB / m
%     alpha = alpha / 8.6858896; % Nepers / m
    
    % Updated formula from JKPS Eq. 1.34
    alpha = 3.3d-3 + 0.11 * f2 / ( 1.0 + f2 ) + 44.0 * f2 / ( 4100.0 + f2 ) + 3d-4* f2;   % dB/km
    alpha = alpha / 8685.8896; % Nepers / m
else
    alpha = 0.0;
end

% For incoherent RunType, convert intensity to pressure
if ( RunType(1:1) ~= 'C' )
    U = sqrt( real( U ) );
end

% add in attenuation
for ir = 1 : length( r )
    if RunType(4:4) == 'X' % line source
        factor = -4 * sqrt( pi ) * const;
    else
        if ( r( ir ) ~= 0 )
            factor = const * exp( -alpha * r( ir ) ) / sqrt( r( ir ) );
        else
            factor = 0.0;
        end
    end
    U( :, ir ) = factor * U( :, ir );
end