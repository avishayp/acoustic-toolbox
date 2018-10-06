function C = field_noise( ShadeFile, sd, rd, freq )

% Computes a noise covariance matrix by summing columns of the pressure field
% ShadeFile is the name of the shade file
% sd is the vector of source depths
% rd is the vector of receiver depths
%
% Returns:
% C is the nrd x nrd covariance matrix
%
% If you want the noise level in dB use:
%    NL = 10 * log10( diag( C ) ) + rho_SL_dB;

% mike porter Nov. 2017

% read in the pressure array
[ PlotTitle, ~, freqVec, ~, Pos, pressure ] = read_shd( ShadeFile, freq );


%%
% interpolate pressure field to user provided sd, rd

if length( Pos.s.depth ) > 1
   Ptemp = interp1( Pos.s.depth, squeeze( pressure( 1, :, :, : ) ), sd );
   Ptemp = squeeze( Ptemp );
else
   zs   = sd;
   isd  = find( Pos.s.depth >= zs );    % index of source depth
   isd  = isd( 1 );
   Ptemp = squeeze( pressure( 1, isd, :, : ) );
end

P = interp1( Pos.r.depth, squeeze( Ptemp( :, : ) ), rd );

% loop over range, adding the contribution of each vector into the
% covariance matrix

nrd = length( rd );
C = zeros( nrd, nrd );

for ir = 1 : length( Pos.r.range ) - 1
   
    dr = Pos.r.range( ir + 1 ) - Pos.r.range( ir );
    d = P( :, ir ) * sqrt( 2 * pi * Pos.r.range( ir ) * dr );

    C = C + d * d';
end

