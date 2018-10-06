function Fire_pe( filename )

% PE
% (a simple MATLAB Parabolic equation model ...)
% M. Porter
% Uses a simple Crank-Nicholson scheme like the IFD PE
% See write up in JKPS
% psi is the pressure (not reduced pressure)
% sound speeds need to conjugated relative to other models in the Acoustics
% Toolbox
%
% First version written 7/94 at the Scripps Institution of Oceanography
% Extensively updated Feb. 2014

clear global
tic

%%
% Here are some additonal parameters you can set to control the PE run

% Think about the two-way travel through the bottom and the accumulated
% attenuation
% Typically 50 wavelengths is enough to kill of the return from the final
% bottom
Nlambda = 50; % number of wavelengths to extend the bottom depth to

% If you use fewer than 20 points/wavelength, you may see a bit of noise
Nsamples_per_wavelength = 10;   % finite difference sampling

c0      = 1500;   % reference sound speed
StarterType = 'Greene';

%%
if ( isempty( filename ) )
    warndlg( 'No envfil has been selected', 'Warning' );
end

global omega Bdry
global Pos
global xBot

% filenames
envfil = [ filename '.env' ];   % input environmental file
shdfil = [ filename '.shd.mat' ];   % output file name (pressure)
btyfil = [ filename '.bty' ];   % input bathymetry
atifil = [ filename '.ati' ];   % input altimetry
sspfil = [ filename '.ssp' ];   % input ssp file (if range-dependent SSP)

[ PlotTitle, freq, SSP, Bdry, Pos, Beam, ~, ~, fid ] = read_env( envfil, 'BELLHOP' );  % read the environmental file
fclose( fid );                    % close out the envfil

Pos.r.range = 1000.0 * Pos.r.range; % convert km to m

TopATI = Bdry.Top.Opt(5:5);
readati( atifil, TopATI, Bdry.Top.depth, Beam.Box.r );    % read in the altimetry  file

BotBTY = Bdry.Bot.Opt(2:2);
readbty( btyfil, BotBTY, Bdry.Bot.depth, Beam.Box.r );    % read in the bathymetry file

is = 1;    % index of source depth
zs = Pos.s.depth( is );			% source depth

omega  = 2 * pi * freq;
lambda = c0 / freq;
k0     = omega / c0;

d  = Bdry.Bot.depth;			% bottom depth
nz = fix( Nsamples_per_wavelength * d / lambda );    % number of finite difference points
h  = d / nz;	% mesh spacing
z  = linspace( 0, d, nz );      % grid coordinates

rmax = Pos.r.range( end );
nr   = fix( Nsamples_per_wavelength * rmax / lambda );   % number of finite difference points

Nrd  = length( Pos.r.depth );
Nrr  = length( Pos.r.range );
deltar = rmax / nr;

psit = zeros( Nrd, Nrr );   % sampled at receivers


%% optionally read a 2D SSP

c = interp1( SSP.z, conj( SSP.c ), z ).';
iProf = 1;   % counter to track the active profile at the current range step

if ( Bdry.Top.Opt( 1 : 1 ) == 'Q' )
    [ cmat, rProf, NProf, ~ ] = readssp2d( sspfil );
    rProf = 1000 * rProf;   % convert km to meters
    % !!! assuming rProf( 1 ) = 0 here; should search for first profile
    c = interp1( SSP.z, cmat( :, iProf ), z ).';
end

%% Form marching matrix

% generate ssp for absorbing bottom layer

d_HS   = d + Nlambda * lambda;   % make absorbing layer Nlambda wavelengths thick
nz_HS  = fix( Nsamples_per_wavelength * ( d_HS - d ) / lambda );   % number of finite difference points
z_HS   = linspace( d, d_HS, nz_HS ) + .01;   % f.d. grid for absorbing bottom layer
c_HS   = conj( Bdry.Bot.cp ) * ones( length( z_HS ), 1 );
rho_HS = Bdry.Bot.rho;

z_tot  = [ z z_HS ];
nz_tot = nz + nz_HS;

[ B, L, U ] = create_ABC( k0, c0, c, nz, d, c_HS, rho_HS, nz_HS, d_HS, deltar );


%% starter

switch ( StarterType )
    case( 'Gaussian' )
        % \psi(0,z) = \sqrt{ k_0 } \,
        % e^{ -{ k_0^2 \over 2 } ( z - z_s )^2 } (6.100)
        
        % Change fac to broaden the Gaussian, making it more narrow angle
        fac = 1;        % the usual formula has fac = 1
        
        psi = sqrt( k0 / fac ) * ...
            ( exp( -.5 * ( k0 / fac )^2 * ...
            ( ( z_tot - zs * ones( 1, nz_tot ) ) ).^2 )' ); % ... % direct
        %exp( -( k0 / fac )^2 * ...
        %( ( z_tot + zs * ones( 1, nz_tot ) ) ).^2 )' );    % surface image
        
    case( 'Greene' )
        % Greene's source JKPS(6.102)
        psi = sqrt( k0 ) * ( 1.4467 - 0.4201 * k0^2 * ( z_tot.' - zs * ones( nz_tot, 1 ) ).^2 ) .* ...
            ( exp( -( k0 )^2 * ( 1 / 3.0512 ) * ...
            ( ( z_tot.' - zs * ones( nz_tot, 1 ) ) ).^2 ) ); % ... % direct
        
    case( 'delta' )
        % delta function starter:
        % psi = zeros( 1, nz_tot );
        %
        % [ ~, isd ] = min( z_tot - zs );
        % b   = zeros( 1, nz_tot );
        % b( isd ) = 1.;
    otherwise
        error( 'Unknown option for starting field' )
end

%% March out in range

new_matrix = 0;   % flag to indicate whether the f.d. matrix needs updating
ircvr = 1;

for ir = 1 : nr - 1
    
    range_march = ir * deltar;
    if ( mod( ir, fix( nr / 10 ) ) == 0 )
        fprintf( 'Range = %9.5g km \n', range_march / 1000 )
    end
    
    % if range dependent bottom, update the matrices
    if ( BotBTY == '~' )
        d_new = interp1( xBot( 1, : ), xBot( 2, : ), range_march );			% bottom depth
        d_HS  = d + Nlambda * lambda;   % make absorbing layer Nlambda wavelengths thick
        
        % show progress
        if ( mod( ir, fix( nr / 50 ) ) == 0 )
            fprintf( 'Range = %9.5g km   Depth = %9.5g m \n', range_march / 1000, d_new )
        end
        
        %% update the f.d. matrices if the depth changed
        
        if ( abs( d_new - d ) > lambda / Nsamples_per_wavelength )
        %if ( abs( d_new - d ) >= h )

            nz_new = fix( Nsamples_per_wavelength * d_new / lambda );   % # of finite difference points
            h  = d_new / nz_new;  	% mesh spacing
            z_new  = linspace( 0, d_new, nz_new );	% grid coordinates
            
            % adjust d_HS a bit so that the mesh size, h, divides it
            d_HS   = d_new + nz_HS * h;
            z_HS   = linspace( d_new, d_HS, nz_HS ) + .01;   % f.d. grid for absorbing bottom layer
            
            z_tot_new = [ z_new z_HS ];
            
            % interpolate the pressure field onto the new grid
            psi   = interp1( z_tot, psi, z_tot_new, 'linear', 0 ).';
            
            z     = z_new;
            z_tot = z_tot_new;
            d     = d_new;
            nz    = nz_new;
            
            c = interp1( SSP.z, SSP.c, z ).';
            new_matrix = 1;   % set flag to indicate we need to create a new marching matrix
        end
        
        %% update the f.d. matrix if the SSP changed or grid has changed
        
        if ( Bdry.Top.Opt( 1 : 1 ) == 'Q' )
            
            if ( ir * deltar > rProf( iProf + 1 ) || mod( ir, fix( nr / 50 ) ) == 0  || new_matrix )
                
                % !!! assuming here there is a profile to find
                ii = find( ir * deltar < rProf );   % first profile left of current range
                
                if ( numel( ii ) > 0 )
                    iProf = ii( 1 ) - 1;
                    iProf = min( max( iProf, 1 ), NProf - 1 );
                    
                    c1 = interp1( SSP.z, cmat( :, iProf     ), z ).';
                    c2 = interp1( SSP.z, cmat( :, iProf + 1 ), z ).';
                    
                    % proportional distance in range between two SSP profiles
                    alpha = ( ir * deltar - rProf( iProf ) ) / ( rProf( iProf + 1 ) - rProf( iProf ) );
                    
                    c     = ( 1 - alpha ) * c1 + alpha * c2;
                    new_matrix = 1;   % set flag to indicate we need to create a new marching matrix
                end
            end
        end
        
        % Form marching matrix
        if ( new_matrix )
            [ B, L, U ] = create_ABC( k0, c0, c, nz, d, c_HS, rho_HS, nz_HS, d_HS, deltar );
            new_matrix = 0;   % flag to indicate whether the f.d. matrix needs updating  
        end
    end
    
    % equivalent to C psi_new = B * psi;
    y1 = B * psi;
    y  = L \ y1;
    psi_new = U \ y;
    
    %% interpolate the field onto the receiver grid
    % if we bracket a receiver, calculate the field at the receiver
    
    while ir * deltar > Pos.r.range( ircvr )   % bracketted?
        
        deltar_temp = ir * deltar - Pos.r.range( ircvr );  % range step to receiver
        
        if ( deltar_temp > 0 )
            % probably could avoid doing a new LU decomposition here
            
            [ Btemp, Ltemp, Utemp ] = create_ABC( k0, c0, c, nz, d, c_HS, rho_HS, nz_HS, d_HS, deltar );
            
            % equivalent to C psi_new = B * psi;
            y1 = Btemp * psi;
            y  = Ltemp \ y1;
            psit( :, ircvr ) = interp1( z_tot, Utemp \ y, Pos.r.depth, 'linear', 0 );
            
            ircvr = ircvr + 1;
            if ( ircvr > Nrr )
                break
            end   % jump out if marched past last receiver
        end
    end   % check next receiver
    psi = psi_new;
end

toc

%% Scale and save the field

% put back Hankel function
rt = Pos.r.range';
% hank = sqrt( 2 / ( pi * k0 ) ) * exp( 1i * ( k0 * rt - pi / 4 ) ) ...
%    * diag( 1.0 ./ sqrt( rt ) );

% scaled Hankel function matching (6.70) in JKPS
hank = exp( 1i * ( k0 * rt - pi / 4 ) ) ...
    * diag( 1.0 ./ sqrt( rt ) );

is = 1;
pressure( 1, is, :, : ) = full( psit * spdiags( hank.', 0, Nrr, Nrr ) );
PlotTitle = ['Fire Orig -' PlotTitle ];
PlotType  = 'rectilin  ';
atten     = 0;

save( shdfil, 'PlotTitle', 'PlotType', 'freq', 'atten', 'Pos', 'pressure' )

end

function [ B, L, U ] = create_ABC( k0, c0, c, nz, d, c_HS, rho_HS, nz_HS, d_HS, deltar )

% sets up the matrices for the finite difference equations
% SPE: 2ik_0 { \pa \psi \over {\pa r} }
%          + { \pa^2 \psi \over {\pa z^2} }
%          + k_0^2 ( n^2 - 1 ) \psi = 0  (6.8)

% ocean
n = c0 ./ c;
h = d / nz;	% mesh spacing

rho1 = 1;
rho2 = rho_HS;

% Absorbing bottom layer
n_HS   = c0 ./ c_HS;
nz_tot = nz + nz_HS;

% standard PE:
% E  = spdiags( ones( nz_tot, 1 ), -1, nz_tot, nz_tot );  % sub diagonal
% D1 = -2 * speye( nz_tot );  %     diagonal
% D2 = spdiags( [ h^2 * k0^2 * ( n.^2    - ones( nz,    1 ) );
%                 h^2 * k0^2 * ( n_HS.^2 - ones( nz_HS, 1 ) ) ], 0, nz_tot, nz_tot );
%
% A = D1 + D2 + E + E';
% B = 2 * 1i * k0 * h^2 / deltar * speye( size( A ) ) - A / 2;
% C = 2 * 1i * k0 * h^2 / deltar * speye( size( A ) ) + A / 2;
%

% Claerbout wide-angle PE:
a0 = 1.00;
a1 = 0.75;
b0 = 1.00;
b1 = 0.25;

w1  = b0 + .5 * 1i * k0 * deltar * ( a0 - b0 );
w1s = b0 - .5 * 1i * k0 * deltar * ( a0 - b0 );   % called w1* in jkps
w2  = b1 + .5 * 1i * k0 * deltar * ( a1 - b1 );
w2s = b1 - .5 * 1i * k0 * deltar * ( a1 - b1 );   % called w2* in jkps

r1 = w1  / w2;
r2 = w1s / w2s;

n2_tot = [ n.^2; n_HS.^2 ];

u    = 2 * ( k0^2 * h^2 / 2 * r2 - 1 ) + k0^2 * h^2 * ( n2_tot - 1 );
uhat = 2 * ( k0^2 * h^2 / 2 * r1 - 1 ) + k0^2 * h^2 * ( n2_tot - 1 );
nu = ones( nz_tot, 1 );

% jump at water/bottom interface
u( nz )    = ( rho1 + rho2 ) / rho2 * ( k0^2 * h^2 / 2 * r2 - 1 ) + ...
    .5 * k0^2 * h^2 * ( ( n2_tot( nz ) - 1 ) + rho1 / rho2 * ( n2_tot( nz + 1 ) - 1 ) );
uhat( nz ) = ( rho1 + rho2 ) / rho2 * ( k0^2 * h^2 / 2 * r1 - 1 ) + ...
    .5 * k0^2 * h^2 * ( ( n2_tot( nz ) - 1 ) + rho1 / rho2 * ( n2_tot( nz + 1 ) - 1 ) );
nu( nz + 1 )   = rho1 / rho2;   % offset by 1 to conform to Matlab spdiags convention

B = spdiags( ones( nz_tot, 1 ), -1, nz_tot, nz_tot ) + ...  % sub   diagonal
    spdiags( nu,                +1, nz_tot, nz_tot ) + ...  % super diagonal
    spdiags( uhat,               0, nz_tot, nz_tot );

B = w2 / w2s * B;

C = spdiags( ones( nz_tot, 1 ), -1, nz_tot, nz_tot ) + ...  % sub   diagonal
    spdiags( nu,                +1, nz_tot, nz_tot ) + ...  % super diagonal
    spdiags( u,                  0, nz_tot, nz_tot );

[ L, U ] = lu( C );   % factor C

end