function make_bdry( interp_type )

% generate Elliptic bottom bathymetry
% Based on Diana McCammon's test case and her code
% mbp Feb. 2009

a = 20177.67;
b = 10044.32;

h = 17500.;
r1 = 0.;
z1 = 5000.;

gam = 400.;   % chord length

rmax = 37600.;

np = 100;
eps = 0.005;
eps = 0.0000001;

z = np;
r = np;
chord = np;
r( 1 ) = r1;
z( 1 ) = z1;

for i = 2 : np
   % finding each point
   rn    = r ( i - 1 ) + gam;
   delta = 100.;
   
   while abs( delta ) > eps    %iterate for solution to combinational equation
      fofr  = ( b * sqrt( 1 - ( rn - h )^2 / a^2 ) - z( i - 1 ) )^2 + ( rn - r( i - 1 ) )^2 - gam^2;
      fprim = b / 2 * 1. / sqrt( 1. - ( rn - h )^2 / a^2 ) + rn;
      delta = -fofr / fprim;
      rn    = rn + delta;
   end
   
   r( i )         = rn;
   z( i )         = b * sqrt( 1. - ( r( i )- h )^2 / a^2 );
   chord( i - 1 ) = sqrt( ( r( i ) - r( i - 1 ) ) ^2 + ( z( i ) - z( i - 1 ) )^2 );
   
   % jump out if we get a range outside rmax
   if r( i ) > rmax
      break
   end
end

% write the alitmetry file
fid = fopen( 'Ellipse.ati', 'w' );
fprintf( fid, '''%c'' \n', interp_type );
fprintf( fid, '%i \n', length( r ) );
for ii = 1: length( r )
   fprintf( fid, ' %f %f \n', r( ii )/1000, -z( ii ) );
end
fclose( fid );

% write the bathymetry file
fid = fopen( 'Ellipse.bty', 'w' );
fprintf( fid, '''%c'' \n', interp_type );
fprintf( fid, '%i \n', length( r ) );
for ii = 1: length( r )
   fprintf( fid, ' %f %f \n', r( ii )/1000, z( ii ) );
end
fclose( fid );

