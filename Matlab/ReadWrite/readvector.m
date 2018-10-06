function [ x, Nx ] = readvector( fid )

% Reads a row-vector
%
% This routine emulates a Fortran language capability that allows a '/'
% to be used to terminate input
%
% A user can type:
%   5     ! Number of values
%   0 1000 /
% as a shorthand for
%   5
%   0 250 500 750 1000
% Here the '/' terminates the reading and the code then creates 4 equally
% spaced values between 0 and 1000.
% However, one also has the option of using specific (non-equispaced
% values)

Nx = fscanf( fid, '%i', 1 );
fgetl( fid );
tline = fgetl( fid );

if ( contains( tline, '/' ) )  % does it contain '/'?
   x = sscanf( tline, '%f', 2 )';
   
   if Nx > 2
      x = linspace( x( 1 ), x( 2 ), Nx ); % generate the vector
   end
else
   x = sscanf( tline, '%f', Nx )'; % read the vector (transpose to return a row vector)
end
