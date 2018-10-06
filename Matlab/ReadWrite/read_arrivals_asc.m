function [ Arr, Pos ] = read_arrivals_asc( ARRFile, Narrmx )

% Read the arrival time/amplitude data computed by BELLHOP
%
% usage:
%[ Arr, Pos ] = read_arrivals_asc( ARRFile, Narrmx );
%
% Arr is a structure containing all the arrivals information
% Pos is a structure containing the positions of source and receivers
%
% ARRFile is the name of the Arrivals File
% Narrmx is the maximum number of arrivals allowed
% mbp 9/96

if nargin == 1
   Narrmx = 200;
end

fid = fopen( ARRFile, 'r');	% open the file
if ( fid == -1 )
   error( 'read_arrivals_asc: Arrivals file cannot be opened' )
end

% read the header info

freq = fscanf( fid, '%f',  1  );

Nsd  = fscanf( fid, '%i',  1  );   % number of source   depths
Nrd  = fscanf( fid, '%i',  1  );   % number of receiver depths
Nrr  = fscanf( fid, '%i',  1  );   % number of receiver ranges

Pos.s.depth   = fscanf( fid, '%f', Nsd );   % source   depths
Pos.r.depth   = fscanf( fid, '%f', Nrd );   % receiver depths
Pos.r.range   = fscanf( fid, '%f', Nrr );   % receiver ranges

% loop to read all the arrival info (delay and amplitude)

Arr.A         = zeros( Nrr, Narrmx, Nrd, Nsd );
Arr.delay     = zeros( Nrr, Narrmx, Nrd, Nsd );
Arr.SrcAngle  = zeros( Nrr, Narrmx, Nrd, Nsd );
Arr.RcvrAngle = zeros( Nrr, Narrmx, Nrd, Nsd );
Arr.NumTopBnc = zeros( Nrr, Narrmx, Nrd, Nsd );
Arr.NumBotBnc = zeros( Nrr, Narrmx, Nrd, Nsd );

for isd = 1 : Nsd
   Narrmx2 = fscanf( fid, '%i', 1 );  % max. number of arrivals to follow
   disp( [ 'Max. number of arrivals for source index ', num2str( isd ), ' is ', num2str( Narrmx2 ) ] );
   for ird = 1:Nrd
      for ir = 1:Nrr
         Narr = fscanf( fid, '%i', 1 );	% number of arrivals
         Arr.Narr( ir, ird, isd ) = Narr;
         
         if Narr > 0   % do we have any arrivals?
            da = fscanf( fid, '%f', [ 8, Narr ] );
            Narr = min( Narr, Narrmx ); % we'll keep no more than Narrmx values
            Arr.Narr( ir, ird, isd ) = Narr;

            Arr.A(         ir, 1:Narr, ird, isd ) = da( 1, 1:Narr ) .* exp( 1i * da( 2, 1:Narr ) * pi/180);
            Arr.delay(     ir, 1:Narr, ird, isd ) = da( 3, 1:Narr ) + 1i * da( 4, 1:Narr );
            Arr.SrcAngle(  ir, 1:Narr, ird, isd ) = da( 5, 1:Narr );
            Arr.RcvrAngle( ir, 1:Narr, ird, isd ) = da( 6, 1:Narr );
            Arr.NumTopBnc( ir, 1:Narr, ird, isd ) = da( 7, 1:Narr );
            Arr.NumBotBnc( ir, 1:Narr, ird, isd ) = da( 8, 1:Narr );
         end
      end		% next receiver range
   end		% next receiver depth
end	% next source depth

fclose( fid );
